#!/usr/bin/env python3
"""
Compute PET regional statistics (CerebrA), MMI/NCC to template, and LEFT-JOIN with
Slicer-exported metadata into a single wide CSV.

- Reads PET NIfTI from a directory (stripped + registered to MNI grid)
- Resamples CerebrA labels to PET grid (nearest)
- For each PET: per-label median/mean with optional intensity floor -> NaN
- Template metrics vs. GM template (Mattes MI, NCC) in GM mask
- LEFT-JOIN with metadata CSV by SeriesInstanceUID (extracted from filename) or fallback keys
- Handles 4D/5D inputs by taking the first frame (axis -1 repeatedly)
- Optional multiprocessing

Example:
  python3 scripts/metrics/compute_pet_regions.py \
    --pet-dir ./outputs/nii/mni_affine \
    --atlas-seg ./atlases/mni_icbm152_CerebrA_tal_nlin_sym_09c_expanded.nii \
    --atlas-labels ./atlases/CerebrA_LabelDetails.csv \
    --template ./atlases/MNIGMS.nii \
    --metadata ./outputs/metadata/all_series_metadata.csv \
    --out ./outputs/tables/pet_regions_mmi_ncc.csv \
    --intensity-nan-below 5000 --gm-threshold 0.20 --mmi-bins 50 --sample-percent 0.20 --jobs 6
"""
import os, re, sys, argparse, glob, math
import numpy as np
import pandas as pd
import SimpleITK as sitk
from typing import Dict, List, Tuple, Optional
from functools import partial

# ---------------- Helpers ----------------
def _abspath(p: str) -> str:
    return os.path.abspath(os.path.expanduser(p))

def ensure_dir(p: str):
    d = os.path.dirname(p)
    if d:
        os.makedirs(d, exist_ok=True)

UID_RX = re.compile(r"(?:^|_)(1(?:\.\d+){3,})(?:_|\.|$)")

def extract_series_uid_from_name(path: str) -> Optional[str]:
    base = os.path.basename(path)
    m_all = list(UID_RX.finditer(base))
    if not m_all:
        return None
    # choose the longest UID match
    uid = max((m.group(1) for m in m_all), key=len)
    return uid

# CerebrA label CSV reader (robust)
def _find_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    m = {re.sub(r"[^a-z0-9]", "", c.lower()): c for c in df.columns}
    for cand in candidates:
        key = re.sub(r"[^a-z0-9]", "", cand.lower())
        if key in m:
            return m[key]
    return None

def read_cerebra_labels(csv_path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(csv_path, sep=None, engine="python", encoding="utf-8")
    except Exception:
        df = pd.read_csv(csv_path, sep=None, engine="python", encoding="latin1")
    # Generic
    id_col   = _find_col(df, ["Label", "Index", "Id", "Value"]) 
    name_col = _find_col(df, ["LabelName", "Name", "Region", "Structure", "Label_Name"]) 
    if id_col and name_col:
        out = df[[id_col, name_col]].copy()
        out.columns = ["Label", "LabelName"]
        out["Label"] = out["Label"].astype(int)
        out["LabelName"] = out["LabelName"].astype(str)
        return out
    # CerebrA special (Label Name + RH/LH)
    name_col = _find_col(df, ["Label Name", "LabelName", "Name"]) 
    rh_col   = _find_col(df, ["RH Label", "RH_Label", "Right Label", "Right_Label"]) 
    lh_col   = _find_col(df, ["LH Labels", "LH Label", "LH_Labels", "Left Label", "Left_Label"]) 
    if not (name_col and rh_col and lh_col):
        raise ValueError(f"{csv_path}: expected CerebrA columns, got: {list(df.columns)}")
    rows = []
    for _, r in df.iterrows():
        name = str(r[name_col]).strip().replace(" ", "_")
        rh = r[rh_col]; lh = r[lh_col]
        try:
            if pd.notna(rh): rows.append({"Label": int(rh), "LabelName": name + "_Right"})
        except Exception: pass
        try:
            if pd.notna(lh): rows.append({"Label": int(lh), "LabelName": name + "_Left"})
        except Exception: pass
    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("CerebrA labels parsing produced empty result")
    return out

# Template metrics (MMI/NCC)
class TemplateMetrics:
    def __init__(self, template_path: str, gm_thresh: float = 0.20, mmi_bins: int = 50, sample_pct: float = 0.20, seed: int = 42):
        self.fixed_img = sitk.ReadImage(template_path, sitk.sitkFloat32)
        self.fixed_mask = sitk.BinaryThreshold(self.fixed_img, lowerThreshold=gm_thresh, upperThreshold=1e9, insideValue=1, outsideValue=0)
        self.mmi_bins = int(mmi_bins)
        self.sample_pct = float(max(0.01, min(0.99, sample_pct)))
        self.seed = int(seed)
    def evaluate(self, moving_img: sitk.Image) -> Tuple[float, float]:
        reg = sitk.ImageRegistrationMethod()
        reg.SetInterpolator(sitk.sitkLinear)
        try:
            reg.SetMetricFixedMask(self.fixed_mask)
        except Exception:
            pass
        if hasattr(reg, "SetMetricSamplingStrategy"):
            try:
                reg.SetMetricSamplingStrategy(reg.RANDOM)
            except Exception:
                pass
            if hasattr(reg, "SetMetricSamplingPercentage"):
                reg.SetMetricSamplingPercentage(self.sample_pct)
            if hasattr(reg, "SetMetricSamplingSeed"):
                reg.SetMetricSamplingSeed(self.seed)
        reg.SetInitialTransform(sitk.Transform(3, sitk.sitkIdentity), inPlace=False)
        reg.SetMetricAsMattesMutualInformation(numberOfHistogramBins=self.mmi_bins)
        mmi = float(reg.MetricEvaluate(self.fixed_img, moving_img))
        reg.SetMetricAsCorrelation()
        ncc = float(reg.MetricEvaluate(self.fixed_img, moving_img))
        return mmi, ncc

# IO helpers

def read_image(path: str) -> sitk.Image:
    return sitk.ReadImage(path, sitk.sitkFloat32)


def to_3d_first_frame(img: sitk.Image) -> sitk.Image:
    # convert SITK -> numpy, squeeze extra dims by taking index 0 on trailing axes until 3D
    arr = sitk.GetArrayFromImage(img)
    while arr.ndim > 3:
        arr = np.take(arr, 0, axis=0)  # SimpleITK array is [z, y, x, (t), ...] (time is front in numpy view)
    arr = np.asarray(arr, dtype=np.float32)
    out = sitk.GetImageFromArray(arr)
    out.SetSpacing(img.GetSpacing())
    out.SetOrigin(img.GetOrigin())
    out.SetDirection(img.GetDirection())
    return out


def resample_labels_to_pet(label_img: sitk.Image, pet_img: sitk.Image) -> sitk.Image:
    return sitk.Resample(
        label_img, pet_img,
        sitk.Transform(3, sitk.sitkIdentity),
        sitk.sitkNearestNeighbor,
        0, sitk.sitkUInt16,
    )


def region_stats(pet_img: sitk.Image, label_img_in_pet: sitk.Image, labels_df: pd.DataFrame, stats=("median","mean"), intensity_nan_below: Optional[float]=None) -> Dict[str, float]:
    pet = sitk.GetArrayFromImage(pet_img).astype(np.float32)
    labs = sitk.GetArrayFromImage(label_img_in_pet)
    if intensity_nan_below is not None:
        pet[pet < float(intensity_nan_below)] = np.nan
    out: Dict[str, float] = {}
    for _, r in labels_df.iterrows():
        lab = int(r["Label"]); name = str(r["LabelName"]) 
        mask = (labs == lab)
        if not mask.any():
            med = np.nan; meanv = np.nan
        else:
            vals = pet[mask]
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                med = np.nan; meanv = np.nan
            else:
                med = float(np.nanmedian(vals)) if "median" in stats else None
                meanv = float(np.nanmean(vals)) if "mean"   in stats else None
        if "median" in stats: out[f"Median_{name}"] = med
        if "mean"   in stats: out[f"Mean_{name}"]   = meanv
    return out

# --------------- Worker ---------------

def process_one(pet_path: str, atlas_labels_df: pd.DataFrame, atlas_seg_img: sitk.Image, tmpl: TemplateMetrics, intensity_nan_below: Optional[float], what_stats: Tuple[str,...]) -> Dict[str, object]:
    pet_img = read_image(pet_path)
    if sitk.GetArrayFromImage(pet_img).ndim > 3:
        pet_img = to_3d_first_frame(pet_img)
    lab_in_pet = resample_labels_to_pet(atlas_seg_img, pet_img)
    stats_map = region_stats(pet_img, lab_in_pet, atlas_labels_df, what_stats, intensity_nan_below)
    mmi, ncc = tmpl.evaluate(pet_img)
    row: Dict[str, object] = {
        "PET_Path": pet_path,
        "SeriesInstanceUID": extract_series_uid_from_name(pet_path),
        "MMI_to_MNIGMS": float(mmi),
        "NCC_to_MNIGMS": float(ncc),
    }
    row.update(stats_map)
    return row

# --------------- Main ---------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pet-dir', required=True)
    ap.add_argument('--atlas-seg', required=True)
    ap.add_argument('--atlas-labels', required=True)
    ap.add_argument('--template', required=True)
    ap.add_argument('--metadata', required=True, help='CSV from export_series_metadata.py')
    ap.add_argument('--out', required=True)
    ap.add_argument('--intensity-nan-below', type=float, default=5000)
    ap.add_argument('--gm-threshold', type=float, default=0.20)
    ap.add_argument('--mmi-bins', type=int, default=50)
    ap.add_argument('--sample-percent', type=float, default=0.20)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--stats', default='median,mean')
    ap.add_argument('--jobs', type=int, default=1)

    args = ap.parse_args()

    pet_dir = _abspath(args.pet_dir)
    atlas_seg = _abspath(args.atlas_seg)
    atlas_labels = _abspath(args.atlas_labels)
    template = _abspath(args.template)
    metadata_csv = _abspath(args.metadata)
    out_csv = _abspath(args.out)

    ensure_dir(out_csv)

    # Load atlas + labels
    labels_df = read_cerebra_labels(atlas_labels)
    seg_img   = read_image(atlas_seg)

    # Template metrics helper
    tmpl = TemplateMetrics(template, gm_thresh=args.gm_threshold, mmi_bins=args.mmi_bins, sample_pct=args.sample_percent, seed=args.seed)

    # Enumerate PET files
    pets = sorted(glob.glob(os.path.join(pet_dir, '*.nii'))) + sorted(glob.glob(os.path.join(pet_dir, '*.nii.gz')))
    if not pets:
        print(f"[ERR] No NIfTI in {pet_dir}")
        sys.exit(1)

    what_stats = tuple([s.strip().lower() for s in args.stats.split(',') if s.strip()])

    # Process (optionally parallel)
    rows: List[Dict[str, object]] = []
    worker = partial(process_one, atlas_labels_df=labels_df, atlas_seg_img=seg_img, tmpl=tmpl, intensity_nan_below=args.intensity_nan_below, what_stats=what_stats)

    if args.jobs and args.jobs > 1:
        import multiprocessing as mp
        with mp.Pool(processes=args.jobs) as pool:
            for res in pool.imap_unordered(worker, pets):
                rows.append(res)
    else:
        for p in pets:
            rows.append(worker(p))

    df_metrics = pd.DataFrame(rows)

    # Read metadata and JOIN
    meta = pd.read_csv(metadata_csv)
    join_key = 'SeriesInstanceUID'
    if join_key not in meta.columns:
        # Fall back to best-effort join using filename stem if needed
        print('[WARN] metadata has no SeriesInstanceUID; join may be partial')
    full = df_metrics.merge(meta, on=join_key, how='left', suffixes=('', '_meta'))
    # ---- Derived columns to match old training inputs ----
    import numpy as np

    def _as_num(s):
        return pd.to_numeric(full.get(s, np.nan), errors='coerce')

    # logAcqDur_s = log(ActualFrameDuration_s)
    if 'ActualFrameDuration_s' in full.columns:
        full['logAcqDur_s'] = np.log(np.clip(_as_num('ActualFrameDuration_s'), 1e-12, None))

    # (volitelné) logVoxelVol = log(VoxelVolume_mm3)
    if 'VoxelVolume_mm3' in full.columns:
        full['logVoxelVol'] = np.log(np.clip(_as_num('VoxelVolume_mm3'), 1e-12, None))

    # binarizace flagů na 0/1 (stejně jako v původním skriptu)
    def _to_bin(col):
        if col not in full.columns: 
            return
        x = full[col].astype(str).str.lower().str.strip()
        full[col] = x.isin(['1','true','t','yes','y']).astype(float)

    for nm in [
        'HasTOF','HasPSF','HasQClear',
        'corr_ATTN','corr_SCAT','corr_NORM','corr_DECY','corr_DTIM','corr_RAN','corr_DRFT','corr_GEO','corr_PVC'
    ]:
        _to_bin(nm)

    # Uložení
    full.to_csv(out_csv, index=False)
    print(f"[OK] wrote {out_csv} (rows={len(full)}, cols={len(full.columns)})")

#   full.to_csv(out_csv, index=False)
 #   print(f"[OK] wrote {out_csv} (rows={len(full)}, cols={len(full.columns)})")

if __name__ == '__main__':
    main()
