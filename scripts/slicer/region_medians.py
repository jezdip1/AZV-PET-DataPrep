# Run with Slicer: Slicer --no-splash --python-script region_medians.py -- --petdir DIR --labels-csv CSV --seg-nifti NII --metadata CSV --threshold 5000 --out CSV
import os, sys, argparse, csv
import numpy as np

def load_nifti(path):
    import SimpleITK as sitk
    img = sitk.ReadImage(path)
    arr = sitk.GetArrayFromImage(img)  # z,y,x
    return img, arr

def to_world_affine(img):
    import SimpleITK as sitk
    # Build an affine from SITK image meta
    spacing = np.array(list(img.GetSpacing()))[::-1]
    origin  = np.array(list(img.GetOrigin()))[::-1]
    direction = np.array(list(img.GetDirection())).reshape(3,3)[::-1, ::-1]
    A = np.eye(4)
    A[:3,:3] = direction @ np.diag(spacing)
    A[:3, 3] = origin
    return A

def nanmedian_by_label(data, labels, label_ids):
    out = {}
    for lid in label_ids:
        m = np.nanmedian(data[labels==lid]) if np.any(labels==lid) else np.nan
        out[int(lid)] = float(m) if np.isfinite(m) else np.nan
    return out

def run(args):
    import pandas as pd
    # Load label table
    labels_tbl = pd.read_csv(args.labels_csv)
    if 'Label' in labels_tbl.columns:
        id_col = 'Label'
    elif 'index' in labels_tbl.columns:
        id_col = 'index'
    else:
        id_col = labels_tbl.columns[0]
    label_ids = labels_tbl[id_col].astype(int).values

    # Load segmentation NIfTI
    seg_img, seg_arr = load_nifti(args.seg_nifti)
    # Enumerate PET files
    pet_files = [f for f in sorted(os.listdir(args.petdir)) if f.lower().endswith('.nii') or f.lower().endswith('.nii.gz')]

    rows = []
    for fname in pet_files:
        pet_path = os.path.join(args.petdir, fname)
        pet_img, pet_arr = load_nifti(pet_path)
        arr = pet_arr.astype(np.float32)
        # threshold to NaN
        arr[arr < args.threshold] = np.nan
        vals = nanmedian_by_label(arr, seg_arr, label_ids)
        rec = {'File': fname}
        for lid in label_ids:
            rec[f'Label_{int(lid)}_median'] = vals[int(lid)]
        rows.append(rec)
        print('[OK] stats for', fname)

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(args.out, index=False)
    print('[OK] wrote', args.out)

    # Try merge with metadata if available
    if args.metadata and os.path.exists(args.metadata):
        meta = pd.read_csv(args.metadata)
        # Derive Prefix UNIS_YYYYMMDD from filename start (split by first underscore * 2 parts*).
        # Adjust here if your file naming differs.
        def prefix_from_file(f):
            base = os.path.splitext(f)[0]
            parts = base.split('_')
            return '_'.join(parts[:2]) if len(parts)>=2 else base
        stats_df['Prefix'] = stats_df['File'].apply(prefix_from_file)
        if 'Prefix' in meta.columns:
            merged = stats_df.merge(meta, on='Prefix', how='left')
        else:
            merged = stats_df
        outm = os.path.splitext(args.out)[0] + '_merged.csv'
        merged.to_csv(outm, index=False)
        print('[OK] wrote', outm)

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--petdir', required=True)
    ap.add_argument('--labels-csv', required=True)
    ap.add_argument('--seg-nifti', required=True)
    ap.add_argument('--metadata', default='')
    ap.add_argument('--threshold', type=float, default=5000.0)
    ap.add_argument('--out', required=True)
    args, extra = ap.parse_known_args()
    run(args)
