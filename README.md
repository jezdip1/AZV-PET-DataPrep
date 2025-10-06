# AZV‑PET Data Preparation Pipeline (beta)

End‑to‑end, reproducible preparation of PET volumes and region‑level features for downstream modelling/classification. This repository standardizes the scattered scripts we’ve been using into a single, documented pipeline.

---

## At a glance

**Inputs**
- A populated **3D Slicer** DICOM database (patients with PET/CT series)
- Atlas: *CerebrA* on MNI ICBM152 2009c symmetrical space
- Two shell helpers for stripping & registration (affine to MNI)

**Outputs**
- `outputs/metadata/all_series_metadata.csv` — series‑level DICOM metadata
- `outputs/nii/raw/*.nii` — exported volumes from Slicer
- `outputs/nii/stripped/*.nii` — skull‑stripped volumes
- `outputs/nii/mni_affine/*.nii` — volumes registered to MNI (affine)
- `outputs/tables/region_medians_*.csv` — region medians per volume (+ merged clinical columns)

---

## Repository layout (proposed)

```
AZV-PET-DataPrep/
├─ config/
│  └─ paths.example.yaml          # Edit & copy to paths.yaml
├─ scripts/
│  ├─ slicer/
│  │  ├─ export_all_metadata.py   # Export DICOM series metadata (CSV)
│  │  ├─ export_nii_from_csv.py   # Export NIfTI by PatientID/Series from CSV
│  │  └─ region_medians.py        # Region stats (medians) + merge clinical
│  └─ utils/
│     └─ check_env.py             # Simple checks (Slicer, paths, atlas)
├─ pipelines/
│  ├─ Makefile                    # One‑command orchestration
│  └─ run_all.sh                  # Bash orchestration (if you prefer)
├─ atlases/
│  ├─ CerebrA_LabelDetails.csv
│  └─ mni_icbm152_CerebrA_tal_nlin_sym_09c_expanded.nii
├─ tools/
│  ├─ stripping_batch             # Your existing skull‑strip batch
│  └─ register_affine_to_MNI      # Your existing affine registration batch
├─ outputs/
│  ├─ metadata/
│  ├─ nii/
│  │  ├─ raw/
│  │  ├─ stripped/
│  │  └─ mni_affine/
│  └─ tables/
└─ README.md
```

> **Why this split?**
> - All Slicer‑driven steps live under `scripts/slicer/` (run with Slicer’s Python).
> - File‑system side‑effects (`outputs/*`) are cleanly separated.
> - `config/paths.yaml` collects all paths so nothing stays hard‑coded in scripts.

---

## Requirements

- **3D Slicer 5.8+** with a populated local DICOM database
- Python inside Slicer (we run scripts with `Slicer --no-splash --python-script ...`)
- Atlas files: CerebrA labels CSV and segmentation NIfTI (see `atlases/`)
- Your batch helpers for **stripping** and **affine registration** available in `tools/`

Optional: a regular Python 3.10+ with `pandas` for quick CSV inspection (not required by the pipeline if you run everything inside Slicer).

---

## Quick start — one command

If your `config/paths.yaml` is set, the whole pipeline can be kicked off via Make:

```bash
make -C pipelines all
```

This runs: 1) metadata export → 2) NIfTI export → 3) stripping → 4) affine registration → 5) region medians (+merge) → final tables in `outputs/tables/`.

---

## Configuration

Create `config/paths.yaml` (copy from `paths.example.yaml`) and set your absolute paths:

```yaml
# config/paths.yaml
slicer:
  dicom_db: /path/to/SlicerDICOMDatabase
  exec: /path/to/Slicer-5.8.0-linux-amd64/Slicer

atlases:
  labels_csv: ./atlases/CerebrA_LabelDetails.csv
  seg_nifti:  ./atlases/mni_icbm152_CerebrA_tal_nlin_sym_09c_expanded.nii

io:
  metadata_csv: ./outputs/metadata/all_series_metadata.csv
  nii_raw_dir:  ./outputs/nii/raw
  nii_strip_dir: ./outputs/nii/stripped
  nii_mni_dir:  ./outputs/nii/mni_affine
  region_csv:   ./outputs/tables/region_medians_mni_affine.csv

options:
  threshold: 5000     # voxels below this set to NaN before computing medians
  pediatric_id_has_slash: false  # set true if PatientID uses Czech RC with “/”
```

---

## Step 1 — Export series metadata (Slicer)

> Output: `outputs/metadata/all_series_metadata.csv`

Run from your shell (Slicer CLI):

```bash
$SLICER_EXEC --no-splash --python-script scripts/slicer/export_all_metadata.py \
  -- --output ./outputs/metadata/all_series_metadata.csv
```

**Notes**
- The script walks **every series** in the Slicer DICOM DB and writes a row with patient, study/series descriptors, geometry, PET specifics (decay correction, energy windows, total dose, half‑life, uptake time), corrected‑image flags, and heuristic recon features (PSF/TOF/iterations/subsets/filter/Q.Clear), plus basic CT‑AC fields.
- If you need only PET emission series, you can filter by `Modality` and `ImageType` downstream or add a switch.

---

## Step 2 — Export NIfTI volumes (Slicer)

> Output: `outputs/nii/raw/*.nii`

Use the CSV from **Step 1** to drive export by `PatientID` and `Series`.

```bash
$SLICER_EXEC --no-splash --python-script scripts/slicer/export_nii_from_csv.py \
  -- --metadata ./outputs/metadata/all_series_metadata.csv \
     --outdir ./outputs/nii/raw \
     --patient-id-col PatientID \
     --uid-col StudyInstanceUID \
     [--pediatric-id-has-slash]
```

**What it does**
- For each row, calls `DICOMUtils.loadPatientByPatientID(...)` to load the series into the MRML scene and saves **scalar volume nodes** to `.nii`.
- Filenames follow: `UNIS_YYYYMMDD_SeriesNumber_Modality_SeriesDescription.nii` (unsafe characters sanitized). This aligns with the *Prefix* used later for merging clinical columns.
- Pediatric variant: if your PatientID is a Czech birth number, enable `--pediatric-id-has-slash` to automatically normalize IDs with a slash (RC `######/####`).

> Tip: You can filter which series to export (e.g., PET only, particular recon) using selectors in the script (SeriesDescription regex, ImageType has “EMISSION”, etc.).

---

## Step 3 — Skull stripping (shell)

> Output: `outputs/nii/stripped/*.nii`

Run your existing batch (edit paths in the script once):

```bash
./tools/stripping_batch ./outputs/nii/raw ./outputs/nii/stripped
```

Ensure the result preserves geometry and NIfTI headers.

---

## Step 4 — Affine registration to MNI (shell)

> Output: `outputs/nii/mni_affine/*.nii`

Run your existing affine registration batch (ANTS/FLIRT/Slicer transforms — pick one and stick to it):

```bash
./tools/register_affine_to_MNI ./outputs/nii/stripped ./outputs/nii/mni_affine
```

Recommendation: store transforms alongside outputs for traceability (e.g., `*.mat`/`*.h5`).

---

## Step 5 — Region medians (+ merge clinical) (Slicer)

> Output: `outputs/tables/region_medians_mni_affine.csv` and `_merged.csv`

```bash
$SLICER_EXEC --no-splash --python-script scripts/slicer/region_medians.py \
  -- --petdir ./outputs/nii/mni_affine \
     --labels-csv ./atlases/CerebrA_LabelDetails.csv \
     --seg-nifti  ./atlases/mni_icbm152_CerebrA_tal_nlin_sym_09c_expanded.nii \
     --metadata   ./outputs/metadata/all_series_metadata.csv \
     --threshold  5000 \
     --out ./outputs/tables/region_medians_mni_affine.csv
```

What happens:
- Loads the segmentation NIfTI (CerebrA expanded). For every PET volume in `--petdir`, computes **per‑hemisphere medians** for all labeled regions, after replacing values **below** `--threshold` by `NaN` (so masked background won’t skew medians).
- Writes the raw region table.
- Derives `Prefix = UNIS_YYYYMMDD` from the NIfTI filename and merges with the metadata CSV to attach clinical/study columns. Saves `*_merged.csv`.

> Swap to **means** by changing `np.nanmedian` to `np.nanmean` if needed. Medians are more robust to residual outliers.

---

## Orchestration via Makefile

`pipelines/Makefile` example:

```makefile
SLICER?=/path/to/Slicer-5.8.0-linux-amd64/Slicer

all: meta nii strip reg stats

meta:
	$(SLICER) --no-splash --python-script scripts/slicer/export_all_metadata.py -- --output outputs/metadata/all_series_metadata.csv

nii:
	$(SLICER) --no-splash --python-script scripts/slicer/export_nii_from_csv.py -- --metadata outputs/metadata/all_series_metadata.csv --outdir outputs/nii/raw

strip:
	./tools/stripping_batch outputs/nii/raw outputs/nii/stripped

reg:
	./tools/register_affine_to_MNI outputs/nii/stripped outputs/nii/mni_affine

stats:
	$(SLICER) --no-splash --python-script scripts/slicer/region_medians.py -- --petdir outputs/nii/mni_affine --labels-csv atlases/CerebrA_LabelDetails.csv --seg-nifti atlases/mni_icbm152_CerebrA_tal_nlin_sym_09c_expanded.nii --metadata outputs/metadata/all_series_metadata.csv --threshold 5000 --out outputs/tables/region_medians_mni_affine.csv
```

Run with:

```bash
make -C pipelines -j1
```

---

## Repro & naming guarantees

- **Filenames**: `UNIS_YYYYMMDD_*` everywhere → stable `Prefix` for joins.
- **Transforms**: keep `.mat`/`.h5` if you need to reapply or audit.
- **Logs**: redirect stdout/stderr from each step to `outputs/logs/*` if running on a server.

---

## Troubleshooting

- **No rows in metadata**: check Slicer points to the intended DICOM DB; ensure permissions.
- **Export fails for some patients**: missing series or PatientID mismatch. If pediatric IDs include `/`, enable the pediatric switch.
- **Registration off by a few mm**: verify orientation codes in NIfTI, and confirm the atlas is in the same MNI variant used by registration.
- **Merge mismatch**: `Prefix` requires consistent `UNIS` and **study date**. If your metadata’s date is in Excel serials, make sure it’s parsed with origin `1899‑12‑30`.

---

## Next steps (out of scope here)

- Compute **SUL_LOG** normalisation per region
- Derive **asymmetry indices** for paired regions
- Build modelling/classification bundles

---

## License & citation

Internal research use (AZV PET). Add an appropriate license if this becomes public.

