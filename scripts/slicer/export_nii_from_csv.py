#!/usr/bin/env python3
"""
Export selected NIfTI volumes from a curated metadata CSV.

- Zero hard-coded paths — everything via CLI args
- Works on a user-filtered CSV (keep only rows to export OR use --filter-col)
- Robust series selection:
    1) load patient by PatientID
    2) wait for Subject Hierarchy to settle
    3) try strict match by SeriesInstanceUID via SH (walk to parents)
    4) fallback match by (StudyDate, SeriesNumber, Modality, SeriesDescription) with normalization
    5) last resort: load series files from DB and re-check
- Clean headless use: exits Slicer with slicer.util.exit(...)
"""

import os, re, argparse, sys
import pandas as pd

# ---- helpers ----
def sanitize(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name)).strip("_._-")

def ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)

def boolish(v) -> bool:
    if v is None:
        return False
    s = str(v).strip().lower()
    return s in {"1","true","yes","y","t"}

# ---- Subject Hierarchy helpers ----
def sh_attr_upwards(node, key):
    """Return SH attribute value by walking from node's SH item up to parents until found."""
    import slicer
    sh = slicer.mrmlScene.GetSubjectHierarchyNode()
    item = sh.GetItemByDataNode(node)
    while item:
        val = sh.GetItemAttribute(item, key)
        if val not in (None, ""):
            return val
        item = sh.GetItemParent(item)
    return None

# ---- export core ----
def export_series_list(df: pd.DataFrame, outdir: str,
                       pid_col: str, study_uid_col: str, series_uid_col: str,
                       name_fields: list[str], debug: bool=False) -> list[tuple[str,str]]:
    """Load and export series listed in df. Returns list of (seriesUID_or_PID, outpath)."""
    import slicer
    from DICOMLib import DICOMUtils

    db = slicer.dicomDatabase
    written = []

    def norm_spaces(s: str) -> str:
        return re.sub(r"\s+"," ", str(s or "")).strip().lower()

    def norm_date(s: str) -> str:
        return re.sub(r"[^0-9]","", str(s or ""))  # 2023-10-23 -> 20231023

    def norm_series_no(s) -> str:
        try:
            return str(int(str(s).strip()))
        except Exception:
            return str(s or "").strip()

    for i, row in df.iterrows():
        series_uid = str(row.get(series_uid_col, '') or '').strip()
        patient_id = str(row.get(pid_col, '') or '').strip()
        row_series_desc = str(row.get('SeriesDescription','') or '')
        row_series_no   = str(row.get('SeriesNumber','') or '')
        row_modality    = str(row.get('Modality','') or '')
        row_study_date  = str(row.get('StudyDate','') or '')

        # Derive PatientID if missing
        if not patient_id and series_uid:
            try:
                files = db.filesForSeries(series_uid)
                if files:
                    pid = db.fileValue(files[0], '0010,0020')
                    if pid:
                        patient_id = pid
            except Exception:
                pass

        if not patient_id:
            print(f"[WARN] Row {i}: missing PatientID and cannot derive for SeriesInstanceUID={series_uid}. Skipping.", file=sys.stderr)
            continue

        # Load patient
        try:
            print(f"[INFO] Row {i}: loading PatientID={patient_id}")
            _ = DICOMUtils.loadPatientByPatientID(patient_id) or []
        except Exception as e:
            print(f"[WARN] Row {i}: loadPatientByPatientID failed for {patient_id}: {e}", file=sys.stderr)
            slicer.mrmlScene.Clear(0)
            continue

        slicer.app.processEvents()  # let SH settle

        # Consider all volume nodes currently in scene (some loaders don't return IDs)
        vol_nodes = [n for n in slicer.util.getNodesByClass('vtkMRMLVolumeNode')]

        if debug:
            print(f"[DBG] Scene has {len(vol_nodes)} volume nodes after loading PID={patient_id}")
            for dn in vol_nodes[:10]:
                print("      -",
                      sh_attr_upwards(dn, 'DICOM.SeriesInstanceUID'),
                      sh_attr_upwards(dn, 'DICOM.StudyDate'),
                      sh_attr_upwards(dn, 'DICOM.SeriesNumber'),
                      sh_attr_upwards(dn, 'DICOM.Modality'),
                      sh_attr_upwards(dn, 'DICOM.SeriesDescription'))

        # Pass 1: strict UID match (via parent-walk)
        nodes_to_save = []
        if series_uid:
            for node in vol_nodes:
                uid_attr = sh_attr_upwards(node, 'DICOM.SeriesInstanceUID')
                if uid_attr == series_uid:
                    nodes_to_save.append(node)

        # Pass 2: relaxed metadata match
        if not nodes_to_save:
            ns_date = norm_date(row_study_date)
            ns_no   = norm_series_no(row_series_no)
            ns_mod  = norm_spaces(row_modality)
            ns_desc = norm_spaces(row_series_desc)

            for node in vol_nodes:
                a_uid   = sh_attr_upwards(node, 'DICOM.SeriesInstanceUID') or ''
                a_date  = sh_attr_upwards(node, 'DICOM.StudyDate') or ''
                a_no    = sh_attr_upwards(node, 'DICOM.SeriesNumber') or ''
                a_mod   = sh_attr_upwards(node, 'DICOM.Modality') or ''
                a_desc  = sh_attr_upwards(node, 'DICOM.SeriesDescription') or ''

                if series_uid and a_uid and a_uid != series_uid:
                    continue
                if ns_date and norm_date(a_date) != ns_date:
                    continue
                if ns_no and norm_series_no(a_no) != ns_no:
                    continue
                if ns_mod and norm_spaces(a_mod) != ns_mod:
                    continue
                if ns_desc and norm_spaces(a_desc) != ns_desc:
                    continue
                nodes_to_save.append(node)

        # Pass 3: load a file from the series and re-check
        if series_uid and not nodes_to_save:
            files = db.filesForSeries(series_uid)
            if files:
                try:
                    print(f"[INFO] Row {i}: fallback – loading series from filesForSeries (n={len(files)})")
                    # Trigger DICOM plugin; even single file should create the series node(s)
                    slicer.util.loadVolume(files[0])
                    slicer.app.processEvents()
                    vol_nodes = [n for n in slicer.util.getNodesByClass('vtkMRMLVolumeNode')]
                    for node in vol_nodes:
                        uid_attr = sh_attr_upwards(node, 'DICOM.SeriesInstanceUID')
                        if uid_attr == series_uid:
                            nodes_to_save.append(node)
                except Exception as e:
                    print(f"[WARN] Row {i}: fallback load by files failed for {series_uid}: {e}", file=sys.stderr)

        if series_uid and not nodes_to_save:
            print(f"[WARN] Row {i}: no matching volume nodes for SeriesInstanceUID={series_uid}", file=sys.stderr)

        # filename
        parts = []
        for f in name_fields:
            v = row.get(f, '')
            if pd.isna(v):
                v = ''
            parts.append(sanitize(v))
        stem = "_".join([p for p in parts if p]) or sanitize(patient_id or series_uid or f"row{i}")

        # save
        for k, node in enumerate(nodes_to_save, start=1):
            suffix = f"_{k}" if len(nodes_to_save) > 1 else ""
            ofile = os.path.join(outdir, f"{stem}{suffix}.nii.gz")
            try:
                slicer.util.saveNode(node, ofile)
                print(f"[OK] Saved {ofile}")
                written.append((series_uid or f"PID:{patient_id}", ofile))
            except Exception as e:
                import traceback; traceback.print_exc()
                print(f"[ERR] save failed for {ofile}: {e}", file=sys.stderr)

        slicer.mrmlScene.Clear(0)

    return written

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--metadata', required=True, help='Curated metadata CSV (filtered by user)')
    ap.add_argument('--outdir',   required=True, help='Output folder for NIfTI files')

    # Column names present in CSV (defaults match export_all_metadata.py)
    ap.add_argument('--patient-id-col',   default='PatientID')
    ap.add_argument('--study-uid-col',    default='StudyInstanceUID')
    ap.add_argument('--series-uid-col',   default='SeriesInstanceUID')

    # Optional filtering on a column (e.g., Export==1)
    ap.add_argument('--filter-col',       default='')
    ap.add_argument('--filter-values',    default='', help='Comma-separated allowed values (case-insensitive). If empty, truthy values are used (1/yes/true).')

    # Filename template fields (ordered). Any of these should be column names in the CSV.
    ap.add_argument('--name-fields',      default='PatientID,StudyDate,SeriesInstanceUID,Modality,SeriesDescription')

    # Optional dry-run / debug
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--debug',   action='store_true')

    args, _ = ap.parse_known_args()

    # Load CSV
    df = pd.read_csv(args.metadata)

    # Apply filter if requested
    if args.filter_col:
        if args.filter_values:
            allow = {s.strip().lower() for s in args.filter_values.split(',') if s.strip()}
            df = df[df[args.filter_col].astype(str).str.lower().isin(allow)]
        else:
            df = df[df[args.filter_col].apply(boolish)]

    # Ensure outdir
    outdir = os.path.abspath(args.outdir)
    ensure_dir(outdir)

    # Name fields list
    name_fields = [s.strip() for s in args.name_fields.split(',') if s.strip()]

    if args.dry_run:
        print(f"[DRY] Would export {len(df)} rows → {outdir}")
        for i,row in df.head(10).iterrows():
            print("   ", row.get(args.series_uid_col, ''), row.get('SeriesDescription',''))
        import slicer; slicer.util.exit(0); return

    # Do the export
    try:
        written = export_series_list(df, outdir,
                                     pid_col=args.patient_id_col,
                                     study_uid_col=args.study_uid_col,
                                     series_uid_col=args.series_uid_col,
                                     name_fields=name_fields,
                                     debug=args.debug)
        print(f"[OK] Written {len(written)} NIfTI files to {outdir}")
        import slicer; slicer.util.exit(0)
    except Exception as e:
        import traceback; traceback.print_exc()
        print(f"[ERR] export_nii_from_csv failed: {e}", file=sys.stderr)
        import slicer; slicer.util.exit(1)

if __name__ == '__main__':
    main()
