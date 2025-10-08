#!/usr/bin/env python3
"""
Export all DICOM Study/Series metadata from the current 3D Slicer DICOM database
into a single CSV. This version RESPECTS the --output path passed from the CLI
(e.g., via Makefile) instead of saving to a hard-coded filename.

Usage (inside Slicer):
  Slicer --no-splash --python-script /path/to/export_series_metadata.py -- \
    --output /absolute/path/to/outputs/metadata/all_series_metadata.csv
"""
import re
import os
import argparse
import slicer
import pandas as pd
import pydicom
from pydicom.multival import MultiValue

# ---------- helpers ----------
def fmt_date(date_str):
    # YYYYMMDD → YYYY-MM-DD
    if date_str and isinstance(date_str, str) and len(date_str) == 8:
        return f"{date_str[0:4]}-{date_str[4:6]}-{date_str[6:8]}"
    return date_str or ""

def as_str(x):
    if x is None:
        return ""
    if isinstance(x, (list, tuple, MultiValue)):
        return " | ".join([str(y) for y in x])
    return str(x)

def as_float(x, default=None):
    try:
        if x is None or x == "":
            return default
        if isinstance(x, (list, tuple, MultiValue)):
            if len(x) == 0:
                return default
            return float(x[0])
        return float(x)
    except Exception:
        return default

def tm_to_seconds(tm):
    """DICOM TM (HHMMSS.frac) -> seconds since midnight."""
    if not tm:
        return None
    s = str(tm)
    try:
        hh = int(s[0:2]) if len(s) >= 2 else 0
        mm = int(s[2:4]) if len(s) >= 4 else 0
        ss = float(s[4:]) if len(s) > 4 else 0.0
        return hh * 3600 + mm * 60 + ss
    except Exception:
        return None

def minutes_diff(t_start, t_end):
    """Difference in minutes between two TM; if wrap past midnight, add 24h."""
    s1 = tm_to_seconds(t_start)
    s2 = tm_to_seconds(t_end)
    if s1 is None or s2 is None:
        return None
    d = s2 - s1
    if d < -12 * 3600:  # wrap over midnight
        d += 24 * 3600
    return d / 60.0

def parse_corrected_flags(ds):
    """0028,0051 CorrectedImage → one-hot flags."""
    flags = {"corr_ATTN": False, "corr_SCAT": False, "corr_NORM": False,
             "corr_DECY": False, "corr_DTIM": False, "corr_RAN": False,
             "corr_DRFT": False, "corr_GEO": False, "corr_PVC": False}
    ci = ds.get("CorrectedImage", None)
    if ci:
        toks = [t.strip().upper() for t in as_str(ci).replace("/", "\\").split("\\")]
        for k in list(flags.keys()):
            code = k.split("_", 1)[1]  # ATTN, SCAT, ...
            flags[k] = code in toks
    return flags

def parse_energy_windows(ds):
    """0054,0011 EnergyWindowInformationSequence → lower/upper keV (first window)."""
    low = None
    up  = None
    try:
        ewis = ds.get((0x0054, 0x0011), None)  # EnergyWindowInformationSequence
        if ewis and len(ewis) > 0:
            ewr = ewis[0].get((0x0054, 0x0012), None)  # EnergyWindowRangeSequence
            if ewr and len(ewr) > 0:
                low = as_float(ewr[0].get((0x0054, 0x0013), None))  # Lower
                up  = as_float(ewr[0].get((0x0054, 0x0014), None))  # Upper
    except Exception:
        pass
    return low, up

def parse_recon_text_fields(ds):
    """Concatenate text from several fields (Series/Protocol/Comments/Software/Private)."""
    parts = []
    for key in ["SeriesDescription", "ProtocolName", "ImageComments",
                "Manufacturer", "ManufacturerModelName", "SoftwareVersions",
                "ReconstructionMethod", "ConvolutionKernel", "ReconstructionAlgorithm"]:
        v = ds.get(key, None)
        if v:
            parts.append(as_str(v))
    # fallback: scan ds for private text with recon keywords
    hay = []
    try:
        for name in ds.dir():
            if re.search(r"(recon|iter|subset|filter|fwhm|tof|psf|q\.?clear|beta)", name, re.I):
                hay.append(f"{name}={as_str(getattr(ds, name))}")
    except Exception:
        pass
    if hay:
        parts.extend(hay)
    txt = " | ".join(parts)
    return txt

def parse_recon_features(txt):
    """Heuristic parsing of PSF/TOF, iterations, subsets, filters/FWHM, Q.Clear beta."""
    if not txt:
        txt = ""
    low = txt.lower()
    has_tof  = bool(re.search(r"\btof\b|time[-\s]?of[-\s]?flight|\btf\b", low))
    has_psf  = bool(re.search(r"\bpsf\b|truex|resolution\s*recovery", low))
    has_qclr = bool(re.search(r"q\.?clear", low))
    it = None
    m = re.search(r"(\d+)\s*(?:i\b|iter|iterations)", low)
    if m: it = int(m.group(1))
    sb = None
    m = re.search(r"(\d+)\s*(?:s\b|subset|subsets)", low)
    if m: sb = int(m.group(1))
    filter_type = None
    if re.search(r"gauss|gaussian", low): filter_type = "Gaussian"
    elif re.search(r"butterworth", low):  filter_type = "Butterworth"
    elif re.search(r"hann|hanning", low): filter_type = "Hann"
    fwhm = None
    m = re.search(r"fwhm[^0-9]*?(\d+(\.\d+)?)", low)
    if not m:
        m = re.search(r"gauss[^0-9]*?(\d+(\.\d+)?)", low)
    if m:
        fwhm = float(m.group(1))
    qbeta = None
    m = re.search(r"q\.?clear[^0-9]*?(\d+(\.\d+)?)", low)
    if m:
        qbeta = float(m.group(1))
    return has_psf, has_tof, has_qclr, it, sb, filter_type, fwhm, qbeta

# ---------- core export ----------
def export_patient_study_series_metadata(output_csv_path):
    db = slicer.dicomDatabase
    rows = []

    for patientUID in db.patients():
        for studyUID in db.studiesForPatient(patientUID):
            for seriesUID in db.seriesForStudy(studyUID):
                files = db.filesForSeries(seriesUID)
                if not files:
                    continue
                fname = files[0]
                ds = pydicom.dcmread(fname, stop_before_pixels=True)

                # Patient basics
                pn = ds.get("PatientName", None)
                family = getattr(pn, "family_name", "") if pn else ""
                given  = getattr(pn, "given_name", "")  if pn else ""

                patient_id    = ds.get("PatientID", "")
                birth_date    = fmt_date(ds.get("PatientBirthDate", ""))
                sex           = ds.get("PatientSex", "")

                # Study-level
                study_desc    = ds.get("StudyDescription", "")
                study_date    = fmt_date(ds.get("StudyDate", ""))
                study_time    = ds.get("StudyTime", "")
                study_id      = ds.get("StudyID", "")
                accession     = ds.get("AccessionNumber", "")
                referring     = ds.get("ReferringPhysicianName", "")

                # Series-level
                modality      = ds.get("Modality", "")
                series_desc   = ds.get("SeriesDescription", "")
                series_date   = fmt_date(ds.get("SeriesDate", "") or ds.get("StudyDate", ""))
                series_time   = ds.get("SeriesTime", "")
                acq_time      = ds.get("AcquisitionTime", "")
                station_name  = ds.get("StationName", "")
                manuf         = ds.get("Manufacturer", "")
                model_name    = ds.get("ManufacturerModelName", "")
                sw_versions   = ds.get("SoftwareVersions", "")
                image_type    = as_str(ds.get("ImageType", ""))

                # Geometry
                rows_i        = as_float(ds.get("Rows", None))
                cols_i        = as_float(ds.get("Columns", None))
                pxs           = ds.get("PixelSpacing", None)
                px_x          = as_float(pxs[0], None) if pxs is not None else None
                px_y          = as_float(pxs[1], None) if pxs is not None else None
                sl_thick      = as_float(ds.get("SliceThickness", None))
                sp_btw        = as_float(ds.get("SpacingBetweenSlices", None))
                recon_diam    = as_float(ds.get("ReconstructionDiameter", None))
                voxel_vol     = None
                if px_x and px_y and sl_thick:
                    voxel_vol = px_x * px_y * sl_thick  # mm^3

                # PET specifics
                units         = ds.get("Units", "")
                decay_corr    = ds.get("DecayCorrection", "")
                act_frame_dur = as_float(ds.get("ActualFrameDuration", None))  # ms
                frame_ref_ms  = as_float(ds.get("FrameReferenceTime", None))   # ms
                ewin_low, ewin_up = parse_energy_windows(ds)

                # Patient size/weight
                patient_size  = as_float(ds.get("PatientSize", None))
                patient_weight= as_float(ds.get("PatientWeight", None))

                # Radiopharm seq → dose / half-life / start time
                radiopharm            = ""
                radiopharm_start_time = ""
                total_dose            = None
                half_life             = None
                if hasattr(ds, "RadiopharmaceuticalInformationSequence"):
                    info = ds.RadiopharmaceuticalInformationSequence[0]
                    radiopharm            = as_str(getattr(info, "Radiopharmaceutical", ""))
                    radiopharm_start_time = as_str(getattr(info, "RadiopharmaceuticalStartTime", ""))
                    total_dose            = as_float(getattr(info, "RadionuclideTotalDose", None))
                    half_life             = as_float(getattr(info, "RadionuclideHalfLife", None))

                # Uptake time (min) – from injection start to SeriesTime (fallback AcquisitionTime)
                series_tm = series_time or acq_time
                uptake_min = minutes_diff(radiopharm_start_time, series_tm)

                # CorrectedImage → flags
                corr_flags = parse_corrected_flags(ds)

                # Recon text & features
                recon_text = parse_recon_text_fields(ds)
                has_psf, has_tof, has_qclr, iters, subsets, filt_type, fwhm, qbeta = parse_recon_features(recon_text)

                # CT-AC parameters (if CT series)
                kvp           = as_float(ds.get("KVP", None))
                conv_kernel   = as_str(ds.get("ConvolutionKernel", ""))

                # QA scaling
                rescale_slope     = as_float(ds.get("RescaleSlope", None))
                rescale_intercept = as_float(ds.get("RescaleIntercept", None))

                rows.append({
                    # patient
                    "PatientFamilyName":        family,
                    "PatientGivenName":         given,
                    "PatientID":                patient_id,
                    "PatientBirthDate":         birth_date,
                    "PatientSex":               sex,
                    "PatientSize_m":            patient_size,
                    "PatientWeight_kg":         patient_weight,

                    # study & series
                    "StudyDescription":         study_desc,
                    "StudyDate":                study_date,
                    "StudyTime":                study_time,
                    "StudyID":                  study_id,
                    "AccessionNumber":          accession,
                    "ReferringPhysician":       as_str(referring),

                    "SeriesDescription":        series_desc,
                    "SeriesDate":               series_date,
                    "SeriesTime":               series_time,
                    "AcquisitionTime":          acq_time,
                    "Modality":                 modality,
                    "ImageType":                image_type,

                    # device/software identity
                    "Manufacturer":             as_str(manuf),
                    "ManufacturerModelName":    as_str(model_name),
                    "SoftwareVersions":         as_str(sw_versions),
                    "StationName":              as_str(station_name),

                    # PET specifics
                    "Units":                    as_str(units),
                    "DecayCorrection":          as_str(decay_corr),
                    "ActualFrameDuration_s":    (act_frame_dur/1000.0 if act_frame_dur else None),
                    "FrameReferenceTime_ms":    frame_ref_ms,
                    "EnergyWindowLower_keV":    ewin_low,
                    "EnergyWindowUpper_keV":    ewin_up,

                    # Radiopharm
                    "Radiopharmaceutical":      radiopharm,
                    "RadiopharmStartTime":      radiopharm_start_time,
                    "RadionuclideTotalDose":    total_dose,
                    "RadionuclideHalfLife_s":   half_life,
                    "UptakeTime_min":           uptake_min,

                    # geometry
                    "Rows":                     rows_i,
                    "Columns":                  cols_i,
                    "PixelSpacingX_mm":         px_x,
                    "PixelSpacingY_mm":         px_y,
                    "SliceThickness_mm":        sl_thick,
                    "SpacingBetweenSlices_mm":  sp_btw,
                    "ReconstructionDiameter_mm":recon_diam,
                    "VoxelVolume_mm3":          voxel_vol,

                    # CorrectedImage flags
                    **corr_flags,

                    # recon heuristics
                    "ReconText":                recon_text,
                    "HasTOF":                   has_tof,
                    "HasPSF":                   has_psf,
                    "HasQClear":                has_qclr,
                    "Iterations":               iters,
                    "Subsets":                  subsets,
                    "FilterType":               filt_type,
                    "FilterFWHM_mm":            fwhm,
                    "QClearBeta":               qbeta,

                    # CT-AC (if Modality==CT, else empty)
                    "CT_KVP":                   kvp,
                    "CT_ConvolutionKernel":     conv_kernel,

                    # QA scale
                    "RescaleSlope":             rescale_slope,
                    "RescaleIntercept":         rescale_intercept,

                    # UIDs for joins
                    "StudyInstanceUID":         as_str(ds.get("StudyInstanceUID","")),
                    "SeriesInstanceUID":        as_str(ds.get("SeriesInstanceUID",""))
                })

    df = pd.DataFrame(rows)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', required=True, help='Output CSV path')
    args, _ = ap.parse_known_args()

    out_path = os.path.abspath(args.output)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    df = export_patient_study_series_metadata(out_path)
    df.to_csv(out_path, index=False, encoding='utf-8')
    print(f"Exportováno {len(df)} řádků do {out_path}")


if __name__ == "__main__":
    main()
    slicer.util.exit()