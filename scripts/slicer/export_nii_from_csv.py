# Run with Slicer: Slicer --no-splash --python-script export_nii_from_csv.py -- --metadata CSV --outdir DIR [--patient-id-col PatientID] [--uid-col StudyInstanceUID] [--pediatric-id-has-slash]
import sys, os, argparse, re
try:
    import pandas as pd
except Exception as e:
    print('pandas is recommended for CSV parsing inside Slicer.', file=sys.stderr); raise

def norm_patient_id(pid, pediatric):
    pid = str(pid)
    return pid.replace('/','') if pediatric else pid

def sanitize(name):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(name))

def run(args):
    import DICOMLib
    from DICOMLib import DICOMUtils
    import slicer

    df = pd.read_csv(args.metadata)
    pid_col = args.patient_id_col
    uid_col = args.uid_col

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    failed = []

    for i,row in df.iterrows():
        pid = norm_patient_id(row[pid_col], args.pediatric_id_has_slash)
        try:
            with DICOMUtils.TemporaryDICOMDatabase() as db:
                # Use existing Slicer DB configured in settings; just load by PatientID
                ok = DICOMUtils.loadPatientByPatientID(pid)
                if not ok:
                    print(f'[WARN] Could not load PatientID={pid}', file=sys.stderr)
                    failed.append((pid, 'load'))
                    continue
                # save all scalar volumes in scene
                vols = slicer.util.getNodesByClass('vtkMRMLScalarVolumeNode')
                if not vols:
                    print(f'[WARN] No scalar volumes for PatientID={pid}', file=sys.stderr)
                    failed.append((pid, 'no_volumes')); continue
                for v in vols:
                    dn = v.GetName()
                    fname = sanitize(dn) + '.nii'
                    path = os.path.join(outdir, fname)
                    slicer.util.saveNode(v, path)
                    print('[OK] Saved', path)
                slicer.mrmlScene.Clear(0)
        except Exception as e:
            failed.append((pid, repr(e)))
            import traceback; traceback.print_exc()

    if failed:
        print('Failed entries:', failed, file=sys.stderr)

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--metadata', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--patient-id-col', default='PatientID')
    ap.add_argument('--uid-col', default='StudyInstanceUID')
    ap.add_argument('--pediatric-id-has-slash', action='store_true')
    args, extra = ap.parse_known_args()
    run(args)
