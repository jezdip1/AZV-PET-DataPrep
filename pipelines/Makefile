SLICER?=/home/jezdip1/Slicer-5.8.0-linux-amd64/Slicer
ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/..)

all: meta nii strip reg stats

meta:
	$(SLICER) --no-splash --python-script $(ROOT)/scripts/slicer/export_all_metadata.py -- --output $(ROOT)/outputs/metadata/all_series_metadata.csv

nii:
	$(SLICER) --no-splash --python-script $(ROOT)/scripts/slicer/export_nii_from_csv.py -- --metadata $(ROOT)/outputs/metadata/all_series_metadata.csv --outdir $(ROOT)/outputs/nii/raw

strip:
	$(ROOT)/tools/stripping_batch $(ROOT)/outputs/nii/raw outputs/nii/stripped

reg:
	$(ROOT)/tools/register_affine_to_MNI $(ROOT)/outputs/nii/stripped outputs/nii/mni_affine

stats:
	$(SLICER) --no-splash --python-script $(ROOT)/scripts/slicer/region_medians.py -- --petdir $(ROOT)/outputs/nii/mni_affine --labels-csv $(ROOT)/atlases/CerebrA_LabelDetails.csv --seg-nifti $(ROOT)/atlases/mni_icbm152_CerebrA_tal_nlin_sym_09c_expanded.nii --metadata $(ROOT)/outputs/metadata/all_series_metadata.csv --threshold 5000 --out $(ROOT)/outputs/tables/region_medians_mni_affine.csv
