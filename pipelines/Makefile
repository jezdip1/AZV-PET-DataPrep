SLICER?=/home/jezdip1/Slicer-5.8.0-linux-amd64/Slicer
ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/..)

all: meta nii strip reg stats

meta:
	$(SLICER) --no-splash --no-main-window --python-script $(ROOT)/scripts/slicer/export_all_metadata.py -- --output $(ROOT)/outputs/metadata/all_series_metadata.csv

nii:
	$(SLICER) --no-splash --no-main-window --python-script $(ROOT)/scripts/slicer/export_nii_from_csv.py -- --metadata $(ROOT)/outputs/metadata/all_series_metadata.csv --outdir $(ROOT)/outputs/nii/raw

strip:
	$(ROOT)/tools/stripping_batch $(ROOT)/outputs/nii/raw $(ROOT)/outputs/nii/stripped --jobs 4

reg:
	@SLICER=$(SLICER) $(ROOT)/tools/register_affine_to_MNI \
		$(ROOT)/outputs/nii/stripped \
		$(ROOT)/outputs/nii/mni_affine \
		--fixed $(ROOT)/atlases/MNIGMS.nii \
		--jobs 8

stats:
	python3 $(ROOT)/scripts/metrics/compute_pet_regions.py \
	  --pet-dir $(ROOT)/outputs/nii/mni_affine \
	  --atlas-seg $(ROOT)/atlases/mni_icbm152_CerebrA_tal_nlin_sym_09c_expanded.nii \
	  --atlas-labels $(ROOT)/atlases/CerebrA_LabelDetails.csv \
	  --template $(ROOT)/atlases/MNIGMS.nii \
	  --metadata $(ROOT)/outputs/metadata/all_series_metadata.csv \
	  --out $(ROOT)/outputs/tables/pet_regions_mmi_ncc.csv \
	  --intensity-nan-below 5000 \
	  --gm-threshold 0.20 \
	  --mmi-bins 50 \
	  --sample-percent 0.20 \
	  --jobs 6
