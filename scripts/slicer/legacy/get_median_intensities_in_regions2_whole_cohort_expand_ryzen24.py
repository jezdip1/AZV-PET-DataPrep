#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import slicer
from slicer.util import (
    loadSegmentation,
    arrayFromSegmentBinaryLabelmap,
    loadVolume,
    arrayFromVolume
)

# ——— Upravené cesty podle vaší struktury ———
atlas_csv        = '/media/jezdip1/UBURYZEN24_48T1/motol/AZV_PET/atlases/mni_icbm152_nlin_sym_09c_CerebrA_nifti/CerebrA_LabelDetails.csv'
seg_nifti        = '/media/jezdip1/UBURYZEN24_48T1/motol/AZV_PET/atlases/mni_icbm152_nlin_sym_09c_CerebrA_nifti/mni_icbm152_CerebrA_tal_nlin_sym_09c_expanded.nii'
pet_normy_csv    = '/media/jezdip1/UBURYZEN24_48T1/motol/AZV_PET/normy/all_series_metadata_export_12082025_pje.csv'
petdir           = '/media/jezdip1/UBURYZEN24_48T1/motol/AZV_PET/normy/whole_cohort_12082025/stripped/affine_reg_slicer/'
output_csv       = '/media/jezdip1/UBURYZEN24_48T1/motol/AZV_PET/normy/whole_cohort_12082025/stripped/affine_reg_slicer/region_medians_whole_cohort_affine_slicer_expanded_12082025.csv'

# Práh intenzity: hodnoty nad touto hranicí budou nastaveny na NaN
threshold_value  = 5000  # adjust as needed

# 1) Načíst atlasové popisky
# CSV obsahuje sloupce: Mindboggle ID, Label Name, RH Label, LH Labels
df_atlas = pd.read_csv(atlas_csv)
regions = []  # seznam (label_name, safe_name, rh_label, lh_label)
for _, row in df_atlas.iterrows():
    label_name = row['Label Name']
    try:
        rh = int(row['RH Label'])
        lh = int(row['LH Labels'])
    except (KeyError, ValueError):
        continue
    safe_name = label_name.strip().replace(' ', '_').replace('-', '_')
    regions.append((label_name, safe_name, rh, lh))
if not regions:
    raise RuntimeError('Nenalezeny platné regiony v atlasu.')

# 2) Načíst segmentation node pro maskování
segNode = loadSegmentation(seg_nifti)

# 3) Výpočet průměrných intenzit
results = []
for fname in sorted(os.listdir(petdir)):
    if not fname.lower().endswith('.nii'):
        continue
    print(f"Processing {fname}")
    vol_node = loadVolume(os.path.join(petdir, fname))
    # načíst data a převést na float pro podporu NaN
    vol_arr = arrayFromVolume(vol_node).astype(np.float32)
    # nastavit hodnoty pod prahem na NaN, aby nebyly započítány do průměru
    vol_arr[vol_arr < threshold_value] = np.nan

    row = {'Filename': fname}
    for label_name, safe_name, rh_label, lh_label in regions:
        seg_rh = arrayFromSegmentBinaryLabelmap(segNode, f"Segment_{rh_label}", vol_node)
        seg_lh = arrayFromSegmentBinaryLabelmap(segNode, f"Segment_{lh_label}", vol_node)
        mask_r = seg_rh > 0
        mask_l = seg_lh > 0
        row[f"Median_{safe_name}_Right"] = float(np.nanmedian(vol_arr[mask_r])) if mask_r.any() else np.nan
        row[f"Median_{safe_name}_Left"]  = float(np.nanmedian(vol_arr[mask_l])) if mask_l.any() else np.nan
    results.append(row)

# sběr výsledků do DataFrame
df_means = pd.DataFrame(results)
# Uložit mezivýsledek
pd.DataFrame(df_means).to_csv(output_csv, index=False)
print(f"Mean values saved to {output_csv}")

# 4) Vytvořit Prefix pro spojení: UNIS_YYYYMMDD
split_cols = df_means['Filename'].str.split('_', expand=True)
df_means['Prefix'] = split_cols[0].astype(str).str.strip() + '_' + split_cols[1]

# 5) Načíst a připravit klinická data
# CSV má 'UNIS' a 'Study Date' (Excel serial date)
df_pet = pd.read_csv(pet_normy_csv)
if 'UNIS' not in df_pet.columns or 'Study Date' not in df_pet.columns:
    raise KeyError('PET Normy CSV musí obsahovat sloupce UNIS a Study Date.')
df_pet['UNIS_str'] = df_pet['UNIS'].astype(str).str.strip()
df_pet['DateStr'] = pd.to_datetime(
    df_pet['Study Date'], unit='D', origin='1899-12-30'
).dt.strftime('%Y%m%d')
df_pet['Prefix'] = df_pet['UNIS_str'] + '_' + df_pet['DateStr']

# 6) Merge podle Prefix
df_merged = pd.merge(
    df_means,
    df_pet,
    on='Prefix',
    how='left',
    suffixes=('', '_orig'),
    indicator=True
)
unmatched = df_merged[df_merged['_merge'] != 'both']
if not unmatched.empty:
    print(f"Warning: {len(unmatched)} records did not match clinical data by Prefix.")
# Odstranit pomocné sloupce
df_merged.drop(columns=['_merge','UNIS_str','DateStr'], inplace=True, errors='ignore')

# 7) Uložit spojenou tabulku
merged_csv = output_csv.replace('.csv','_merged.csv')
df_merged.to_csv(merged_csv, index=False)
print(f"Merged table saved to {merged_csv}")
