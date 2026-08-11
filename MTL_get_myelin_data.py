import os
import numpy as np
import pandas as pd
import nibabel as nib

'''
    Extracts subfield-mean, subfield-laminar-mean, and subfield-zone-mean
    qT1 values from cropped/upsampled qT1 maps, using ASHS subfield masks
    and LayNii laminar layer masks. Output: one row per subject, one
    column per (ROI x layer/zone x hemisphere) combination.
'''

base = "/Volumes/SFB/B04/MTLCIRCUITS/Layers"
subs = os.listdir(base)
hemis = ["left", "right"]

rois = {
    "ERC": 9,
    "SUB": 8,
    "CA1": 1,
    "CA2": 2,
    "CA3": 4,
    "DG": 3
}


def assign_zone(layer):
    if 1 <= layer <= 7:
        return "Inner"
    elif 8 <= layer <= 14:
        return "Middle"
    elif 15 <= layer <= 21:
        return "Outer"
    return None


def compute_hemi_means(hemi, roi_label_dict, qt1_array, ashs_array, layers_array):
    """
    Computes subfield mean, subfield-by-layer mean, and subfield-by-zone
    mean qT1 values for one subject, one hemisphere.
    Returns a dict of {column_name: value}.
    """
    row = {}

    layer_vals = np.unique(layers_array[layers_array != 0])

    # Group layers into zones using assign_zone, once per hemisphere
    zone_layers = {"Inner": [], "Middle": [], "Outer": []}
    for layer in layer_vals:
        zone = assign_zone(int(layer))
        if zone is not None:
            zone_layers[zone].append(layer)

    for roi_name, roi_val in roi_label_dict.items():
        roi_mask = ashs_array == roi_val

        # --- Subfield mean (whole ROI, all layers combined) ---
        if roi_mask.sum() == 0:
            row[f"mean_{roi_name}_{hemi}"] = np.nan
        else:
            row[f"mean_{roi_name}_{hemi}"] = np.mean(qt1_array[roi_mask])

        # --- Subfield-by-layer means ---
        for layer in layer_vals:
            layer_mask = roi_mask & (layers_array == layer)
            col = f"mean_{roi_name}_layer{int(layer)}_{hemi}"
            row[col] = np.mean(qt1_array[layer_mask]) if layer_mask.sum() > 0 else np.nan

        # --- Subfield-by-zone means (Inner/Middle/Outer) ---
        for zone_name, zone_layer_list in zone_layers.items():
            zone_mask = roi_mask & np.isin(layers_array, zone_layer_list)
            col = f"mean_{roi_name}_{zone_name}_{hemi}"
            row[col] = np.mean(qt1_array[zone_mask]) if zone_mask.sum() > 0 else np.nan

    return row


rows = []

for sub in subs:
    row = {"ID": sub}

    for hemi in hemis:
        qt1_file = os.path.join(base, sub, f"{sub}_{hemi}_cropped_upsampled.nii.gz")
        ashs_file = os.path.join(base, sub, f"{sub}_ASHS_{hemi}_upsampled.nii.gz")
        layer_file = os.path.join(base, sub, "layer_files", f"sub-{sub}_{hemi}_N21_layers_equivol.nii")

        if not (os.path.exists(qt1_file) and os.path.exists(ashs_file) and os.path.exists(layer_file)):
            print(f"Skipping {sub} ({hemi}): missing one or more input files.")
            continue

        qt1_img = nib.load(qt1_file).get_fdata()
        ashs_img = nib.load(ashs_file).get_fdata()
        layers_img = nib.load(layer_file).get_fdata()

        hemi_row = compute_hemi_means(hemi, rois, qt1_img, ashs_img, layers_img)
        row.update(hemi_row)

    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv("/Volumes/SFB/B04/MTLCIRCUITS/Tables/MTLCIRCUITS_7T.csv", index=False)
print(df.head())
		
        
			
				
