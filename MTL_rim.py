import nibabel as nib
import numpy as np
import os
from pathlib import Path

# -----------------------
# Paths
# -----------------------
mtl_base_path = "/Volumes/SFB/B04/MTLCIRCUITS/Layers"
ashs_base_path = "/Volumes/SFB/B04/MTLCIRCUITS/Layers"
output_base_path = "/Volumes/SFB/B04/MTLCIRCUITS/Layers"

# Get all subject folders
subject_folders = [d for d in os.listdir(mtl_base_path) if d.startswith("sub-")]


# -----------------------
# Function to label rim passes
# -----------------------
def rim_pass(volume, roi_mask, pass_type="inf_sup", label=1, is_right=False):
    """
    volume: the full rim+body array to label
    roi_mask: boolean mask of rim voxels for this ROI
    pass_type: "inf_sup", "sup_inf", "med_lat", "lat_med"
    label: label to assign to the voxel
    is_right: if True, swap lat_med and med_lat (mirror for right hemisphere)
    """
    S = volume.shape[2]

    # Swap lat_med and med_lat for right hemisphere
    if is_right:
        if pass_type == "lat_med":
            pass_type = "med_lat"
        elif pass_type == "med_lat":
            pass_type = "lat_med"

    for z in range(S):
        slice_2d = roi_mask[:, :, z]
        if pass_type == "inf_sup":
            rows = range(slice_2d.shape[0])
        elif pass_type == "sup_inf":
            rows = reversed(range(slice_2d.shape[0]))
        elif pass_type == "med_lat":
            for col in range(slice_2d.shape[1]):
                rows_in_col = np.where(slice_2d[:, col])[0]
                if rows_in_col.size == 0:
                    continue
                first_voxel = rows_in_col[0]
                if volume[first_voxel, col, z] == 4:
                    volume[first_voxel, col, z] = label
            continue
        elif pass_type == "lat_med":
            for col in reversed(range(slice_2d.shape[1])):
                rows_in_col = np.where(slice_2d[:, col])[0]
                if rows_in_col.size == 0:
                    continue
                first_voxel = rows_in_col[-1]
                if volume[first_voxel, col, z] == 4:
                    volume[first_voxel, col, z] = label
            continue
        else:
            raise ValueError("Unknown pass_type")

        for row in rows:
            cols = np.where(slice_2d[row, :] != 0)[0]
            if cols.size == 0:
                continue
            if pass_type in ["inf_sup", "med_lat"]:
                first_voxel = cols[0]
            else:  # sup_inf
                first_voxel = cols[-1]  # last voxel along row
            if volume[row, first_voxel, z] == 4:
                volume[row, first_voxel, z] = label
    return volume


def correct_ca1_ca2_sub_pattern(volume, roi_mask):
    """
    For CA1+CA2+SUB, ensure rows with 2 runs follow pattern: 2,1,1,2 (inf-sup pass)
    Only relabels the boundary voxels (first and last of each run) - voxels bordering zeros.
    """
    S = volume.shape[2]

    for z in range(S):
        slice_2d = volume[:, :, z].copy()
        roi_slice = roi_mask[:, :, z]

        # Iterate through rows (inferior to superior, scanning left to right)
        for row in range(slice_2d.shape[0]):
            # Find all runs in this row (any non-zero voxel in roi_slice)
            runs = []
            col = 0
            while col < slice_2d.shape[1]:
                if roi_slice[row, col]:
                    run_left = col
                    while col < slice_2d.shape[1] and roi_slice[row, col]:
                        col += 1
                    run_right = col - 1
                    runs.append((run_left, run_right))
                else:
                    col += 1

            # If exactly 2 runs, apply alternating pattern to boundary voxels only
            if len(runs) == 2:
                run1_left, run1_right = runs[0]
                run2_left, run2_right = runs[1]

                # Pattern: 2,1,1,2 (only for the boundary voxels)
                slice_2d[row, run1_left] = 2  # First voxel of run 1
                slice_2d[row, run1_right] = 1  # Last voxel of run 1
                slice_2d[row, run2_left] = 1  # First voxel of run 2
                slice_2d[row, run2_right] = 2  # Last voxel of run 2
        volume[:, :, z] = slice_2d

    return volume



# -----------------------
# Process each subject and side
# -----------------------
for subject in subject_folders:
    for side in ["left", "right"]:
        is_right = (side == "right")
        #side_label = "left" if side == "LH" else "right"
        #output_side = "L" if side == "LH" else "R"

        # Build file paths
        mtl_path = os.path.join(mtl_base_path, subject, f"{subject}_{side}_MTL_binary_mask_no_DG.nii.gz")
        
        
        ashs_path = os.path.join(ashs_base_path, subject, f"{subject}_ASHS_{side}_upsampled.nii.gz")
        output_path = os.path.join(output_base_path, subject, f"{subject}_ASHS_{side}_upsampled_rim.nii")

        # Check if files exist
        #if not os.path.exists(mtl_path):
            #print(f"Warning: MTL file not found: {mtl_path}")
            #continue
        if not os.path.exists(ashs_path):
            print(f"Warning: ASHS file not found: {ashs_path}")
            continue

        print(f"Processing {subject} {side}...")

        # -----------------------
        # Load data
        # -----------------------
        nii_ashs = nib.load(ashs_path)
        ashs = nii_ashs.get_fdata().astype(np.uint8)
        
        #Quickly remove DG
        binary_mask = (ashs != 0) & (ashs != 3)
        mtl = binary_mask.astype(np.uint8)
        
        
        S = mtl.shape[2]

        # -----------------------
        # Generate 1-voxel-thick rim
        # -----------------------
        rim = np.zeros_like(mtl, dtype=np.uint8)
        for z in range(S):
            s l = mtl[:, :, z]
            up = np.zeros_like(sl);
            up[1:] = sl[:-1]
            down = np.zeros_like(sl);
            down[:-1] = sl[1:]
            left = np.zeros_like(sl);
            left[:, 1:] = sl[:, :-1]
            right = np.zeros_like(sl);
            right[:, :-1] = sl[:, 1:]
            neighbor_sum = up + down + left + right
            border = (sl == 1) & (neighbor_sum < 4)
            rim[:, :, z][sl == 1] = 3  # interior
            rim[:, :, z][border] = 4  # rim

        # -----------------------
        # Initialize final volume
        # -----------------------
        final_rim = rim.copy()

        # -----------------------
        # BA35 + BA36 + ERC
        # -----------------------
        roi_labels = [10, 11, 9]  # BA35, BA36, ERC
        roi_mask = np.isin(ashs, roi_labels) & (rim == 4)

        # Inferior→Superior: label 1
        final_rim = rim_pass(final_rim, roi_mask, "inf_sup", 1, is_right)
        # Superior→Inferior: label 2
        final_rim = rim_pass(final_rim, roi_mask, "sup_inf", 2, is_right)
        # Medial→Lateral: label 1
        final_rim = rim_pass(final_rim, roi_mask, "med_lat", 1, is_right)

        # -----------------------
        # SUB +  ERC (flip labels)
        # -----------------------
        roi_labels = [8, 9]  # SUB
        roi_mask = np.isin(ashs, roi_labels) & (rim == 4)
        # Medial→Lateral: label 1
        final_rim = rim_pass(final_rim, roi_mask, "med_lat", 1, is_right)

        # -----------------------
        # CA1 + SUB
        # -----------------------
        roi_labels = [1, 8]  # CA1, SUB
        roi_mask = np.isin(ashs, roi_labels) & (rim == 4)

        final_rim = rim_pass(final_rim, roi_mask, "inf_sup", 2, is_right)

        # -----------------------
        # BA36
        # -----------------------
        roi_labels = [11]  # BA36
        roi_mask = np.isin(ashs, roi_labels) & (rim == 4)
        # Lateral→Medial: label 1
        final_rim = rim_pass(final_rim, roi_mask, "lat_med", 1, is_right)

        # -----------------------
        # BA35
        # -----------------------
        roi_labels = [10]  # BA35
        roi_mask = np.isin(ashs, roi_labels) & (rim == 4)
        # Lateral→Medial: label 1
        final_rim = rim_pass(final_rim, roi_mask, "inf_sup", 1, is_right)

        # -----------------------
        # PHC
        # -----------------------
        roi_labels = [12]  # PHC
        roi_mask = np.isin(ashs, roi_labels) & (rim == 4)

        final_rim = rim_pass(final_rim, roi_mask, "inf_sup", 1, is_right)
        final_rim = rim_pass(final_rim, roi_mask, "sup_inf", 2, is_right)

        # -----------------------
        # Tail
        # -----------------------
        roi_labels = [5]  # PHC
        roi_mask = np.isin(ashs, roi_labels) & (rim == 4)

        final_rim = rim_pass(final_rim, roi_mask, "inf_sup", 2, is_right)
        final_rim = rim_pass(final_rim, roi_mask, "sup_inf", 1, is_right)

        # -----------------------
        # SUB (flip labels)
        # -----------------------
        roi_labels = [8]  # SUB
        roi_mask = np.isin(ashs, roi_labels) & (rim == 4)
        # Superior→Inferior: label 1
        final_rim = rim_pass(final_rim, roi_mask, "sup_inf", 1, is_right)
        # Inferior→Superior: label 2
        final_rim = rim_pass(final_rim, roi_mask, "inf_sup", 2, is_right)

        # -----------------------
        # CA1 (first pass: Lateral→Medial)
        # -----------------------
        roi_labels = [1]  # CA1
        roi_mask = np.isin(ashs, roi_labels) & (rim == 4)

        # Lateral→Medial: label 2
        final_rim = rim_pass(final_rim, roi_mask, "lat_med", 2, is_right)
        # Superior→Inferior: label 2
        roi_labels = [1, 4]  # CA1 + CA3
        roi_mask = np.isin(ashs, roi_labels) & (rim == 4)
        final_rim = rim_pass(final_rim, roi_mask, "sup_inf", 2, is_right)

        # Medial→Lateral: label 1
        roi_labels = [1]
        roi_mask = np.isin(ashs, roi_labels) & (rim == 4)
        final_rim = rim_pass(final_rim, roi_mask, "med_lat", 1, is_right)
        # Inferior→Superior: label 2
        final_rim = rim_pass(final_rim, roi_mask, "inf_sup", 1, is_right)

        # -----------------------
        # CA2  (first pass: Inferior→Superior)
        # -----------------------
        roi_labels = [2]  # CA2
        roi_mask = np.isin(ashs, roi_labels) & (rim == 4)

        # Inferior→Superior: label 2
        final_rim = rim_pass(final_rim, roi_mask, "inf_sup", 1, is_right)
        # Superior→Inferior: label 2
        final_rim = rim_pass(final_rim, roi_mask, "sup_inf", 2, is_right)

        # -----------------------
        # CA3  (first pass: Inferior→Superior)
        # -----------------------
        roi_labels = [4]  # CA3
        roi_mask = np.isin(ashs, roi_labels) & (rim == 4)

        # Inferior→Superior: label 2
        final_rim = rim_pass(final_rim, roi_mask, "inf_sup", 1, is_right)
        # Superior→Inferior: label 2
        final_rim = rim_pass(final_rim, roi_mask, "sup_inf", 2, is_right)

        # -----------------------
        # Post-processing: CA1 + SUB alternating pattern correction
        # -----------------------


        final_rim[final_rim == 4] = 3

        print("Now post-processing CA1 + SUB")
        roi_labels = [1, 2, 5, 8 ]  # CA1, SUB
        roi_mask = np.isin(ashs, roi_labels)
        final_rim = correct_ca1_ca2_sub_pattern(final_rim, roi_mask)



        # -----------------------
        # Save final combined volume
        # -----------------------

        # -----------------------
        # Save final combined volume
        # -----------------------
        nib.save(nib.Nifti1Image(final_rim, nii_ashs.affine, nii_ashs.header), output_path)
        print(f"  Saved: {output_path}")

print("All subjects processed!")
'''
if combine_layers:
    base_dir = "/Volumes/SFB/B04/MTLRSFC/data/nifti/derivatives/LAYERS/Layer_Files"
    for subj in subjects:
        for side in sides:
            file = os.path.join(base_dir, subj, f"{subj}_HC_{side}_layers_layers_equidist.nii")
            img_data = nib.load(file)
            img= img_data.get_fdata()
            affine = img_data.affine
            header = img_data.header
            halfs = np.where((img >= 1) & (img <= 10), 1,
                             np.where((img >= 11) & (img <= 20), 2, 0))

            # Save corrected halfs
            out_file = os.path.join(base_dir, subj, f"{subj}_HC_{side}_layers_N2.nii")
            nib.save(nib.Nifti1Image(halfs, affine, header), out_file)
            print(f"Saved two-half file: {out_file}")
'''





