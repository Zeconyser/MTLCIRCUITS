"""
crop_upsample_mtl.py

For each subject/hemisphere:
  1. Load ASHS segmentation and qT1 map.
  2. Clean ASHS labels to keep only the desired MTL subfields.
  3. Crop qT1 (and ASHS) tightly to the bounding box of the cleaned labels.
  4. Upsample both the cropped qT1 and cropped ASHS to 0.2mm isotropic,
     updating the affine consistently so both stay aligned.
  5. Save outputs.
"""

import os
import numpy as np
import nibabel as nib
from nibabel.affines import voxel_sizes
from scipy.ndimage import zoom

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

BASE = "/Volumes/SFB/B04/MTLCIRCUITS"
ASHS = os.path.join(BASE, "Masks")
ANAT = os.path.join(BASE, "Anat")
LAYERS = os.path.join(BASE, "Layers") 
hemis = ["left", "right"]

crop_prefix = "rqT1_mtl_crop_"
TARGET_RES = 0.2  # mm, isotropic
PAD = 0            # voxels of padding around bounding box, in original resolution

labels_left = {
    "CA1": 1,
    "CA2": 2,
    "DG": 3,
    "CA3": 4,
    "SUB": 8,
    "ERC": 9,
}
labels_right = {key: int(val) + 100 for key, val in labels_left.items()}
labels = [labels_left, labels_right]


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def crop_to_bbox(img, affine, mask, pad=0):
    """
    Tightly crop img to the bounding box of mask>0, and return the
    affine updated to reflect the new (smaller) origin.

    Parameters
    ----------
    img : np.ndarray
        3D array to crop.
    affine : np.ndarray
        4x4 affine associated with img.
    mask : np.ndarray (bool)
        Same shape as img; True where labels are present.
    pad : int
        Voxels of padding to add around the bounding box (clipped to bounds).

    Returns
    -------
    cropped : np.ndarray
        Tightly cropped array (smaller than img, not zero-padded).
    new_affine : np.ndarray
        4x4 affine updated so the cropped array's origin matches the
        correct physical location.
    bbox : tuple of slices
        The bounding box used, in case you need to apply it to another
        array of the same original shape (e.g. the segmentation).
    """
    if not mask.any():
        raise ValueError("No nonzero labels found in segmentation")

    coords = np.array(np.nonzero(mask))  # shape (3, N)
    mins = coords.min(axis=1)
    maxs = coords.max(axis=1)

    mins = np.maximum(mins - pad, 0)
    maxs = np.minimum(maxs + pad, np.array(img.shape) - 1)

    bbox = tuple(slice(mn, mx + 1) for mn, mx in zip(mins, maxs))
    cropped = img[bbox]

    # Shift the affine's origin to the first voxel of the crop
    new_affine = affine.copy()
    offset_vox = np.array([mins[0], mins[1], mins[2], 1])
    new_affine[:3, 3] = (affine @ offset_vox)[:3]    
    return cropped, new_affine, bbox


def upsample_to_iso(img, affine, target_res=0.2, order=3):
    """
    Resample img to isotropic target_res (mm) and return the correctly
    rescaled affine.

    Parameters
    ----------
    img : np.ndarray
        3D array to resample.
    affine : np.ndarray
        4x4 affine associated with img.
    target_res : float
        Desired isotropic voxel size in mm.
    order : int
        Interpolation order for scipy.ndimage.zoom.
        Use order=3 (cubic) for continuous data (e.g. qT1).
        Use order=0 (nearest neighbor) for label/segmentation data.

    Returns
    -------
    resampled : np.ndarray
    new_affine : np.ndarray
        4x4 affine reflecting the new voxel size.
    """
    current_res = voxel_sizes(affine)  # (vx, vy, vz) in mm
    zoom_factors = current_res / target_res  # >1 means upsampling

    resampled = zoom(
        img,
        zoom_factors,
        order=order,
        mode="nearest" if order == 0 else "mirror",
    )

    # Recompute exact per-axis zoom actually applied (guards against
    # scipy's internal rounding of output shape) so the affine stays
    # perfectly self-consistent with the resampled array's shape.
    exact_zoom = np.array(resampled.shape) / np.array(img.shape)

    new_affine = affine.copy()
    new_affine[:3, :3] = affine[:3, :3] / exact_zoom

    return resampled, new_affine


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

def main():
    subjects = os.listdir(ANAT)

    for sub in subjects:
        print(f"Now Processing Subject: {sub}...")

        qt1_path = os.path.join(ANAT, sub, f"r{sub}_T1map.nii")
        if not os.path.exists(qt1_path):
            print(f"  Skipping {sub}: qT1 not found at {qt1_path}")
            continue
        qt1_nii = nib.load(qt1_path)
        qt1_img = qt1_nii.get_fdata()
        qt1_affine = qt1_nii.affine

        for hemi in range(2):
            print(f"\n  Hemisphere: {hemis[hemi]}...")

            ashs_path = os.path.join(ASHS, sub, f"{sub}_ASHS_{hemis[hemi]}.nii")
            if not os.path.exists(ashs_path):
                print(f"    Skipping: ASHS file not found at {ashs_path}")
                continue
            ashs_nii = nib.load(ashs_path)
            ashs_img = ashs_nii.get_fdata().astype(np.uint8)

            if qt1_img.shape != ashs_img.shape:
                raise ValueError(
                    f"{sub} {hemis[hemi]}: qT1 shape {qt1_img.shape} != "
                    f"ASHS shape {ashs_img.shape}"
                )

            # --- Step 1: clean labels for this hemisphere ---
            print("    Removing unwanted subfields/labels...")
            rois = list(labels_left.values())
            mask = np.isin(ashs_img, rois)
            ashs_clean = np.where(mask, ashs_img, 0)

            # --- Step 2: crop qT1 and ASHS to bounding box of labels ---
            print("    Cropping to label bounding box...")
            qt1_crop, crop_affine, bbox = crop_to_bbox(
                qt1_img, qt1_affine, mask, pad=PAD
            )
            ashs_crop = ashs_clean[bbox]

            # --- Step 3: upsample both to target isotropic resolution ---
            print(f"    Upsampling to {TARGET_RES}mm isotropic...")
            qt1_hi, qt1_hi_affine = upsample_to_iso(
                qt1_crop, crop_affine, target_res=TARGET_RES, order=3
            )
            ashs_hi, ashs_hi_affine = upsample_to_iso(
                ashs_crop, crop_affine, target_res=TARGET_RES, order=0
            )

            # Sanity check: both outputs should share the same affine
            if not np.allclose(qt1_hi_affine, ashs_hi_affine):
                raise RuntimeError(
                    f"{sub} {hemis[hemi]}: affines diverged after upsampling "
                    "-- qT1 and ASHS crop shapes may not have matched."
                )

            sub_out = os.path.join(LAYERS, sub)
            if not os.path.exists(sub_out): 
                os.makedirs(sub_out) 		
            # --- Step 4: save ---
            out_qt1 = os.path.join(
                sub_out, f"{sub}_{hemis[hemi]}_cropped_upsampled.nii.gz"
            )

	  
            out_ashs = os.path.join(
                sub_out, f"{sub}_ASHS_{hemis[hemi]}_upsampled.nii.gz"
            )
		

            print(f"    Saving qT1 crop to {out_qt1}")
            nib.save(nib.Nifti1Image(qt1_hi, qt1_hi_affine), out_qt1)

            print(f"    Saving ASHS crop to {out_ashs}")
            nib.save(
                nib.Nifti1Image(ashs_hi.astype(np.uint8), ashs_hi_affine),
                out_ashs,
            )

        print(f"Done with subject {sub}.\n")


if __name__ == "__main__":
    main()
