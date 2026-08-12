import sys
import numpy as np
import nibabel as nib
from scipy.stats import pearsonr, spearmanr

def load_data(path):
    return nib.load(path).get_fdata()

def dice_coefficient(a, b, threshold=0):
    a_bin = a > threshold
    b_bin = b > threshold
    intersection = np.logical_and(a_bin, b_bin).sum()
    denom = a_bin.sum() + b_bin.sum()
    if denom == 0:
        return np.nan
    return 2.0 * intersection / denom

def main():
    tdi_a_path, tdi_b_path, label_a, label_b, out_csv = sys.argv[1:6]

    tdi_a = load_data(tdi_a_path)
    tdi_b = load_data(tdi_b_path)

    if tdi_a.shape != tdi_b.shape:
        raise ValueError(f"Shape mismatch: {tdi_a.shape} vs {tdi_b.shape}. "
                          f"Both TDIs must be generated with -template on the same reference image.")

    dice = dice_coefficient(tdi_a, tdi_b, threshold=0)

    # Voxel-wise correlation, restricted to voxels where at least one map is nonzero
    mask = (tdi_a > 0) | (tdi_b > 0)
    a_vals = tdi_a[mask]
    b_vals = tdi_b[mask]

    pearson_r, pearson_p = pearsonr(a_vals, b_vals)
    spearman_r, spearman_p = spearmanr(a_vals, b_vals)

    n_voxels_a = int((tdi_a > 0).sum())
    n_voxels_b = int((tdi_b > 0).sum())
    n_voxels_union = int(mask.sum())
    n_voxels_intersection = int(((tdi_a > 0) & (tdi_b > 0)).sum())

    print(f"\n--- Comparison: {label_a} vs {label_b} ---")
    print(f"Nonzero voxels ({label_a}): {n_voxels_a}")
    print(f"Nonzero voxels ({label_b}): {n_voxels_b}")
    print(f"Union: {n_voxels_union}, Intersection: {n_voxels_intersection}")
    print(f"Dice coefficient: {dice:.4f}")
    print(f"Pearson r: {pearson_r:.4f} (p={pearson_p:.3g})")
    print(f"Spearman rho: {spearman_r:.4f} (p={spearman_p:.3g})")

    with open(out_csv, "w") as f:
        f.write("comparison,n_voxels_a,n_voxels_b,n_voxels_union,n_voxels_intersection,dice,pearson_r,pearson_p,spearman_r,spearman_p\n")
        f.write(f"{label_a}_vs_{label_b},{n_voxels_a},{n_voxels_b},{n_voxels_union},{n_voxels_intersection},"
                f"{dice:.6f},{pearson_r:.6f},{pearson_p:.6g},{spearman_r:.6f},{spearman_p:.6g}\n")

if __name__ == "__main__":
    main()
