import sys
import os
import numpy as np

'''
Reads a single tck2connectome connectome.csv (5x5 matrix) for one
timepoint/version/hemisphere, and computes MSP and TSP streamline counts
and proportions.

Called once per (timepoint, version, hemisphere) combination; appends one
row to a running summary CSV.

Node order after remapping (node-creation step): 1=CA1, 2=CA3, 3=DG, 4=SUB, 5=ERC
-> zero-indexed matrix positions: CA1=0, CA3=1, DG=2, SUB=3, ERC=4

Usage:
    python3 quantify_pathways.py <connectome_csv> <timepoint> <version> <hemisphere> <summary_csv>
'''

IDX = {"CA1": 0, "CA3": 1, "DG": 2, "SUB": 3, "ERC": 4}

# version "1" = native-space-warped-to-T1 pipeline, "2" = T1-space/ACT pipeline
VERSION_LABELS = {"1": "native", "2": "t1space_ACT"}


def load_connectome(path):
    return np.loadtxt(path, delimiter=',')


def main():
    conn_csv, timepoint, version, hemi, summary_csv = sys.argv[1:6]

    if not os.path.exists(conn_csv):
        print(f"Connectome not found: {conn_csv} -- skipping quantification.")
        return

    m = load_connectome(conn_csv)

    msp_count = m[IDX["ERC"], IDX["CA1"]]
    tsp_dg_ca3 = m[IDX["DG"], IDX["CA3"]]
    tsp_ca3_ca1 = m[IDX["CA3"], IDX["CA1"]]
    tsp_total = tsp_dg_ca3 + tsp_ca3_ca1

    # Matrix is symmetric (from -symmetric); sum/2 avoids double-counting each edge
    total_mtl = m.sum() / 2.0

    msp_prop = msp_count / total_mtl if total_mtl > 0 else np.nan
    tsp_prop = tsp_total / total_mtl if total_mtl > 0 else np.nan

    space_label = VERSION_LABELS.get(version, version)

    write_header = not os.path.exists(summary_csv)
    with open(summary_csv, "a") as out_f:
        if write_header:
            out_f.write(
                "timepoint,version,space,hemisphere,msp_count,tsp_dg_ca3_count,tsp_ca3_ca1_count,"
                "tsp_total_count,total_mtl_streamlines,msp_proportion,tsp_proportion\n"
            )
        out_f.write(
            f"{timepoint},{version},{space_label},{hemi},{msp_count:.1f},{tsp_dg_ca3:.1f},{tsp_ca3_ca1:.1f},"
            f"{tsp_total:.1f},{total_mtl:.1f},{msp_prop:.6f},{tsp_prop:.6f}\n"
        )

    print(f"[{timepoint}_{version} {hemi}] MSP={msp_count:.0f}, TSP={tsp_total:.0f} "
          f"(DG-CA3={tsp_dg_ca3:.0f}, CA3-CA1={tsp_ca3_ca1:.0f}), "
          f"total MTL streamlines={total_mtl:.0f}, "
          f"MSP prop={msp_prop:.4f}, TSP prop={tsp_prop:.4f}")


if __name__ == "__main__":
    main()
