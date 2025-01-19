#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from qmplot import manhattanplot
from scipy.stats import chi2

# Script Parameters (Modify as Necessary)
INPUT_FILE = "../out/party1/assoc.txt"  # Path to the input file containing SF-GWAS results
POS_FILE = "../example_data/party1/snp_pos.txt"  # File with a list of genomic positions for input SNPs
GKEEP_FILE = "../cache/party1/gkeep.txt"  # Binary vector indicating which variants passed QC for inclusion in the final output
NUM_INDS = 2000  # Total number of individuals across all parties; required for conversion of statistics
NUM_COV = 5  # Total number of covariates in the analysis; required for conversion of statistics
OUTPUT_FILE = "sfgwas_results.jpg"  # Path to save the generated plot figure

def postprocess_assoc(
    new_assoc_file: str,
    assoc_file: str,
    pos_file: str,
    gkeep_file: str,
    num_ind_total: int,
    num_cov: int,
) -> None:
    # new_assoc_file: Name of new assoc file (processed)
    # assoc_file: Name of original assoc file
    # pos_file: Path to pos.txt
    # gkeep_file: Path to gkeep.txt
    # num_ind_total: Total number of individuals
    # num_cov: Number of covariates

    # Load SNP filter
    gkeep = np.loadtxt(gkeep_file, dtype=bool)

    # Load and check dimension of output association stats
    assoc = np.loadtxt(assoc_file)
    assert len(assoc) == gkeep.sum()

    # Calculate p-values
    t2 = (assoc**2) * (num_ind_total - num_cov) / (1 - assoc**2 + 1e-10)
    log10p = np.log10(chi2.sf(t2, df=1))

    # Append SNP position information and write to a new file
    lineno = 0
    assoc_idx = 0

    with open(new_assoc_file, "w") as out:
        out.write("\t".join(["#CHROM", "POS", "R", "LOG10P"]) + "\n")

        for line in open(pos_file):
            pos = line.strip().split()

            if gkeep[lineno]:
                out.write(pos[0] + "\t" + pos[1] + "\t" + str(assoc[assoc_idx]) + "\t" + str(log10p[assoc_idx]) + "\n")
                assoc_idx += 1

            lineno += 1


def plot_assoc(plot_file: str, new_assoc_file: str) -> None:
    # Load postprocessed assoc file and convert p-values
    tab = pd.read_table(new_assoc_file)
    tab["P"] = 10 ** tab["LOG10P"]

    # Create a Manhattan plot
    plt.figure()
    manhattanplot(
        data=tab,
        suggestiveline=None,  # type: ignore
        genomewideline=None,  # type: ignore
        marker=".",
        xticklabel_kws={"rotation": "vertical"},  # set vertical or any other degrees as you like.
    )
    plt.tight_layout()
    plt.savefig(plot_file)

def main():
    print("Plotting script called...")
    print(f"SF-GWAS output file: {INPUT_FILE}")
    print(f"SNP position file: {POS_FILE}")
    print(f"QC filter file: {GKEEP_FILE}")
    print(f"Number of individuals: {NUM_INDS}")
    print(f"Number of covariates: {NUM_COV}")
    
    processed_input = INPUT_FILE + ".processed"
    postprocess_assoc(processed_input, INPUT_FILE, POS_FILE, GKEEP_FILE, NUM_INDS, NUM_COV)
    plot_assoc(OUTPUT_FILE, processed_input)

    print(f"Plot saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()