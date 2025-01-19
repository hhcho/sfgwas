import numpy as np
import sys
import os

pgen_filename_template = sys.argv[1] # e.g., "gwas_data_chr%d_wgs"
                                     # Note the formatting string "%d",
                                     # which will be replaced with "1", "2", ..., "22"
sample_keep_file = sys.argv[2] # Samples to keep, in a format expected by PLINK2 with the "--keep" flag
out_dir = sys.argv[3] # Output directory, will be created if it does not exist

os.system(f"mkdir -p {out_dir}")

all_fname = os.path.join(out_dir, "all.gcount")
all_file = open(all_fname, "w")

for chr in range(1,23):
    pgen_prefix = pgen_filename_template % chr
    out_prefix = os.path.join(out_dir, f"chr{chr}")
    
    os.system(f"plink2 --threads 1 --pfile {pgen_prefix} --keep {sample_keep_file} --geno-counts --out {out_prefix}")

    print(f"Geno counts computed for chromosome {chr}")

    out_file = f"{out_prefix}.gcount"

    with open(out_file, "r") as fp:
        first_line = True
        for line in fp:
            if first_line:
                first_line = False
                continue
            tok = line.split("\t")
            all_file.write("\t".join(tok[4:10]) + "\n")

    print(f"Geno counts for chromosome {chr} added to: {all_fname}")    

all_file.close()

print(f"Processing all geno counts: {all_fname}")

x = np.loadtxt(all_fname, dtype=np.uint32)

print("Dimensions:", x.shape, x.dtype)

x=x.transpose()

print("After transpose:", x.shape, x.dtype)

all_bin_fname = all_fname + ".transpose.bin"
x.tofile(all_bin_fname)

print("Saved output binary file to:", all_bin_fname)