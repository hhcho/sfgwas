import sys
import os

pgen_filename_template = sys.argv[1] # e.g., "gwas_data_chr%d_wgs"
                                     # Note the formatting string "%d",
                                     # which will be replaced with "1", "2", ..., "22"
out_dir = sys.argv[2] # Output directory where snp_pos.txt and snp_ids.txt will be saved

pos_fname = os.path.join(out_dir, "snp_pos.txt")
ids_fname = os.path.join(out_dir, "snp_ids.txt")
csize_fname = os.path.join(out_dir, "chrom_sizes.txt")

pos_file = open(pos_fname, 'w')
ids_file = open(ids_fname, 'w')
csize_file = open(csize_fname, 'w')

for chr in range(1,23):
    pvar_fname = (pgen_filename_template % chr) + ".pvar"

    with open(pvar_fname, "r") as fp:
        first_line = True
        count = 0
        for line in fp:
            if first_line:
                first_line = False
                continue
            tok = line.split()[:3]
            pos_file.write(f"{tok[0]}\t{tok[1]}\n")
            ids_file.write(f"{tok[2]}\n")
            count += 1
        csize_file.write(f"{count}\n")
 
    print(f"Processed chromosome {chr}")

pos_file.close()
ids_file.close()
csize_file.close()

print("Saved output files:")
print(pos_fname)
print(ids_fname)
print(csize_fname)