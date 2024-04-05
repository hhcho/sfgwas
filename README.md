# SF-GWAS

Software for secure and federated genome-wide association studies, as described in:

**Secure and Federated Genome-Wide Association Studies for Biobank-Scale Datasets**\
Hyunghoon Cho, David Froelicher, Jeffrey Chen, Manaswitha Edupalli, Apostolos Pyrgelis, Juan R. Troncoso-Pastoriza, Jean-Pierre Hubaux, Bonnie Berger\
Under review, 2024

This repository provides the code for PCA-based association analysis. For LMM-based GWAS, see [here](https://github.com/hhcho/sfgwas-lmm). We are working on providing both workflows in one package for convenience.


## Installation

### Dependencies

SF-GWAS requires that `go`, `python3`, and `plink2` are available in the exec path in shell. Here are the links for installation:

- [Go](https://go.dev/doc/install) (>=1.18.3)
- Python (>=3.9.2) with [NumPy](https://numpy.org/install/)
- [PLINK2](https://www.cog-genomics.org/plink/2.0/)

### Go libraries for secure computation

1. SF-GWAS uses the [Lattigo](https://github.com/tuneinsight/lattigo) library for multiparty homomorphic encryption. To install a forked version used by SF-GWAS (branch: `lattigo_pca`), run:
```
git clone https://github.com/hcholab/lattigo.git
cd lattigo
git checkout lattigo_pca
cd ..
```

2. Next, SF-GWAS also uses our own library of secret sharing-based multiparty computation routines. This can be obtained by running:
```
git clone https://github.com/hhcho/mpc-core
```

### Install SF-GWAS

To install SF-GWAS, clone the repository and try building as follows.
```
git clone https://github.com/hhcho/sfgwas
cd sfgwas
go get github.com/hhcho/sfgwas
go build
```

Note that, if `lattigo` and `mpc-core` repos from the previous steps are cloned to a different location,
update `../lattigo` and `../mpc-core` in the following lines of `sfgwas/go.mod`
to point to the correct folders. The paths are relative, starting from the root directory of `sfgwas` repo where the `go.mod` file is located.

```
replace github.com/ldsec/lattigo/v2 => ../lattigo
replace github.com/hhcho/mpc-core => ../mpc-core
```

If `go build` produces an error, run any commands suggested by Go and try again. If the build
finishes without any output, the package has been successfully configured.

## Usage

### Input data

We provide an example synthetic dataset in `example_data/`, generated using the [genotype data simulation routine](https://zzz.bwh.harvard.edu/plink/simulate.shtml) in PLINK1.9
and converted to the PLINK2 PGEN format.
The example data is split between two parties. Each party's local data is stored in
`party1` and `party2` directories. Note that SF-GWAS can be run with more than two parties.

Main input data files include:
- `geno/chr[1-22].[pgen|psam|pvar]`: [PGEN files](https://www.cog-genomics.org/plink/2.0/input#pgen) for each chromosome. 
- `pheno.txt`: each line includes the phenotype under study for each sample in the `.psam` file.
- `cov.txt`: each line includes a tab-separated list of covariates for each sample in the `.psam` file.
- `sample_keep.txt`: a list of sample IDs from the `.psam` file to include in the analysis; to be used with the `--keep` flag in PLINK2 (see [here](https://www.cog-genomics.org/plink/2.0/filter#sample) for file specification).

### Preparing additional input files

We provide two preprocessing scripts in `scripts/` for producing additional input files needed for SF-GWAS. 

1. `createSnpInfoFiles.py` processes the provided `.pvar` files to create a number of files specifying variant information. It can be run as follows:

`python3 createSnpInfoFiles.py PGEN_PREFIX OUTPUT_DIR`

Note that `PGEN_PREFIX` is expected to be a format string including `%d` in place of the chromosome number (e.g., `geno/chr%d` for the example dataset), which the script sequentially replaces with the chromosome numbers 1-22. 

This command generates the following three files in `OUTPUT_DIR`:
- `chrom_sizes.txt`: the number of SNPs for each chromosome
- `snp_ids.txt`: a list of variant IDs
- `snp_pos.txt`: a tab-separated, two-column matrix specifying the genomic position of each variant (chromosome number followed by base position)

2. `computeGenoCounts.py` runs PLINK2's [genotype counting routine](https://www.cog-genomics.org/plink/2.0/basic_stats#geno_counts) (`--geno-counts`) on each chromosome to obtain genotype, allele, and missingness counts for each variant. It is run as follows:

`python3 computeGenoCounts.py PGEN_PREFIX SAMPLE_KEEP OUTPUT_DIR`

`PGEN_PREFIX` is the same as before. `SAMPLE_KEEP` points to the `sample_keep.txt` described above as a main input file, including a list of sample IDs to be included in the analysis in a PLINK2-recognized format.

This script generates `all.gcount.transpose.bin` in `OUTPUT_DIR`, which needs to be provided to SF-GWAS. It is a binary file encoding a 6-by-m matrix of precomputed allele, genotype, and missingness counts for all m variants in the dataset. 

We provide both variant information files and the genotype counts for the example dataset in `example_data/`.

### Setting the configuration

Example config files are provided in `config/`. There are both global config parameters shared by all parties and party-specific parameters.

### Running the program

An example script `run_example.sh` shows how to run SF-GWAS. 

The script spawns 3 processes on the same machine---one for each of the two data-contributing parties (`PID=1` and `PID=2`)and the third for the auxiliary party coordinating the computation (`PID=0`). In practice, each party runs their process on their own machine and provides the IP addresses of other parties in the configuration for network communication. 

We also provide `stop.sh` for terminating the jobs launched by `run_example.sh`.

### Output

Once SF-GWAS finishes, it generates `assoc.txt` in the output directory specified in the configuration. This file includes the Pearson correlation coefficient for each variant passing the quality control filters. Binary vector indicating the inclusion of each variant in the final output is specified in `gkeep.txt` in the cache directory.

## Contact for questions

David Froelicher, dfroelic@mit.edu;
Hoon Cho, hoon.cho@yale.edu
