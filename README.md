# SF-GWAS

Software for secure and federated genome-wide association studies, described in:

**Secure and Federated Genome-Wide Association Studies for Biobank-Scale Datasets**\
Hyunghoon Cho, David Froelicher, Manaswitha Edupalli, Apostolos Pyrgelis, Juan R. Troncoso-Pastoriza, Jean-Pierre Hubaux, Bonnie Berger\
Under review, 2022

## Installation

### Install dependencies

- [Go](https://go.dev/doc/install) (>=1.18.3)
- Python (>=3.9.2) with [NumPy](https://numpy.org/install/)
- [PLINK2](https://www.cog-genomics.org/plink/2.0/)

SF-GWAS assumes that `go`, `python3`, and `plink2` are available in the exec path in shell.

### Obtain required Go libraries

1. To obtain a forked [Lattigo](https://github.com/tuneinsight/lattigo) library used by SF-GWAS, run:
```
git clone https://github.com/hcholab/lattigo.git
cd lattigo
git checkout lattigo_pca
cd ..
```
Note the switch to the `lattigo_pca` branch.

2. Next, obtain our library with core routines for secure multiparty computation:
```
git clone https://github.com/hhcho/mpc-core
```

### Install SF-GWAS
```
git clone https://github.com/hhcho/sfgwas
cd sfgwas
go get github.com/hhcho/sfgwas-private
go build
```

If the `lattigo` and `mpc-core` repos were cloned to a different location,
update `../lattigo` and `../mpc-core` in the following lines of `sfgwas/go.mod`
to point to the correct folders.

```
replace github.com/ldsec/lattigo/v2 => ../lattigo
replace github.com/hhcho/mpc-core => ../mpc-core
```

If `go build` produces an error, run the suggested commands and try again. If the build
finishes without any output, the package has been successfully set up.

## Usage

### Input data

For reference, we provide an example synthetic dataset in `example_data/`, which was generated using
the [genotype data simulation routine](https://zzz.bwh.harvard.edu/plink/simulate.shtml) in PLINK1.9
and converted to the PLINK2 PGEN format.
The example data is split between two parties. Each party's local data is stored in
`party1` and `party2` directories.

Main input data files include:
- `geno/chr[1-22].[pgen|psam|pvar]`: [PGEN files](https://www.cog-genomics.org/plink/2.0/input#pgen) for each chromosome. 
- `pheno.txt`: each line includes the phenotype under study for each sample in the `.psam` file
- `cov.txt`: each line includes a tab-separated list of covariates for each sample in the `.psam` file
- `sample_keep.txt`: a list of sample IDs from the `.psam` file to include in the analysis; to be used with the `--keep` flag in PLINK2 (see [here](https://www.cog-genomics.org/plink/2.0/filter#sample) for specification)

### Preparing additional input files

We provide two scripts in `scripts/` for producing additional input files needed for the program. These can be run as follows.

1. `python3 createSnpInfoFiles.py PGEN_PREFIX OUTPUT_DIR`

Note that `PGEN_PREFIX` is expected to be a format string including `%d` in place of the chromosome number (e.g., `geno/chr%d` for the example dataset), which the script sequentially replaces with 1-22. 

This command generates the following three files in `OUTPUT_DIR`:
- `chrom_sizes.txt`: the number of SNPs for each chromosome
- `snp_ids.txt`: a list of variant IDs
- `snp_pos.txt`: a tab-separated two-column matrix indicating the genomic position of each variant; the first column includes the chromosome number (1-22), and the second column includes the base position

2. `python3 computeGenoCounts.py PGEN_PREFIX SAMPLE_KEEP OUTPUT_DIR`:

`PGEN_PREFIX` is the same as above. `SAMPLE_KEEP` points to the `sample_keep.txt` described before, including a list of sample IDs included in the analysis. 

This script generates `all.gcount.transpose.bin` in `OUTPUT_DIR`, which needs to be passed into SF-GWAS. It is binary file encoding a 6-by-m matrix of allele, genotype, and missingness counts for all m variants in the dataset. Note that the counting is performed only on the local dataset using PLINK2.

These additional files have already been generated for the example dataset.

### Setting the configuration

Example config files are provided in `config/`. There are both global config parameters shared by all parties and party-specific parameters.

### Running the program

An example script `run_example.sh` shows how the SF-GWAS program is run. The script spawns 3 processes in the same machine---one for each of the two data-contributing parties and the third for the auxiliary party coordinating the computation. In practice, each party would run their process on their own machine. We also provided `stop.sh` for terminating the jobs launched by `run_example.sh`.

## Contact

Hoon Cho, hhcho@broadinstitute.org
