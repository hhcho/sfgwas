# SF-GWAS

Software for secure and federated genome-wide association studies, described in:

**Secure and Federated Genome-Wide Association Studies for Biobank-Scale Datasets**\
Hyunghoon Cho, David Froelicher, Manaswitha Edupalli, Apostolos Pyrgelis, Juan R. Troncoso-Pastoriza, Jean-Pierre Hubaux, Bonnie Berger\
Under review, 2022

## Installation

### Dependencies

- [Go](https://go.dev/doc/install) (>=1.18.3)
- Python (>=3.9.2) with [NumPy](https://numpy.org/install/)
- [PLINK2](https://www.cog-genomics.org/plink/2.0/)

SF-GWAS assumes that `go`, `python3`, and `plink2` are available in the exec path in shell.

### Get required Go libraries

1. To obtain a modified Lattigo library used by SF-GWAS, run:
```
git clone https://github.com/hcholab/lattigo.git
cd lattigo
git checkout lattigo_pca
cd ..
```
2. Next, obtain a library with core routines for secure multiparty computation by running:
```
git clone https://github.com/hhcho/mpc-core
```

### Obtain SF-GWAS
```
git clone https://github.com/hhcho/sfgwas
cd sfgwas
go get github.com/hhcho/sfgwas-private
go build
```

Note: If the `lattigo` and `mpc-core` repos were cloned to a different location,
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
###
###

## Contact

Hoon Cho, hhcho@broadinstitute.org
