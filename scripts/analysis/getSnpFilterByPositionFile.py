import pandas as pd
import os, sys

# Creates file for use with plink2 --extract bed1 flag to filter SNPs based on chromosome and position
# plink2 --extract only filters based on rsID
def formatSnpFilterByLocus(inPath, snpInFile, outPath, snpOutFile):
    # Read in .pvar file of SNPs that passed QC
    qcSnps = pd.read_csv(inPath + snpInFile, delim_whitespace=True, usecols=['#CHROM', 'POS'])
    
    qcSnps.rename(columns={'#CHROM': "CHR"}, inplace=True)
    qcSnps['START'] = qcSnps.POS.copy()
    qcSnps['END'] = qcSnps.POS.copy()
    qcSnps['SET_ID'] = range(1, qcSnps.shape[0] + 1)
    qcSnps.drop(columns=['POS'], inplace=True)

    qcSnps = qcSnps[['CHR', 'START', 'END', 'SET_ID']]
    qcSnps.to_csv(outPath + snpOutFile, sep="\t", index=False)

def main():
    inPath=sys.argv[1]
    snpInFile=sys.argv[2]

    outPath=sys.argv[3]
    snpOutFile=sys.argv[4]

    formatSnpFilterByLocus(inPath, snpInFile, outPath, snpOutFile)

main()