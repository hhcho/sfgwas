import pandas as pd
import sys

def getBiAllelicVariants(snpInFile, snpOutFile):
    locusInfo = pd.read_csv(snpInFile, index_col=False, delim_whitespace=True,
                            usecols=['#CHROM', 'POS'])
    
    # Remove any mult-allelic SNPs (duplicate position values)
    locusInfo.drop_duplicates(subset=['POS'], keep=False, inplace=True, ignore_index=True)

    locusInfo.rename(columns={"#CHROM": "CHR", "POS": "START"}, inplace=True)
    locusInfo['END'] = locusInfo.loc[:, 'START'].copy()
    locusInfo['SET_ID'] = locusInfo.index + 1
    locusInfo.to_csv(snpOutFile, sep="\t", index=False)


def main():
    snpInFile = sys.argv[1]
    snpOutFile = sys.argv[2]

    getBiAllelicVariants(snpInFile, snpOutFile)

main()