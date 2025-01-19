import pandas as pd
import os,sys

def mergePhenoCovWithPCs(phenoCovFile, pcFile, outFile):
    phenoCov = pd.read_csv(phenoCovFile, delim_whitespace=True)
    print(phenoCov.columns)
    pc = pd.read_csv(pcFile, delim_whitespace=True)
    pc.rename(columns={'#FID': 'FID'}, inplace=True)
    
    phenoCovPCs = pd.merge(phenoCov, pc, on=['IID', 'FID'])
    phenoCovPCs = phenoCovPCs[['FID', 'IID', 'AGE', 'SEX', 'AGE_SQUARED', 'CENTER', 'PC1' ,'PC2' ,'PC3','PC4', 'PC5' ,'PC6' ,'PC7' ,'PC8' ,'PC9', 'PC10', 'BMI']] 
    phenoCovPCs.to_csv(outFile, sep="\t", index=False)
    print("Written", outFile)

def main():
    phenoCovFile = sys.argv[1]
    pcFile = sys.argv[2]
    outFile = sys.argv[3]

    mergePhenoCovWithPCs(phenoCovFile, pcFile, outFile)

main()