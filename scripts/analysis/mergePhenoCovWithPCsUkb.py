import pandas as pd
import os,sys

def mergePhenoCovWithPCs(phenoCovFile, pcFile, outFile):
    phenoCov = pd.read_csv(phenoCovFile, delim_whitespace=True)
    print(phenoCov.columns)
    pc = pd.read_csv(pcFile, delim_whitespace=True)
    pc.rename(columns={'#FID': 'FID'}, inplace=True)
    
    phenoCov = phenoCov.iloc[:, 1:]
    pc = pc.iloc[:, 1:]

    phenoCovPCs = pd.merge(phenoCov, pc, on=['IID'])
    phenoCovPCs['FID'] = phenoCovPCs.IID.copy()

    phenoCovPCs = phenoCovPCs[['FID', 'IID', 'Center', 'Age', 'Sex', 'Age_Sex' ,'Age_Squared', 'Age_Squared_Sex' ,'PC1' ,'PC2' ,'PC3','PC4', 'PC5' ,'PC6' ,'PC7' ,'PC8' ,'PC9', 'PC10', 'BMI']] 
    phenoCovPCs.to_csv(outFile, sep="\t", index=False)

def main():
    phenoCovFile = sys.argv[1]
    pcFile = sys.argv[2]
    outFile = sys.argv[3]

    mergePhenoCovWithPCs(phenoCovFile, pcFile, outFile)

main()
