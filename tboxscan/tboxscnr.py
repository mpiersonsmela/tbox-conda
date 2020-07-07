#! /usr/bin/env python
#Tbox-scan
#Version 1.0.0.
#Updated June 24, 2020
#By Merrick Pierson Smela and Jorge A. Marchand
#Harvard Medical School - Department of Genetics
#George M. Church Lab

import os
import pandas as pd
from os import path
import sys
import pkg_resources


#Input parameters
infile = sys.argv[1]
outfile = sys.argv[2]
infernal = sys.argv[4]
logfile = sys.argv[5]
verbose = sys.argv[6]
silence = sys.argv[7]
cutoff = sys.argv[8]

#Covariance model choice
if sys.argv[3]=='1':
    cm = pkg_resources.resource_filename('tboxscan', 'data/RF00230.cm')
    mode='classI'
elif sys.argv[3]=='2':
    cm = pkg_resources.resource_filename('tboxscan', 'data/TBDB001.cm')
    mode='classII'
else:
    print('Invalid covariance model supplied. Try using -m 1 or -m 2 as flags')
    exit()
    
#Initialize
print('\nRunning T-box scanner on '+infile+' using covariance model for '+mode+' t-boxes'+'\n')

#Check files exist
if path.exists(infile)==False:
    print('Error: Input file '+infile+' does not exist.')

if path.exists(cm)==False:
    print('Error: Covariance model '+cm+' does not exist.')

if path.exists(pkg_resources.resource_filename('tboxscan', 'data/rccodonLUT.csv'))==False:
    print('Error: Missing amino acid LUT file '+pkg_resources.resource_filename('tboxscan', 'data/rccodonLUT.csv')+'.')


#Look up table for amino acid family predictions
aalut=pd.read_csv(pkg_resources.resource_filename('tboxscan', 'data/rccodonLUT.csv'))
 
#Remove previous out file
os.system('rm '+outfile+' >/dev/null 2> /dev/null')


#Run CMsearch and feature identification depending on chosen model
if mode == "classI":
    os.system('cmsearch --notrunc --notextw '+cm+' '+infile+' > '+infernal)
    print("INFERNAL output saved to " +infernal+ " Extracting features . . .")
    os.system('python -m tboxscan.pipeline_master '+infernal+' '+outfile+' '+infile+' $3 > '+logfile)
elif mode == "classII":
    os.system('cmsearch --notrunc --notextw '+cm+' '+infile+' > '+infernal)
    print("INFERNAL output saved to " +infernal+ " Extracting features . . .")
    os.system('python -m tboxscan.pipeline_translational '+infernal+' '+outfile+' '+infile+' $3 > '+logfile)


#Read output file
try:
    out = pd.read_csv(outfile)
    #Perform aa family lookup using LUT
    aalist=[None]*len(out)
    for i in range (0,len(out)):
        try:
            aalist[i]='T-box '+aalut.loc[aalut['AC']==out['codon'].iloc[i],['AA']].values[0][0]
        except IndexError: #codon not in LUT
            aalist[i]="Unknown"
    out['AA']=aalist

    #Parse for locus from output file, generate print out response
    out['Locus'] = out['Name'].str.split(':').str[-1]
    outmess=out[['Locus','Score','AA', 'codon_region', 'codon', 'discriminator']]
    outmess.columns=['Locus', 'Score', 'AA Family', 'Spec_Region', 'Specifier', 'T-box Seq']


    #Apply cutoff
    out=out.loc[out['Score'] >= float(cutoff)]
    outmess=outmess.loc[outmess['Score'] >= float(cutoff)]

    #Write output
    if verbose == 'True':
        out.to_csv(outfile)
    elif verbose == 'False':
        outmess.to_csv(outfile)

    #Print output
    if silence == 'False':
        print(outmess)
        print('\n\n')
        
        
except:
    print('Error: Failed to detect T-boxes in '+infile+' using '+cm)


