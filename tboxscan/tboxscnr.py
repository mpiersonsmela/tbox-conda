#! /usr/bin/env python
#Tbox-scan
#Version 0.1.3
#Updated July 8, 2020
#By Merrick Pierson Smela and Jorge A. Marchand
#Harvard Medical School - Department of Genetics
#George M. Church Lab

import os
import pandas as pd
import sys
#import pkg_resources


#Input parameters
parent_path = os.path.dirname(sys.argv[0])
infile = sys.argv[1]
outfile = sys.argv[2]
infernal = sys.argv[4]
logfile = sys.argv[5]
verbose = sys.argv[6]
silence = sys.argv[7]
cutoff = sys.argv[8]

#Covariance model choice
if sys.argv[3]=='1':
    cm = os.path.join(parent_path,'data','RF00230.cm')
    mode='classI'
elif sys.argv[3]=='2':
    cm = os.path.join(parent_path,'data','TBDB001.cm')
    mode='classII'
else:
    print('Invalid covariance model supplied. Try using -m 1 or -m 2 as flags')
    sys.exit(1)
    
#Initialize
print('\nRunning T-box scanner on '+infile+' using covariance model '+cm+' for '+mode+' t-boxes'+'\n')

#Check files exist
if os.path.exists(infile)==False:
    print('Error: Input file '+infile+' does not exist.')
    sys.exit(1)

if os.path.exists(cm)==False:
    print('Error: Covariance model '+cm+' does not exist.')
    sys.exit(1)

aa_lut_path = os.path.join(parent_path,'data','rccodonLUT.csv')

if not os.path.exists(aa_lut_path):
    print('Error: Missing amino acid LUT file: '+aa_lut_path)
    sys.exit(1)

#Look up table for amino acid family predictions
aalut=pd.read_csv(aa_lut_path)
 
#Remove previous out file
os.system('rm '+outfile+' >/dev/null 2> /dev/null')


#Run CMsearch and feature identification depending on chosen model
if mode == "classI":
    os.system('cmsearch --notrunc --notextw '+cm+' '+infile+' > '+infernal)
    print("INFERNAL output saved to " +infernal+ " Extracting features . . .")
    os.system('python '+parent_path+'/pipeline_master.py '+infernal+' '+outfile+' '+infile+' $3 > '+logfile)
elif mode == "classII":
    os.system('cmsearch --notrunc --notextw '+cm+' '+infile+' > '+infernal)
    print("INFERNAL output saved to " +infernal+ " Extracting features . . .")
    os.system('python '+parent_path+'/pipeline_translational.py '+infernal+' '+outfile+' '+infile+' $3 > '+logfile)


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
    sys.exit(1)
