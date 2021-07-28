#!/usr/bin/env python
# filter_gatk_vcf.py
# manually filter vcf file (remove variants of low VAF)
# 

import pandas as pd
import os.path
import sys
import gzip
from re import search

sid = snakemake.wildcards
varfile = snakemake.input[0]
outfile = snakemake.output[0]
vcut = float(snakemake.config["vcut"])

cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','VAR']

def my_filt(tx):
    # filter low VAF variants
    # 0/1:1598,1511:3109:99:35467,0,36781
    # 1/2:47,194,114:355:99:7061,1316,4103,3719,0,5394
    #print(tx['#CHROM'],tx['POS'],tx['REF'],tx['ALT'])
    gt = [int(x) for x in tx['VAR'].split(':')[0].split('/')]
    rd = [int(x) for x in tx['VAR'].split(':')[1].split(',')]
    # if no reads detected
    if sum(rd) == 0:
        return 'NoCov'
    vaf = sum([rd[x] for x in gt if x != 0])/sum(rd)
    filt = tx['FILTER']
    if filt == '.':
        if vaf > vcut:
            return 'PASS'
        else:
            return 'LowVAF'
    else:
        if vaf > vcut:
            return filt
        else:
            return '{};LowVAF'.format(filt)

with gzip.open(varfile,'r') as fvar, open(outfile,'w') as fout:
    # save headers
    headers = [x.decode("utf-8") for x in fvar.readlines() if search('^#',x.decode("utf-8"))]
    #fout.write(''.join(headers).encode())
    fout.write(headers[0])
    fout.write('##FILTER=<ID=LowVAF,Description="VAF<35%">\n')
    fout.write('##FILTER=<ID=NoCov,Description="No reads covered">\n')
    fout.write(''.join(headers[1:]))
    # load variants
    variants = pd.read_table(varfile, header=None, sep='\t', names=cols, low_memory=False, comment='#')
    print('raw variants: {}'.format(variants.shape))
    # filter based on VAF
    ##variants.drop('FILTER', axis=1, inplace=True)
    variants['FILTER'] = variants.apply(my_filt, axis=1)
    print('filtered variants: {}'.format(variants['FILTER'].str.contains('PASS').sum()))
    # write to file
    variants.to_csv(fout, sep='\t', columns=cols, header=False, index=False)

print('Complete!')

