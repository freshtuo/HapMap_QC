#!/usr/bin/env python
# summarize_gatk_joint_call.py
# summarize benchmark results on GATK joint calls
# 

import sys
import re
import gzip
import pybedtools
import pandas as pd
from contextlib import redirect_stdout

# sample ids
samples = snakemake.config["bams"].keys()
# query bed file
query_bed_file = snakemake.config["bed"]
# reference bed files (dict)
ref_bed_files = snakemake.config["refbeds"]
# original query vcf files
query_vcf_files = snakemake.input[0:len(samples)]
# benchmark summary files
summary_files = snakemake.input[len(samples):2*len(samples)]
# benchmark output vcf files
bench_vcf_files = snakemake.input[2*len(samples):]
# log file
log_file = snakemake.log[0]
# output excel file
out_excel_file = snakemake.output[0]

# functions
def calculate_total_bases(bedfileA, bedfileB):
    """calculate total overlapping bases in two bed files"""
    # load first bed file
    bedA = pybedtools.BedTool(bedfileA)
    # get overlapped regions
    overlapbed = bedA.intersect(bedfileB)
    return overlapbed.sort().total_coverage()

def fetch_confusion_matrix(tx, sid, total):
    """fetch TP/FP/TN/FN; used by apply to a DataFrame"""
    P = tx['TRUTH.TOTAL']
    N = total - P
    PP = tx['QUERY.TOTAL'] - tx['QUERY.UNK']
    TP = tx['TRUTH.TP']
    FN = tx['TRUTH.FN']
    FP = tx['QUERY.FP']
    TN = N - FP
    return pd.Series([sid,total,TP,FN,FP,TN,P,N,PP,TP/P,TN/N,TP/(TP+FP)], index=['Sample','Total','TP','FN','FP','TN','P','N','PP','Sensitivity','Specificity','Precision'])

def process_results(sumfile, sid, total):
    """extract TP/FP/TN/FN, etc from the benchmark summary file"""
    # read in summary file
    bsum = pd.read_csv(sumfile, header=0, low_memory=False)
    # focus on PASS variants only
    bsum = bsum[bsum['Filter'] == 'PASS']
    # collect confusion matrix
    cmtx = bsum.apply(fetch_confusion_matrix, axis=1, args=(sid,total))
    return cmtx

def prepare_overall_table(samples, summary_files, query_bed_file, ref_bed_files):
    """prepare 1st excel sheet: overall stats"""
    print('preparing overall table:')
    cmtx_list = []
    for k,sid in enumerate(samples):
        print('{}...'.format(sid), end='')
        # calculate the number of  overlapping bases
        total = calculate_total_bases(query_bed_file, ref_bed_files[sid])
        # process summary
        cmtx_list.append(process_results(summary_files[k], sid, total))
        print('ok.')
    # concate matrix
    return pd.concat(cmtx_list, ignore_index=True)

def fetch_fns(sid, bvcf_file):
    """extract FNs in a sample"""
    # read in benchmark vcf file
    bvcf = pd.read_table(bvcf_file, header=None, names=['Chromosome','Position','ID','Ref','Alt','QUAL','Filter','INFO','Format','Benchmark','Ours'], comment='#', low_memory=False)
    # FN or LowVAF+TP
    filt1 = bvcf['Benchmark'].str.contains('FN')
    filt2 = bvcf['Filter'].str.contains('LowVAF')
    filt3 = bvcf['Benchmark'].str.contains('TP')
    # add a 'sample' column
    bvcf['Sample'] = sid
    # get FN entries, select columns for output
    return bvcf[filt1 | (filt2 & filt3)][['Sample','Chromosome','Position','Ref','Alt','Filter','Format','Benchmark','Ours']]

def prepare_fns_table(samples, bench_vcf_files):
    """prepare 2nd excel sheet: FNs"""
    print('Preparing FN table:')
    fns_list = []
    for k,sid in enumerate(samples):
        print('{}...'.format(sid), end='')
        # extract FNs
        fns_list.append(fetch_fns(sid, bench_vcf_files[k]))
        print('ok.')
    # concate fns table
    return pd.concat(fns_list, ignore_index=True)

def fetch_fps(sid, bvcf_file):
    """extract FPs in a sample"""
    # read in benchmark vcf file
    bvcf = pd.read_table(bvcf_file, header=None, names=['Chromosome','Position','ID','Ref','Alt','QUAL','Filter','INFO','Format','Benchmark','Ours'], comment='#', low_memory=False)
    # FP and not LowVAF
    filt1 = bvcf['Ours'].str.contains('FP')
    filt2 = ~ bvcf['Filter'].str.contains('LowVAF')
    # add a 'sample' column
    bvcf['Sample'] = sid
    # get FP entries, select columns for output
    return bvcf[filt1 & filt2][['Sample','Chromosome','Position','Ref','Alt','Filter','Format','Benchmark','Ours']]

def prepare_fps_table(samples, bench_vcf_files):
    """prepare 3rd excel sheet: FPs"""
    print('Preparing FP table:')
    fps_list = []
    for k,sid in enumerate(samples):
        print('{}...'.format(sid), end='')
        # extract FNs
        fps_list.append(fetch_fps(sid, bench_vcf_files[k]))
        print('ok.')
    # concate fps table
    return pd.concat(fps_list, ignore_index=True)

def write_tables_to_excel(out_excel_file, samples, summary_files, query_bed_file, ref_bed_files, bench_vcf_files):
    """write tables to the excel file"""
    # overall stats
    table1 = prepare_overall_table(samples, summary_files, query_bed_file, ref_bed_files)
    # fns
    table2 = prepare_fns_table(samples, bench_vcf_files)
    # fps
    table3 = prepare_fps_table(samples, bench_vcf_files)
    # create a excel writer and write tables
    with pd.ExcelWriter(out_excel_file) as writer:
        table1.to_excel(writer, index=False, sheet_name='overall')
        table2.to_excel(writer, index=False, sheet_name='FNs')
        table3.to_excel(writer, index=False, sheet_name='FPs')

with open(log_file, 'w') as flog:
    with redirect_stdout(flog):
        write_tables_to_excel(out_excel_file, samples, summary_files, query_bed_file, ref_bed_files, bench_vcf_files)
        print('Complete!')

