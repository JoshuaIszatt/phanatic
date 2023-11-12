#!/usr/bin/env python
#

import os
import sys
import matplotlib.pyplot as plt
import pandas as pd

### Functions

def generate_coverage_graph(covstats, basecov, outdir):

    # Read covstats file
    headers = [
        "ID",
        "Avg_fold",
        "Length",
        "Ref_GC",
        "Covered_percent",
        "Covered_bases",
        "Plus_reads",
        "Minus_reads",
        "Read_GC",
        "Median_fold",
        "Std_Dev"
    ]
    df = pd.read_csv(covstats, sep='\t', comment='#', names=headers)

    # Separating data from covstats
    finished = df[df['Avg_fold'] >= 400]
    complete = df[df['Avg_fold'] >= 100]
    complete = complete[complete['Avg_fold'] <= 400]
    check = df[df['Avg_fold'] >= 50]
    check = check[check['Avg_fold'] <= 100]

    if not len(check) == 0:
        print("ERROR, manual curation required")
    elif len(finished) > 1:
        print("ERROR, >1 finished genome")
    elif len(complete) > 1:
        print("ERROR, >1 complete genome")
    elif len(complete) + len(finished) > 1:
        print("ERROR, A finished and complete genome are present")
    elif len(finished) == 1:
        contig = finished['ID'][0]
        finished = True
    elif len(complete) == 1:
        contig = complete['ID'][0]
        finished = False

    # Checking what we're working with
    if finished:
        print("Genome has > 400 X avg read depth")
        cutoff = 400
    else:
        print("Genome has 100-400 X avg read depth")
        cutoff = 100

    # Read data from the basecov file
    headers = ["ID", "Pos", "Coverage"]
    df = pd.read_csv(basecov, sep='\t', comment='#', names=headers)

    # Separating data from basecov
    coverage = df[df['ID'] == contig]

    # Plot basecoverage
    x_values = coverage['Pos']
    y_values = coverage['Coverage']

    # Check coverage
    below_cutoff_count = sum(1 for y in y_values if y < cutoff)
    if any(y < cutoff for y in y_values):
        print(f"Warning: {below_cutoff_count} coverage values have dipped below {cutoff}")

    # Create a plot for coverage data
    plt.figure(figsize=(15, 8))
    plt.plot(x_values, 
            y_values, 
            marker = ',', 
            markersize = 0.1,
            linestyle = '-', 
            color='b')
    plt.title(f"Per base coverage for {contig}")
    plt.xlabel("Position")
    plt.ylabel("Coverage")
    plt.grid(True)

    # Control line
    plt.axhline(y=400, color='red', linestyle=':')
    plt.axhline(y=100, color='red', linestyle='-')

    # Save the plot as an image file (e.g., PNG)
    outfile = os.path.join(outdir, f"{contig}.png")
    plt.savefig(outfile, dpi = 300)

########################################################################

# Directories
qc_phage = '/assemble/output/mapping_QC_to_phage'
norm_phage = '/assemble/output/mapping_Norm_to_phage'
qc_host = '/assemble/output/mapping_QC_to_host'
dirs = [qc_phage, norm_phage, qc_host]

# Running script
for directory in dirs:
    if os.path.exists(directory):
        for file in os.listdir(directory):
            print(file)
            # Setting function inputs
            dirpath = os.path.join(directory, file)
            cov = os.path.join(dirpath, "covstats.tsv")
            base = os.path.join(dirpath, "basecov.tsv")

            # Running graph
            try:
                generate_coverage_graph(cov, base, dirpath)
            except Exception as e:
                print(e)
    else:
        print("error")

######
'''
The below section may or may not be necessary, it concatenates data for all contigs > 1000bp including:
    > Mapping statistics from cov and scaf files from QC reads
    > checkv statistics from quality summary, completeness, and contamination files
'''
######

def filescan(dir, filetype):
    dfs = []
    for file in os.listdir(dir):
        path = os.path.join(dir, file, filetype)
        if path.endswith('.tsv'):
            df = pd.read_csv(path, sep='\t')
        elif path.endswith('.csv'):
            df = pd.read_csv(path)
        else:
            print(f"Error reading {path}")
            continue
        
        # Adding to list
        df.insert(0, 'sample', file)
        dfs.append(df)
        
    df = pd.concat(dfs)
    return df
    
# Preparing function parameters
outdir = '/assemble/output/'

# CheckV files
checkv_dir = os.path.join(outdir, "checkv")
qual = "quality_summary.tsv"
comp = "completeness.tsv"
cont = "contamination.tsv"

try:
    df = filescan(checkv_dir, qual)
    cols = ['provirus', 'proviral_length', 'viral_genes', 'host_genes', 'provirus', 'proviral_length', 'kmer_freq']
    df = df.drop(columns=cols)
    df2 = filescan(checkv_dir, comp)
    df3 = filescan(checkv_dir, cont)
    
    # Merging
    checkv = df.merge(df2, on=['sample','contig_id', 'contig_length']).merge(df3, on=['sample','contig_id', 'contig_length'])
    
except Exception as e:
    print(f"ERROR {e}")

# Coverage files
qc_dir = os.path.join(outdir, 'mapping_QC_to_phage')
covstat = 'covstats.tsv'
scafstat = 'scafstats.tsv'

try:
    df = filescan(qc_dir, covstat)
    df2 = filescan(qc_dir, scafstat)
    
    # Merging
    df.rename(columns={'#ID' : 'contig_id'}, inplace=True)
    df2.rename(columns={'#name' : 'contig_id'}, inplace=True)
    coverage = df.merge(df2, on=['sample', 'contig_id'])
    
except Exception as e:
    print(f"ERROR {e}")

# Final merge
merge = coverage.merge(checkv, on=['sample', 'contig_id'])

# Saving file 
outfile = os.path.join(outdir, 'raw_data.csv')
merge.to_csv(outfile, index=False)

######
'''

The below summary files need to answer the following questions:
    > Do the mapped QC reads produce a genome of exact same size when re assembled? (Present as yes/no)
    > Do the mapped QC reads assembly contain any other contigs >1000bp ? (Present as number of contigs >1000bp)
    > Do the unmapped reads assemble into anything above 1000bp ? (Present as number of contigs >1000bp)

mapped_assembly.csv should contain the mapped putative genome name, 
the name of the first NODE in the mapped reassembly, and the number of genomes >1000bp

unmapped_assembly.csv should contain the mapped putative genome name, and the number of
genomes >1000bp within the assembly

^^^ Potentially all in one file??? 

'''
######

# Building combined summary file
try:
    # Reading merge files
    contig_sum = os.path.join(outdir, 'sample_summary.csv')
    sample_sum = os.path.join(outdir, 'contig_summary.csv')
    df = pd.read_csv(contig_sum)
    df2 = pd.read_csv(sample_sum)

    # Merging + sorting
    merge = df.merge(df2, on='sample', how='outer')
    merge.sort_values(by='phage_QC_mapped_(%)', ascending=False, inplace=True)

    # Saving
    outfile = os.path.join(outdir, 'combined_summary.csv')
    merge.to_csv(outfile, index=False)

except Exception as e:
    print(f"ERROR: {e}")

# Building reassembly_summary.csv
# Assessing transduction using basecov for host
