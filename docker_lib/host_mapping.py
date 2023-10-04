#!/usr/bin/env python

# Host mapping script by J.J.Iszatt

import os
import sys
import pandas as pd
import configparser
import datetime

# Reading inputs
input = '/assemble/input'
output = '/assemble/output'
logs = "/assemble/output/phanatic_log.tsv"

# Configuring pipeline
config_file = "/assemble/output/config.ini"
config = configparser.ConfigParser()
if os.path.isfile(config_file):
    config.read(config_file)
else:
    config_file = "/assemble/config.ini"
    config.read(config_file)

image = config["phanatic"]["image"]
r1_ext = config["input"]["r1_ext"]
r2_ext = config["input"]["r2_ext"]
threads = int(config["SPAdes"]["threads"])

# Phanatic settings
try:
    enable_normalise = config.getboolean("pipeline", "normalise")
    enable_filter = config.getboolean("pipeline", "filter")
    enable_qc = config.getboolean("pipeline", "fastqc")
    enable_barcodes = config.getboolean("pipeline", "barcode")
    enable_clean = config.getboolean("pipeline", "clean_up")
except ValueError:
    sys.exit("Config file incorrectly set, pipeline values must be booleans")

################################################################################
#_______________________________________________________________________________

# Functions
def logfile(function, text, logfile):
    newline = "\n"
    date_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    report = f"{newline}[{image}]\t[{date_time}]\t[{function}]\t[{text}]"
    with open(logfile, "a") as file:
        file.write(report)

def mapping(genome, reads, outdir):
    
    # Creating temp file and moving directories
    os.makedirs(outdir)
    os.system(f"cp {genome} {outdir}")
    
    # Reusing genome variable
    genome = os.path.join(outdir, os.path.basename(genome)) 
    
    # Indexing and alignment
    aligned_sam = os.path.join(outdir, "aligned.sam")
    os.system(f"bwa index {genome}")
    os.system(f"bwa mem {genome} {reads} > {aligned_sam}")
    
    # Sorting sam file 
    aligned_bam = os.path.join(outdir, "aligned.bam")
    sorted = os.path.join(outdir, "sorted.bam")
    outfile = os.path.join(outdir, "results.txt")
    os.system(f"samtools view -bS {aligned_sam} > {aligned_bam} -t {threads}")
    os.system(f"samtools sort {aligned_bam} > {sorted} -t {threads}")
    os.system(f"samtools flagstat {sorted} > {outfile}")
    
    # Creating reads dir
    reads_dir = os.path.join(outdir, 'reads')
    os.makedirs(reads_dir)
    
    # Separating mapped and unmapped reads
    mapped = os.path.join(reads_dir, 'mapped.bam')
    unmapped = os.path.join(reads_dir, 'unmapped.bam')
    os.system(f"samtools view -F 4 -b {sorted} > {mapped} -t {threads}")
    os.system(f"samtools view -f 4 -b {sorted} > {unmapped} -t {threads}")
    
    # Generating pileup file
    pile = os.path.join(outdir, "pileup.txt")
    os.system(f"samtools mpileup -f {genome} {sorted} > {pile}")

################################################################################
#_______________________________________________________________________________

# CLASSES (UNUSED)
class Mapper(object):
    def __init__(self, host_genome, phage_genome, phage, qc_reads):
        self.host_genome = host_genome;
        self.phage_genome = phage_genome;
        self.phage = phage;
        self.qc1 = qc_reads;

################################################################################
#_______________________________________________________________________________


# Read mapping
logfile("Phanatic host mapping", "collating bacterial hosts", logs)
mapping_file = "/assemble/output/host_mapping.csv"
if os.path.isfile(mapping_file):
    df = pd.read_csv(mapping_file)
    df.dropna(inplace=True)
else:
    logfile("Phanatic host mapping", "ERROR, mapping file not found", logs)
    sys.exit("ERROR, mapping file not found")

# Checking dataframe
if not 'host' in df.columns:
    sys.exit("DF ERROR")
if len(list(df['host'])) == 0:
    sys.exit("DF ERROR")
if not 'read_1' in df.columns:
    sys.exit("DF ERROR")
if not 'read_2' in df.columns:
    sys.exit("DF ERROR")

# Checking index file
for index, row in df.iterrows():
    host = os.path.join(input, row['host'])
    read_1 = os.path.join(input, row['read_1'])
    read_2 = os.path.join(input, row['read_2'])
    cut = len(r1_ext)
    
    if not os.path.exists(host):
        print("Mapping file error: host file does not exist")
        continue
    elif not os.path.exists(read_1):
        print("Mapping file error: r1 does not exist")
        continue
    elif not os.path.exists(read_2):
        print("Mapping file error: r2 does not exist")
        continue
    else:
        name = os.path.basename(read_1[:-cut])
    
    # Finding QC reads 
    print(name)
    qc = os.path.join(output, "deduped", f"{name}.fastq")
    if os.path.exists(qc):
        logfile("Successful host pairing", f"{row['host']} to {name}", logs)

    # Mapping reads to host 
    outdir = os.path.join(output, "host_read_mapping", (row['host'].replace(".fasta", "")+"_"+name))
    mapping(host, qc, outdir)
    logfile("Host mapping complete", host, logs)

    # Obtaining list of phages
    phages = []
    format_dir = os.path.join(output, "format_dir")
    outdir = os.path.join(output, "phage_read_mapping")
    
    # Mapping reads to phage
    for file in os.listdir(format_dir):
        if name in file:
            path = os.path.join(format_dir, file)
            outdir = os.path.join(output, "phage_read_mapping", file)
            
            # Mapping reads
            logfile(f"Mapping reads to phage", f"{os.path.basename(qc)} to {file}", logs)
            mapping(path, qc, outdir)
            
            

logfile(f"Finished host_mapping.py", "---------------------------------", logs)
