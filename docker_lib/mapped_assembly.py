#!/usr/bin/env python

# Mapped assembly script by J.J.Iszatt

import os
import sys
import datetime
import subprocess
import configparser
from Bio import SeqIO
import csv
import functions as ji

# Reading inputs
input = '/assemble/input'
output = '/assemble/output'
logs = "/assemble/output/phanatic_log.tsv"

# Setting directory locations
phage_dir = os.path.join(output, "phage_read_mapping")

# Configuring pipeline
config_file = "/assemble/output/config.ini"
config = configparser.ConfigParser()
if os.path.isfile(config_file):
    config.read(config_file)
else:
    config_file = "/assemble/config.ini"
    config.read(config_file)

# Settings
threads = int(config["SPAdes"]["threads"])

# Read mapping
ji.logfile("Phanatic mapped assembly", "collating bacterial hosts", logs)
mapping_file = "/assemble/output/host_mapping.csv"
if not os.path.isfile(mapping_file):
    ji.logfile("Phanatic mapped assembly", "ERROR, mapping file not found", logs)
    sys.exit("ERROR, mapping file not found")

# Phanatic settings
try:
    enable_mapped = config.getboolean("host_mapping_pipeline", "mapped_assembly")
    enable_unmapped = config.getboolean("host_mapping_pipeline", "unmapped_assembly")
except ValueError:
    sys.exit("Config file incorrectly set, pipeline values must be booleans")

# Mapped assembly (To check contig size against original, indicating no interference from contaminating reads)
if enable_mapped:
    ji.logfile("Mapped assembly enabled", "---", logs)

    # Finding mapped reads
    for file in os.listdir(phage_dir):
        reads_path = os.path.join(phage_dir, file, "reads")
        mapped = os.path.join(reads_path, "mapped.bam")
        unmapped = os.path.join(reads_path, "unmapped.bam")
        
        # Converting mapped reads
        if os.path.isfile(mapped):
            mapped_reads = os.path.join(reads_path, "mapped.fastq")
            os.system(f"samtools fastq -@ {threads} {mapped} > {mapped_reads}")
        
        # Converting unmapped reads
        if os.path.isfile(unmapped):
            unmapped_reads = os.path.join(reads_path, "unmapped.fastq")
            os.system(f"samtools fastq -@ {threads} {unmapped} > {unmapped_reads}")
    
    

# Unmapped assembly (To query any contigs against inphared that may indicate prophage induction)
if enable_unmapped:
    ji.logfile("Unmapped assembly enabled", "---", logs)

