#!/usr/bin/env python

# Coordinator script by J.J.Iszatt

import os
import sys
import datetime
import subprocess
import configparser
from Bio import SeqIO
import csv
import functions as ji

# Reading input directory
input = '/assemble/input'
output = '/assemble/output'

# Configuring pipeline
config_file = "/assemble/output/config.ini"
config = configparser.ConfigParser()
if os.path.isfile(config_file):
    config.read(config_file)
else:
    config_file = "/assemble/config.ini"
    config.read(config_file)

# Phanatic settings
enable_normalise = config["pipeline"]["normalise"]
enable_filter = config["pipeline"]["filter"]
enable_qc = config["pipeline"]["fastqc"]
enable_extract = config["pipeline"]["extract"]

# Reading input files
pairs = ji.find_read_pairs(input)
print(f"Paired read files: {len(pairs)}")

# Setting directories
trim_dir = os.path.join(output, "trimmed")
dedupe_dir = os.path.join(output, "deduped")
merged_dir = os.path.join(output, "merged")
norm_dir = os.path.join(output, "normalised")
spades_dir = os.path.join(output, "spades")
raw_genome_dir = os.path.join(output, "raw_assembly")
checkv_dir = os.path.join(output, "checkv")
qualimap_dir = os.path.join(output, "qualimap")

# Phanatic run 
for pair in pairs:
    
    # Reads processing:
    trim = ji.PE_trim(pair, trim_dir)
    
    if ji.check_filepath(trim):
        deduped = ji.remove_duplicate_reads(trim, dedupe_dir, pair.name)
    else:
        continue
        
    if ji.check_filepath(deduped):
        merged, unmerged = ji.merge_reads(deduped, merged_dir, pair.name)
    else:
        continue

    if ji.check_filepath(merged):
        normalised = ji.normalise_reads(merged, norm_dir, pair.name)
    else:
        continue

    # Assembly
    if ji.check_filepath(normalised):
        assembly = ji.PE_assembly(normalised, unmerged, spades_dir, pair.name)
    else:
        continue
        
    if ji.check_filepath(assembly):
        filtered = ji.filter_genome(assembly, raw_genome_dir, pair.name)
    else:
        continue
    
    # CheckV
    if ji.check_filepath(filtered):
        checkv = ji.checkv(filtered, checkv_dir, pair.name)
    else:
        continue
        
    # Extraction
    if ji.check_filepath(checkv):
        
    
    # Quality checks
    
    # Barcoding
    
    
    
        
    
    


    
    

    
    