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
filtered_dir = os.path.join(output, "filtered_contigs")
checkv_dir = os.path.join(output, "checkv")
extraction_dir = os.path.join(output, "genome_extractions")


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
        filtered = ji.filter_genome(assembly, filtered_dir, pair.name)
    else:
        continue
    
    # CheckV
    if ji.check_filepath(filtered):
        checkv = ji.checkv(filtered, checkv_dir, pair.name)
    else:
        continue

    # Extractions
    complete_genomes = os.path.join(checkv, "complete_genomes.tsv")
    quality_summary = os.path.join(checkv, "quality_summary.tsv")
    
    if ji.check_filepath(complete_genomes):
        complete = ji.find_complete_genomes(complete_genomes, pair.name)
    else:
        continue
    
    if ji.check_filepath(quality_summary):
        hq = ji.find_hq_genomes(quality_summary, pair.name)
    else:
        continue
    
    if len(complete) + len(hq) == 0:
        continue
    else:
        extracted_genomes = ji.extract_genomes(complete, hq, extraction_dir)
    
    
    # Quality checks

    # Coverage calculation
    
    # Barcoding
    