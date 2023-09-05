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

# Phanatic settings
try:
    enable_normalise = config.getboolean("pipeline", "normalise")
    enable_filter = config.getboolean("pipeline", "filter")
    enable_qc = config.getboolean("pipeline", "fastqc")
    enable_barcodes = config.getboolean("pipeline", "barcode")
    enable_clean = config.getboolean("pipeline", "clean_up")
except ValueError:
    sys.exit("Config file incorrectly set, pipeline values must be booleans")

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
format_dir = os.path.join(output, "format_dir")
coverage_dir = os.path.join(output, "coverage_dir")
barcode_dir = os.path.join(output, "barcode_phage")

if enable_qc:
    ji.logfile("pipeline options", "QC enabled", logs)
    qc_dir = os.path.join(output, "reads_quality")
    os.makedirs(qc_dir)

# Phanatic run 
for pair in pairs:
    
    # Trimming reads
    trim = ji.PE_trim(pair, trim_dir)
    
    # Removing duplicates
    if ji.check_filepath(trim):
        deduped = ji.remove_duplicate_reads(trim, dedupe_dir, pair.name)
    else:
        continue
    
    # Merging reads
    if ji.check_filepath(deduped):
        merged, unmerged = ji.merge_reads(deduped, merged_dir, pair.name)
    else:
        continue

    # Normalising reads
    if enable_normalise:
        if ji.check_filepath(merged):
            normalised = ji.normalise_reads(merged, norm_dir, pair.name)
        else:
            continue
        assemble_reads = normalised
    else:
        assemble_reads = merged
    
    # Assembly
    if ji.check_filepath(assemble_reads):
        assembly = ji.PE_assembly(assemble_reads, unmerged, spades_dir, pair.name)
    else:
        continue
    
    if os.path.getsize(assembly) == 0:
        ji.logfile("ERROR: contigs file empty", f"{pair.name}: check SPAdes log", logs)
        continue
    
    if ji.check_filepath(assembly):
        if enable_filter:
            filtered = ji.format_genome(assembly, filtered_dir, pair.name, filter=True)
        else:
            filtered = ji.format_genome(assembly, filtered_dir, pair.name)
    else:
        continue
    
    # CheckV
    if os.path.getsize(filtered) == 0:
        ji.logfile("ERROR: Filtered contigs empty", f"{pair.name}: check contigs file", logs)
        continue
    
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

    headers = complete+hq
    ji.logfile("Expected genomes", f"{pair.name}: {len(headers)}", logs)
    
    if len(headers) == 0:
        ji.logfile("Sample failed", pair.name, logs)
        continue
    elif len(headers) == 1:
        ji.logfile("Clean sample", pair.name, logs)
    elif len(headers) > 1:
        ji.logfile("Contamination detected", pair.name, logs)
    
    if not os.path.exists(extraction_dir):
        os.makedirs(extraction_dir)

    genomes = []
    for header in headers:
        extraction = ji.extract_genome(filtered, 
                                        header, 
                                        extraction_dir, 
                                        pair.name)
        genomes.append(extraction)

    ji.logfile("Genomes extracted", f"{pair.name}: {len(genomes)}", logs)

    formatted_genomes = []
    for genome in genomes:
        name = os.path.basename(genome).replace(".fasta", "")
        format_genome = ji.format_genome(genome, format_dir, name)
        formatted_genomes.append(format_genome)

    # Coverage calculation
    for genome in formatted_genomes:
        name = os.path.basename(genome).replace(".fasta", "")
        ji.coverage_calculation(genome, assemble_reads, coverage_dir, name)
    
    # Quality checks
    if enable_qc:
        if enable_normalise:
            ji.fastqc(normalised, qc_dir)
        else:
            ji.fastqc(merged, qc_dir)  
    
    # Sample finish
    ji.logfile("Sample run complete", pair.name, logs)

# Barcoding
if enable_barcodes:
    ji.logfile("Barcoding", "-----", logs)
    os.makedirs(barcode_dir)
    tags = []
    
    # Tagging
    for file in os.listdir(format_dir):
        filepath = os.path.join(format_dir, file)
        new_tag = ji.generate_unique_tag(tags)
        tags.append(new_tag)
        ji.barcode_phage(filepath, new_tag, barcode_dir)
        
        # Producing index file
        index = os.path.join(output, 'index.csv')
        ji.create_csv(index, "sample_name,phage_ID")
        ji.append_csv(index, f"{file},{new_tag}.fasta")
        
    # Formatting
    for file in os.listdir(barcode_dir):
        filepath = os.path.join(barcode_dir, file)
        genome = ji.format_genome(filepath, barcode_dir, file)
        
        if os.path.exists(genome):
            os.system(f"rm {filepath}")

# Cleaning up
if enable_clean:
    ji.logfile("Cleaning raw data", "-----", logs)
    remove = [trim_dir, 
              dedupe_dir, 
              merged_dir, 
              norm_dir, 
              spades_dir, 
              filtered_dir,
              extraction_dir
              ]
    for data in remove:
        try:
            os.system(f"rm -rf {data}")
        except:
            continue

# Phanatic finish
ji.logfile("Phanatic finished", "-----", logs)
os.system(f"chmod -R 777 {output}/*")
