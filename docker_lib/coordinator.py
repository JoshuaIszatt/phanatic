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
    enable_mapping = config.getboolean("pipeline", "mapping")
    enable_reassembly = config.getboolean("pipeline", "re_assembly")
except ValueError:
    sys.exit("Config file incorrectly set, pipeline values must be booleans")

# Reading input files
pairs = ji.find_read_pairs(input)
print(f"Paired read files: {len(pairs)}")

# Setting directories
trim_dir = os.path.join(output, "trimmed")
dedupe_dir = os.path.join(output, "deduped")
norm_dir = os.path.join(output, "normalised")
spades_dir = os.path.join(output, "initial_assembly")
filtered_dir = os.path.join(output, "filtered_contigs")
checkv_dir = os.path.join(output, "checkv")
extraction_dir = os.path.join(output, "genome_extractions")
format_dir = os.path.join(output, "phage_genomes")
barcode_dir = os.path.join(output, "barcode_phage")

if enable_mapping:
    ji.logfile("pipeline options", "mapping enabled", logs)
    mapped = os.path.join(output, "mapping_QC_to_phage")
    mapped2 = os.path.join(output, "mapping_Norm_to_phage")

if enable_reassembly:
    ji.logfile("pipeline options", "mapped reassembly enabled", logs)
    mapped_assembly = os.path.join(output, "mapping_reassembly")

if enable_qc:
    ji.logfile("pipeline options", "QC enabled", logs)
    qc_dir = os.path.join(output, "reads_quality")
    os.makedirs(qc_dir)

# Host mapping file
enable_host_mapping = False
host_mapping_file = "/assemble/output/host_mapping.csv"
if os.path.exists(host_mapping_file):
    ji.logfile("pipeline options", "host mapping enabled", logs)
    host_mapping_dir = os.path.join(output, "mapping_QC_to_host")
    phage_host_mapping_dir = os.path.join(output, "mapping_phageQC_to_host")
    enable_host_mapping = True

# Phanatic run 
for pair in pairs:
    
    # Trimming reads
    trim = ji.PE_trim(pair, trim_dir)
    
    # Removing duplicates
    if ji.check_filepath(trim):
        deduped = ji.remove_duplicate_reads(trim, dedupe_dir, pair.name)
    else:
        continue

    # Removing trimmed reads
    try:
        os.remove(trim)
    except Exception as e:
        ji.logfile(f"Error: {trim} could not be removed", f"{e}", logs)

    # Normalising reads
    if ji.check_filepath(deduped):
        if enable_normalise:
            normalised = ji.normalise_reads(deduped, norm_dir, pair.name)
            assemble_reads = normalised
        else:
            assemble_reads = deduped

    # Assembly
    if ji.check_filepath(assemble_reads):
        assembly = ji.PE_assembly(assemble_reads, spades_dir, pair.name)
    else:
        continue
    
    # Checking for empty assembly
    if os.path.getsize(assembly) == 0:
        ji.logfile("ERROR: contigs file empty", f"{pair.name}: check SPAdes log", logs)
        continue
    
    # Filtering assembly
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
    
    # Mapping reads to phage contigs
    if enable_mapping:
        
        # QC
        ji.logfile("QC read mapping to phage contigs", f"{pair.name}", logs)
        ji.map_reads(filtered, deduped, mapped, pair.name)

        # Normalised
        if enable_normalise:
            ji.logfile("Normalised / subsampled read mapping to phage contigs", f"{pair.name}", logs)
            ji.map_reads(filtered, assemble_reads, mapped2, pair.name)
    
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
    
    # Making summary file
    sample_file = os.path.join(output, 'sample_summary.csv')
    if not os.path.exists(sample_file):
        ji.create_csv(sample_file, "sample,genomes,sample_status")
    
    if len(headers) == 0:
        ji.logfile("Sample failed", pair.name, logs)
        ji.append_csv(sample_file, f"{pair.name},{len(headers)},failed")
        continue
    elif len(headers) == 1:
        ji.logfile("Clean sample", pair.name, logs)
        ji.append_csv(sample_file, f"{pair.name},{len(headers)},clean")
    elif len(headers) > 1:
        ji.logfile("Potential contamination", pair.name, logs)
        ji.append_csv(sample_file, f"{pair.name},{len(headers)},contaminated")
    
    # Making extraction dir
    if not os.path.exists(extraction_dir):
        os.makedirs(extraction_dir)

    # Looping through contigs
    contig_file = os.path.join(output, 'contig_summary.csv')
    if not os.path.exists(contig_file):
        ji.create_csv(contig_file, "sample,contig_name,norm_coverage,cov_status,phage_QC_mapped_(%),mapped_status")
    genomes = []
    for header in headers:
        
        # Mapping filter: 80% of target coverage
        ji.logfile("Coverage filtering", f"Scanning: {header}", logs)
        covstat = os.path.join(mapped2, pair.name, 'covstats.tsv')
        scafstat = os.path.join(mapped, pair.name, 'scafstats.tsv')
        
        # Coverage filter
        coverage, coverage_target = ji.covstat_filter(header, covstat)
        perc_mapped, perc_target = ji.scafstat_filter(header, scafstat)
        
        # Assessing None values
        if coverage is None:
            ji.logfile("Coverage check failed", f"Check {header} coverage", logs)
        if perc_mapped is None:
            ji.logfile("Mapping percentage check failed", f"Check {header} QC mapping", logs)

        # Logging file data
        if coverage >= coverage_target:
            cov_pass = "PASS"
            ji.logfile("Coverage check PASSED", header, logs)
        else:
            cov_pass = "WARNING"
            ji.logfile("Coverage check: WARNING, low norm cov detected", header, logs)

        # Logging file data
        if perc_mapped >= perc_target:
            perc_pass = "PASS"
            ji.logfile("Percentage mapping check PASSED", header, logs)
        else:
            perc_pass = "FAIL"
            ji.logfile("Percentage mapping check: FAILED", header, logs)
        
        # Recording data
        ji.append_csv(contig_file, f"{pair.name},{header},{coverage},{cov_pass},{perc_mapped},{perc_pass}")
        
        # PASS / FAIL
        if perc_pass == 'FAIL':
            continue
        
        # Genome extraction from filtered contigs
        extraction = ji.extract_genome(filtered, 
                                        header, 
                                        extraction_dir, 
                                        pair.name)
        genomes.append(extraction)

    ji.logfile("Genomes extracted", f"{pair.name}: {len(genomes)}", logs)

    # Host mapping process
    if enable_host_mapping:
        
        # Mapping QC reads to host genome (For contamination check and signs of transduction by % mapped reads)
        host = None
        host = ji.host_csv_scan(host_mapping_file, pair.read_1, pair.read_2)
        if host is not None:
            bacteria_name = os.path.basename(host).replace(".fasta", "")
            ji.logfile("Mapping QC reads to host", f"{pair.name} QC reads mapped to {bacteria_name}", logs)
            ji.map_reads(host, deduped, host_mapping_dir, f"{pair.name}_{bacteria_name}")

    # Looping through checkv genomes
    formatted_genomes = []
    for genome in genomes:
        name = os.path.basename(genome).replace(".fasta", "")
        format_genome = ji.format_genome(genome, format_dir, name)
        formatted_genomes.append(format_genome)

        # Separating reads for mapped reassembly process
        if enable_reassembly:
            qc_map, qc_unmap, outdir = ji.separate_reads(genome, deduped, mapped_assembly, name)
            
            # Assembling mapped and unmapped reads
            if not os.path.getsize(qc_map) == 0:
                mapped_contigs = ji.PE_assembly(qc_map, outdir, "spades_mapped")
            if not os.path.getsize(qc_unmap) == 0:
                unmapped_contigs = ji.PE_assembly(qc_unmap, outdir, "spades_unmapped")    

            # Mapping QC phage_mapped reads to host genome (For transduction check of target phage, looking for regions >5000bp long)
            #ji.logfile(f"Mapping phage reads to host", f"{name} reads mapped to {bacteria_name}", logs)
            #ji.map_reads(host, qc_map, phage_host_mapping_dir, name)
            
 
    # Quality checks
    if enable_qc:
        ji.fastqc(assemble_reads, qc_dir) 
    
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
        if not os.path.exists(index):
            ji.create_csv(index, "sample_name,phage_ID")
        ji.append_csv(index, f"{file},{new_tag}.fasta")
        
    # Formatting
    for file in os.listdir(barcode_dir):
        filepath = os.path.join(barcode_dir, file)
        genome = ji.format_genome(filepath, barcode_dir, file)
        
        if os.path.exists(genome):
            os.system(f"rm {filepath}")

# Phanatic finish
ji.logfile("Phanatic base assembly finished", "-----", logs)
os.system(f"chmod -R 777 {output}/*")


'''
Host mapping genomes:

    > Samples were considered free from contamination 
    if the total assembly size was close to the range 
    of genome sizes within that particular virus family 
    and if less than 5% of reads mapped to host reference 
    genomes.
    
    > The recommended category for viruses used in animal models for vaccine 
    development, and by extension phage therapy, is “Finished”. 
    Finished status is defined as a single consensus sequence 
    representing 100% of the genome with all open reading frames 
    (ORFs) identified and population diversity, or lack of population 
    diversity as an indicator of purity, of the sequence 
    verified via deep coverage.
    
    > Assemblies were considered validated if >90% of reads mapped back 
    to the phage genome
    
    > Additionally, in order for genomes to proceed from this checkpoint, 
    average whole genome coverage and lowest coverage were at least 
    100x for complete genomes and ~400x or finished genomes
    
    > Finished viral genomes have a complete sequence plus a minimum of 
    400-1000x coverage depth to resolve population-level variations
    
    - Philipson et al, (2018), Characterizing Phage Genomes for Therapeutic Applications
    
    

So what metrics do we need for this:

    More than 90% reads mapping to phage contig scafstats.txt from QC to phage mapped
    
    Average whole genome coverage + lowest coverage for any base: 
        = 100X for 'complete'
        = 400X for 'finished'

    Less than 5% reads mapping to host contig scafstats.txt from QC to host mapped (if available)
    
    
THE AIM:
    Finished: This final category represents a special instance in
which, in addition to having a completed consensus genome sequence, 
there has been a population-level characterization of genomic diversity. 
Typically this requires ~400 to 1,000 coverage (see below). This provides 
the most complete picture of a viral population; however, this 
designation will apply only for a single stock. 

Additional characterizations will be necessary for future
passages

- Ladner 2014 Standards for Sequencing Viral Genomes in the Era of HighThroughput Sequencing
'''
