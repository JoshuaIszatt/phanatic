
import os
import sys
import datetime
import configparser
import subprocess
import csv
import random
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

## CONFIGURATION

config_file = "/assemble/output/config.ini"
logs = "/assemble/output/phanatic_log.tsv"

config = configparser.ConfigParser()
if os.path.isfile(config_file):
    config.read(config_file)
else:
    config_file = "/assemble/config.ini"
    config.read(config_file)

image = config["phanatic"]["image"]
memory = config["system"]["RAM"]

r1_ext = config["input"]["r1_ext"]
r2_ext = config["input"]["r2_ext"]

read_length = int(config["trim"]["read_length"])
trim_length = int(config["trim"]["trim_length"])
minimum_length = config["trim"]["minimum_length"]
q_trim = config["trim"]["read_quality"]

minimum_insert = int(config["merge"]["minimum_insert"])
minimum_overlap = int(config["merge"]["minimum_overlap"])

target_coverage = int(config["normalise"]["target_coverage"])

threads = int(config["SPAdes"]["threads"])
memory_gb = int(config["SPAdes"]["memory_gb"])

filter_length = int(config["filter"]["filter_length"])

prefix = config["barcoding"]["prefix"]
barcode_length = int(config["barcoding"]["barcode_length"])

###_______________________________________________________________________________________

## CLASSES

class Pair(object):
    def __init__(self, name, read_1, read_2):
        self.name = name;
        self.read_1 = read_1;
        self.read_2 = read_2;
        
###_______________________________________________________________________________________


## FUNCTIONS

def logfile(function, text, logfile):
    newline = "\n"
    date_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    report = f"{newline}[{image}]\t[{date_time}]\t[{function}]\t[{text}]"
    with open(logfile, "a") as file:
        file.write(report)

def check_filepath(filepath, create=False):
    if os.path.exists(filepath):
        print(f"Found:{filepath}")
        return True
    else:
        if create:
            try:
                os.makedirs(filepath)
                print(f"Creating directory: {filepath}")
            except:
                print(f"Failed to create directory: {filepath}")
                sys.exit(1)
        else:
            print(f"Cannot find {filepath}")
            sys.exit(1)

def generate_unique_tag(existing_tags):
    characters = "abcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    while True:
        password = "".join(random.sample(characters, barcode_length))
        tag = f"{prefix}_{password}"
        if tag not in existing_tags:
            logfile("Barcoding", f"Tag gen: {tag}", logs)
            return tag

def create_csv(filename, headers):
    with open(filename, 'w', newline='') as csvfile:
        csvfile.write(headers+'\n')

def append_csv(filename, data):
    with open(filename, 'a', newline='') as csvfile:
        csvfile.write(data+'\n')

def find_read_pairs(input_dir):
    read_pairs = []
    for file in os.listdir(input_dir):
        filepath_1 = os.path.join(input_dir, file)
        
        if file.endswith(r1_ext):
            cut = len(r1_ext)
            name = file[:-cut]
            filepath_2 = os.path.join(input_dir, name + r2_ext)
            if os.path.isfile(filepath_2):
                logfile("Input file", f"Adding pair: {name}", logs)
                pair = Pair(name, filepath_1, filepath_2)
                read_pairs.append(pair)
            else:
                logfile("Input file", f"No second read file found for: {name}", logs)
        
    logfile("Pairing input files", f"Read pairs = {len(read_pairs)}", logs)
    return read_pairs

def PE_trim(read_pair, outdir):
    read_1 = read_pair.read_1
    read_2 = read_pair.read_2
    outfile = f"{outdir}/{read_pair.name}.fastq"
    tl = 0+trim_length
    tr = read_length-trim_length
    command = [
        "bbduk.sh",
        f"-Xmx{memory}",
        "tpe",
        "tbo",
        f"in1={read_1}",
        f"in2={read_2}",
        f"out={outfile}",
        f"ftl={tl}",
        f"ftr={tr}",
        f"minavgquality={q_trim}",
        f"minlength={minimum_length}"
    ]
    try:
        subprocess.run(command, check=True)
        logfile("Trimming", f"{read_pair.name}: success", logs)
        return outfile
    except subprocess.CalledProcessError:
        logfile("Trimming", f"{read_pair.name}: failed", logs)

def remove_duplicate_reads(infile, outdir, name):
    outfile = f"{outdir}/{name}.fastq"
    command = [
        "dedupe.sh",
        f"-Xmx{memory}",
        "ac=f",
        "s=5",
        "e=2",
        f"in={infile}",
        f"out={outfile}"
    ]
    try:
        subprocess.run(command, check=True)
        logfile("Dedupe", f"{name}: success", logs)
        return outfile
    except subprocess.CalledProcessError:
        logfile("Dedupe", f"{name}: failed", logs)
        
def merge_reads(infile, outdir, name):
    outfile_merged = f"{outdir}/{name}_merged.fastq"
    outfile_unmerged = f"{outdir}/{name}_unmerged.fastq"
    command = [
        "bbmerge.sh",
        f"-Xmx{memory}",
        f"mininsert={minimum_insert}",
        f"minoverlap={minimum_overlap}",
        f"in={infile}",
        f"out={outfile_merged}",
        f"outu={outfile_unmerged}"
    ]
    try:
        subprocess.run(command, check=True)
        logfile("Merge", f"{name}: success", logs)
        return outfile_merged, outfile_unmerged
    except subprocess.CalledProcessError:
        logfile("Merge fail", f"{name}: failed", logs)
    
    
def normalise_reads(infile, outdir, name):
    outfile = f"{outdir}/{name}.fastq"
    command = [
        "bbnorm.sh",
        f"-Xmx{memory}",
        "min=5",
        f"target={target_coverage}",
        f"in={infile}",
        f"out={outfile}"
    ]
    try:
        subprocess.run(command, check=True)
        logfile("Normalise", f"{name}: success", logs)
        return outfile
    except subprocess.CalledProcessError:
        logfile("Normalise", f"{name}: failed", logs)

def PE_assembly(reads, outdir, name):
    
    if memory_gb < 24:
        logfile("Warning", f"{memory_gb} GB is low memory for SPAdes", logs)
    
    command = [
        "spades.py",
        "-t", f"{threads}",
        "-m", f"{memory_gb}",
        "--only-assembler",
        "--careful",
        "-k", "55,77,99,127",
        "-o", f"{outdir}/{name}",
        "--12", f"{reads}"
    ]
    try:
        subprocess.run(command, check=True)
        logfile("Assembly", f"{name}: success", logs)
        return f"{outdir}/{name}/contigs.fasta"
    except subprocess.CalledProcessError:
        logfile("Assembly", f"{name}: failed", logs)

def format_genome(infile, outdir, name, filter=False):
    outfile = f"{outdir}/{name}.fasta"
    if filter:
        command = [
            "reformat.sh",
            f"in={infile}",
            f"out={outfile}",
            f"minlength={filter_length}"
        ]
        note = "Filtering"
    else:
        command = [
            "reformat.sh",
            f"in={infile}",
            f"out={outfile}",
            "minlength=0"
        ]
        note = "Formatting"
    try:
        subprocess.run(command, check=True)
        logfile(note, f"{name}: success", logs)
        return outfile
    except subprocess.CalledProcessError:
        logfile(note, f"{name}: failed", logs)

def checkv(infile, outdir, name):
    outfile = f"{outdir}/{name}"
    command = [
        "checkv", "end_to_end",
        f"{infile}",
        f"{outdir}/{name}"
    ]
    try:
        subprocess.run(command, check=True)
        logfile("CheckV", f"{name}: success", logs)
        return outfile
    except subprocess.CalledProcessError:
        logfile("CheckV", f"{name}: failed", logs)

def find_complete_genomes(checkv, name):
    genomes = []
    with open(checkv) as tsvfile:
        read = csv.reader(tsvfile, delimiter='\t')
        next(read)
        for row in read:
            genomes.append(row[0])
    genomes = tuple(genomes)
    logfile("Finding complete genomes", name, logs)
    return genomes

def find_hq_genomes(checkv, name):
    genomes = []
    with open(checkv) as tsvfile:
        read = csv.reader(tsvfile, delimiter='\t')
        next(read)
        for row in read:
            if 'High-quality' == row[7]:
                genomes.append(row[0])
    genomes = tuple(genomes)
    logfile("Finding high-quality genomes", name, logs)
    return genomes
    
def extract_genome(contigs, header, outdir, name):
    outfile = f"{outdir}/{name}_{header}.fasta"
    phage = f"{name}_{header}"
    handle = open(contigs)
    textfile = open(outfile, 'w')
    entries = list(SimpleFastaParser(handle)) 
    for name, seq in entries: 	
        if header in name: 		
            print(f">{phage}\n{seq}\n", file=textfile)
            logfile("Genome extracted", phage, logs)
            return outfile

def map_reads(genome, reads, outdir, name):
    
    # Output dir
    out = os.path.join(outdir, name)
    os.makedirs(out)
    
    # Output files
    covstats = os.path.join(out, "covstats.tsv")
    basecov = os.path.join(out, "basecov.tsv")
    scafstats = os.path.join(out, "scafstats.tsv")
    mapped = os.path.join(out, "mapped.fastq.gz")
    unmapped = os.path.join(out, "unmapped.fastq.gz")
    
    # Command    
    command = [
        "bbmap.sh", 
        f"-Xmx{memory}",
        f"ref={genome}",
        f"in={reads}",
        f"covstats={covstats}",
        f"basecov={basecov}",
        f"scafstats={scafstats}"
    ]
    try:
        subprocess.run(command, check=True)
        logfile("Read mapping", f"{name}: success", logs)
    except subprocess.CalledProcessError:
        logfile("Read mapping", f"{name}: failed", logs)

def separate_reads(genome, reads, outdir, name):
    
    # Output dir
    out = os.path.join(outdir, name)
    os.makedirs(out)
    
    # Output files
    mapped = os.path.join(out, "mapped.fastq.gz")
    unmapped = os.path.join(out, "unmapped.fastq.gz")
    
    # Command    
    command = [
        "bbmap.sh", 
        f"-Xmx{memory}",
        f"ref={genome}",
        f"in={reads}",
        f"outm={mapped}",
        f"outu={unmapped}"
    ]
    try:
        subprocess.run(command, check=True)
        logfile("Read extraction", f"{name}: success", logs)
    except subprocess.CalledProcessError:
        logfile("Read extraction", f"{name}: failed", logs)
    
    return mapped, unmapped, out

def fastqc(reads, outdir):
    command = [
        "fastqc",
        f"{reads}",
        "-o",
        f"{outdir}"
    ]
    try:
        subprocess.run(command, check=True)
        logfile("Reads QC", "success", logs)
    except subprocess.CalledProcessError:
        logfile("Reads QC", "failed", logs)

def barcode_phage(original, tag, outdir):
    original_name = os.path.basename(original)
    new = os.path.join(outdir, tag)
    os.system(f"cp {original} {new}")
    logfile("BARCODE", f"{original_name}:{tag}", logs)

def host_csv_scan(mapping_file, read_1, read_2):
    read_1 = os.path.basename(read_1)
    read_2 = os.path.basename(read_2)
    rows = []
    with open(mapping_file, newline='') as file:
        reader = csv.reader(file, delimiter=',', quotechar='|')
        for row in reader:
            if row[0] == 'host':
                continue
            elif os.path.basename(read_1) == row[1] and os.path.basename(read_2) == row[2]:
                rows.append(row)
            elif not os.path.exists(os.path.join("/assemble/input", row[0])):
                continue
            elif not os.path.exists(os.path.join("/assemble/input", row[1])):
                continue
            elif not os.path.exists(os.path.join("/assemble/input", row[2])):
                continue

    
    if len(rows) > 1:
        logfile("ERROR", f"Duplicate read host pairings found", logs)
        return None
    elif len(rows) == 1:
        logfile(f"Host Identified for {read_1}", f"{row[0]}", logs)
        return os.path.join("/assemble/input", row[0])
    elif len(rows) == 0:
        logfile("No host identified", f"---", logs)
        return None
    else:
        return None

def covstat_filter(header, covstat):
    cov = None
    cut = int(target_coverage * 0.8)
    # Obtaining covstat
    with open(covstat, newline='') as file:
        reader = csv.reader(file, delimiter='\t', quotechar='|')
        for row in reader:
            if row[0] == header:
                cov = float(row[1])
                break
    # Return value
    if cov is None:
        return None
    
    if cov:
        return cov, cut
        
def scafstat_filter(header, scafstat):
    perc_mapped = None
    # Obtaining unambig % mapped reads
    with open(scafstat, newline='') as file:
        reader = csv.reader(file, delimiter='\t', quotechar='|')
        for row in reader:
            if row[0] == header:
                perc_mapped = float(row[1])
    # Return value
    if perc_mapped is None:
        return None
    
    if perc_mapped:
        return perc_mapped, 90
