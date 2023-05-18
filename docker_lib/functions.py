
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

def create_csv(filename, data):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)

def append_csv(filename, data):
    with open(filename, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(data)

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
        "qhdist=1",
        "qtrim=10",
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

def PE_assembly(infile_1, infile_2, outdir, name):
    command = [
        "spades.py",
        "-t", f"{threads}",
        "-m", f"{memory_gb}",
        "--only-assembler",
        "--careful",
        "-k", "55,77,99,127",
        "-o", f"{outdir}/{name}",
        "--merged", f"{infile_1}",
        "-s", f"{infile_2}"
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

def coverage_calculation(genome, reads, outdir, name):
    cov_out = f"{outdir}/{name}.tsv"
    command = [
        "bbmap.sh", 
        f"-Xmx{memory}",
        f"ref={genome}",
        f"in={reads}",
        f"covstats={cov_out}"
    ]
    try:
        subprocess.run(command, check=True)
        logfile("Coverage", f"{name}: success", logs)
    except subprocess.CalledProcessError:
        logfile("Coverage", f"{name}: failed", logs)
    
def fastqc(reads, outdir):
    command = [
        "fastqc",
        f"{reads}",
        f"-o {outdir}"
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
