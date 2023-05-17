
import os
import sys
import datetime
import configparser
import subprocess
import csv

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


###_______________________________________________________________________________________

## CLASSES

class Pair(object):
    def __init__(self, name, read_1, read_2, filter=None):
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
        
    logfile("Finish", f"Read pairs = {len(read_pairs)}", logs)
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
        logfile("Trim", f"{read_pair.name}", logs)
        return outfile
    except subprocess.CalledProcessError:
        logfile("Trim fail", f"{read_pair.name}", logs)

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
        logfile("Dedupe", name, logs)
        return outfile
    except subprocess.CalledProcessError:
        logfile("Dedupe fail", name, logs)
        
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
        logfile("Merge", name, logs)
        return outfile_merged, outfile_unmerged
    except subprocess.CalledProcessError:
        logfile("Merge fail", name, logs)
    
    
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
        logfile("Normalise", name, logs)
        return outfile
    except subprocess.CalledProcessError:
        logfile("Normalise", name, logs)

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
        logfile("Assembly", name, logs)
        return f"{outdir}/{name}/contigs.fasta"
    except subprocess.CalledProcessError:
        logfile("Assembly failed", name, logs)

def filter_genome(infile, outdir, name):
    
    outfile = f"{outdir}/{name}.fasta"

    command = [
        "reformat.sh",
        f"in={infile}",
        f"out={outfile}",
        f"minlength={filter_length}"
    ]
    try:
        subprocess.run(command, check=True)
        logfile("Filter", name, logs)
        return outfile
    except subprocess.CalledProcessError:
        logfile("Filter failed", name, logs)

def checkv(infile, outdir, name):
    
    outfile = f"{outdir}/{name}/complete_genomes.tsv"
    
    command = [
        "checkv", "end_to_end",
        f"{infile}",
        f"{outdir}/{name}"
    ]
    try:
        subprocess.run(command, check=True)
        logfile("CheckV", name, logs)
        return outfile
    except subprocess.CalledProcessError:
        logfile("CheckV failed", name, logs)

