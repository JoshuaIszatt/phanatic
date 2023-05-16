#!/usr/bin/env python
# Python script to run the assemble docker container

import os
import sys
import subprocess
import argparse
import random
import pandas as pd

# Did you know prompts
prompts = [
    "Phanatic currently only supports short, paired end, illumina reads (default 'PE_illumina_150').",
    "I (J. Iszatt) made Phanatic to assemble S. aureus phage genomes as part of my PhD project.",
    "This is a short read assembler primarily designed for bacteriophage",
    "Use the --nanopore flag to include long read sequencing in your assembly",
    "My favourite bacteriophage is a Silviavirus named Koomba-kaat_1"
]
random_prompt = random.choice(prompts)

# Creating function to check directory path
def valid_dir(dir_path):
    if not os.path.isdir(dir_path):
        raise argparse.ArgumentTypeError(
            f"{dir_path} is not a valid directory path")
    if not os.access(dir_path, os.R_OK):
        raise argparse.ArgumentTypeError(
            f"{dir_path} is not a readable directory")
    return dir_path

def valid_file(file_path):
    if not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(
            f"{file_path} is not a valid file path")
    if not os.access(file_path, os.R_OK):
        raise argparse.ArgumentTypeError(
            f"{file_path} is not a readable file")
    return file_path

# Parsing arguments
image = 'iszatt/phanatic:2.2.0'
parser = argparse.ArgumentParser(description=f"Easy short read assembly. Joshua J Iszatt: https://github.com/JoshuaIszatt")

# Input/output options
parser.add_argument('-i', '--input', type=valid_dir, help='Input reads files')
parser.add_argument('-o', '--output', type=valid_dir, help='Direct output to this location')
parser.add_argument('-r', '--reads', type=str, choices=['PE_illumina_150'], default='PE_illumina_150', help='Pipeline options')
parser.add_argument('-c', '--config', type=valid_file, help='Use config file to customise assembly')
#parser.add_argument('--nanopore', action="store_true", help='Detect and utilise nanopore reads')
parser.add_argument('--manual', action="store_true", help='Enter container interactively')
args = parser.parse_args()

# Obtaining absolute paths
input_path = os.path.abspath(args.input)
output_path = os.path.abspath(args.output)

# Printing command variables
print(
    f"Program run: {image}",
    f"Input path: {input_path}",
    f"Output path: {output_path}",
    f"Reads type: {args.reads}",
    ">>>\n",
    f"Did you know:",
    f"{random_prompt}",
    ">>>\n",
    sep='\n'
    )

# Copying config file to output dir
if args.config:
    print(f"Using {args.config} file \n")
    os.system(f"cp {args.config} {args.output}/config.ini")

# Running docker
if args.manual: 
    os.system(f"docker exec -it \
        $(docker run -d \
        -v {input_path}:/assemble/input \
        -v {output_path}:/assemble/output \
        {image} sleep 1d) bash")
else:
    command = ["docker run -d -v %s:/assemble/input -v %s:/assemble/output %s /assemble/bin/assemble.sh" %
            (input_path, output_path, image)]
    result = subprocess.Popen(command, shell=True)
    print(command)
