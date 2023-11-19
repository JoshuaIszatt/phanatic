#!/usr/bin/env python
#

import os
import hashlib
import pandas as pd

# Base input/output
input = '/assemble/input'
output = '/assemble/output'

# Hash function
def generate_sha256_hash(path):
    with open(path, 'rb') as file:
        file_contents = file.read()
    hash_object = hashlib.sha256(file_contents)
    hex = hash_object.hexdigest()
    return hex

# Initialising 
hash_files = []

# Find base files
log = os.path.join(output, 'phanatic_log.tsv')
raw = os.path.join(output, 'raw_data.csv')
summary = os.path.join(output, 'combined_summary.csv')

# Adding config and mapping files if they exist
config = os.path.join(output, 'config.ini')
mapping = os.path.join(output, 'host_mapping.csv')

# Checking files exist
files = [log, raw, summary, config, mapping]
for file in files:
    if os.path.exists(file):
        hash_files.append(file)

# Adding format_dir phage files
format_dir = os.path.join(output, 'phage_genomes')
for file in os.listdir(format_dir):
    path = os.path.join(format_dir, file)
    hash_files.append(path)

# Hashing files
hash_entries = []
for file_path in hash_files:
    basename = os.path.basename(file_path)
    hash_key = generate_sha256_hash(file_path)
    hash_entries.append((basename, hash_key))

# Creating dataframe
df = pd.DataFrame(hash_entries, columns=['file_name', 'hash_key'])

# Creating .hash_keys to store hash file data
outfile = os.path.join(output, '.hash_keys')
df.to_csv(outfile, index=False)
