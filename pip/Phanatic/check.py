import os
import sys
import pandas as pd
import hashlib

def generate_sha256_hash(path):
    with open(path, 'rb') as file:
        file_contents = file.read()
    hash_object = hashlib.sha256(file_contents)
    hex = hash_object.hexdigest()
    return hex

def check_task(output):

    # Initialising
    error = False

    # Initialising 
    hash_files = []

    # Find required files
    log = os.path.join(output, 'phanatic_log.tsv')
    raw = os.path.join(output, 'raw_data.csv')
    summary = os.path.join(output, 'combined_summary.csv')

    # Adding config and mapping files if they exist
    config = os.path.join(output, 'config.ini')
    mapping = os.path.join(output, 'host_mapping.csv')

    # Checking required files exist
    files = [log, raw, summary]
    for file in files:
        if os.path.exists(file):
            hash_files.append(file)
        else:
            print(f"ERROR: {os.path.basename(file)} does not exist, cannot perform checks")
            sys.exit(0)

    # Checking optional files exist
    files = [config, mapping]
    for file in files:
        if os.path.exists(file):
            hash_files.append(file)
        else:
            print(f"File: {os.path.basename(file)} does not exist")

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

    # Creating dataframes
    keys = os.path.join(output, '.hash_keys')
    original = pd.read_csv(keys)
    df = pd.DataFrame(hash_entries, columns=['file_name', 'hash_key_check'])

    # Checking size
    if len(original) == len(df):
        print("Same number of files hashed")
    elif len(original) > len(df):
        print("Files are missing from phage_genomes directory")
    elif len(original) < len(df):
        print("Files have been added to phage_genomes directory")

    # Checking files
    missing_files = []
    changed_files = []
    for index, row in original.iterrows():
        name = row['file_name']
        key = row['hash_key']

        # Missing
        if name not in df['file_name'].values:
            missing_files.append(name)
            continue
        
        # Check if the hash key is different
        check_key = df[df['file_name'] == name]['hash_key_check'].values[0]
        if key != check_key:
            changed_files.append(name)

    # Report
    if len(missing_files) > 0:
        for file in missing_files:
            print(f"Warning, file is missing: {file}")
            error = True
    
    if len(changed_files) > 0:
        for file in changed_files:
            print(f"Warning, file hash has changed: {file}")
            error = True

    if not error:
        print("No errors detected")

    return missing_files, changed_files
