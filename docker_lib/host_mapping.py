#!/usr/bin/env python

# Host mapping script by J.J.Iszatt

import os
import sys
import pandas as pd
import configparser
import datetime

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

image = config["phanatic"]["image"]

# Phanatic settings
try:
    enable_normalise = config.getboolean("pipeline", "normalise")
    enable_filter = config.getboolean("pipeline", "filter")
    enable_qc = config.getboolean("pipeline", "fastqc")
    enable_barcodes = config.getboolean("pipeline", "barcode")
    enable_clean = config.getboolean("pipeline", "clean_up")
except ValueError:
    sys.exit("Config file incorrectly set, pipeline values must be booleans")

################################################################################
#_______________________________________________________________________________


# Functions
def logfile(function, text, logfile):
    newline = "\n"
    date_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    report = f"{newline}[{image}]\t[{date_time}]\t[{function}]\t[{text}]"
    with open(logfile, "a") as file:
        file.write(report)



################################################################################
#_______________________________________________________________________________


# Read mapping
logfile("Phanatic mapping process", "collating bacterial hosts", logs)
mapping_file = "/assemble/output/host_mapping.csv"
if os.path.isfile(mapping_file):
    df = pd.read_csv(mapping_file)
else:
    logfile("Phanatic mapping process", "ERROR, mapping file not found", logs)
    sys.exit("ERROR, mapping file not found")

# Building class to hold data
class Mapper(object):
    def __init__(self, host, phage, qc1, qc2, n1, n2):
        self.host = host;
        self.phage = phage;
        self.qc1 = qc1;
        self.qc2 = qc2;
        self.n1 = n1;
        self.n2 = n2;
    

# Checking for bacterial genomes

bacterial_genomes = list(df['host'])
print("Its worked so far...")
