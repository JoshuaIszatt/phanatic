#!/usr/bin/env python

# Coordinator script by J.J.Iszatt

import os
import sys
import datetime
import subprocess
import configparser
import pandas as pd
from Bio import SeqIO
import csv
import functions as ji

# Reading input directory
input = '/assemble/input'
output = '/assemble/output'
index = '/assemble/database/PHROG_index.csv'
logs = '/assemble/output/phanatic_log.tsv'

# Creating logfile
ji.logfile("", "Checking input files", logs)
