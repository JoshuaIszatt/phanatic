#!/usr/bin/env python

# Coordinator script by J.J.Iszatt

import os
import sys
import datetime
import subprocess
import pandas as pd
from Bio import SeqIO
import csv

# Reading input directory
input = '/lab/input'
output = '/lab/output'
index = '/lab/database/PHROG_index.csv'
logs = '/lab/output/docker_log.tsv'