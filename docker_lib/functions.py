
import os
import sys
import datetime
import subprocess
import csv

# Logging function
def logfile(function, text, logfile):
    newline = '\n'
    container = "PhageOrder v0.0.2"
    date_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    report = f"{newline}[{container}]\t[{date_time}]\t[{function}]\t[{text}]"
    with open(logfile, 'a') as file:
        file.write(report)
