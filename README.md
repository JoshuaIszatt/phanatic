# Phanatic v2.2.0
[![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)](https://pypi.org/project/PhageOrder/)
[![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white)](https://hub.docker.com/repository/docker/iszatt/phageorder/general)

This python script (phanatic.py) will run a docker container to assemble genomes 'de novo' using SPAdes (version number below in third party software).

## Associated publications:
1. [Genome Sequences of Two Lytic Staphylococcus aureus Bacteriophages Isolated from Wastewater](https://journals.asm.org/doi/10.1128/mra.00954-22)
2. [Genome Sequence of a Lytic Staphylococcus aureus Bacteriophage Isolated from Breast Milk](https://journals.asm.org/doi/10.1128/mra.00953-22)
3. [Complete Genome Sequences of Four Pseudomonas aeruginosa Bacteriophages: Kara-mokiny 8, Kara-mokiny 13, Kara-mokiny 16, and Boorn-mokiny 1](https://journals.asm.org/doi/10.1128/mra.00960-22)

## FUNCTIONS
The main functions of phanatic are:
* De novo assembly for phages
* Quality checks run using qualimap
* Extraction of 'Complete' and 'High-quality' contigs determined by CheckV
* CSV file with details of the assembly for each 

## Citation:
If you use this software please cite the source code as below and look at the third party software to cite the correct versions of software utilised by this container.
```
Iszatt J.(2023).Phanatic(v2.2.0)[Source code].Github:https://github.com/JoshuaIszatt/Phanatic
```

## Installation
Install using pip
```sh

```

Run the help command to see options 
```sh

```

## Output
* The assembled genome

## Third-party software
| Software | Version | Description | Please cite |
| -------- | -------- | -------- | -------- |
| SPAdes | 3.15.4 | -------- | -------- |
| bbmap | 38.18 | -------- | -------- |
| samtools | 1.6 | -------- | -------- |
| biopython | 1.78 | -------- | -------- |
| checkv | 1.0.1 | -------- | -------- |
| checkv-db | 1.5 | -------- | -------- |
| fastqc | 0.11.9 | Quality control for reads | -------- |

## Docker tags
https://hub.docker.com/r/iszatt

## License
[GNU AGPLv3](https://github.com/JoshuaIszatt/phanatic/blob/master/LICENSE.md)

## Config file
This is the default config file, copy this and specify its location to use your own with adjustments.
```
[phanatic]
image = iszatt/phanatic:2.2.0
author = 'Joshua J Iszatt'
citation = 'pending'

[pipeline]
normalise = True
filter = True
fastqc = True
barcode = True
clean_up = False

[system]
RAM = 24000m

[input]
r1_ext = _R1.fastq.gz
r2_ext = _R2.fastq.gz

[trim]
read_length = 150
trim_length = 12
minimum_length = 100

[merge]
minimum_insert = 120
minimum_overlap = 20

[normalise]
target_coverage = 250

[SPAdes]
memory_gb = 24
threads = 24

[filter]
filter_length = 4000

[barcoding]
prefix = phanatic
barcode_length = 5
```
