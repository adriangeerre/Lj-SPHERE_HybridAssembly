# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  	: Adrián Gómez Repollés
# Email			: adrian.gomez@mbg.au.dk
# Created Date	: 15/11/2022
# version 		: '1.0'
# ---------------------------------------------------------------------------
""" Python pipeline to perform hybrid assembly of nanopore and illumina sequen
-cing data. It is based on Benjamin Perry pipeline and uses flye and Unicy
-cler as the nanopore and hybrid assemblers, respectively.""" 
# ---------------------------------------------------------------------------
# Imports 
# ---------------------------------------------------------------------------

# External
import os
import sys

# Internal
from src import indiv

# Color Scheme
#-------------
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# Functions
def read_file(sample_file):
    # Check if file exists
    if os.path.exists(sample_file) and os.stat(sample_file).st_size != 0:
        # Create list of lists
        samples = []

        # Read file
        f = open(sample_file, "r")
        for row in f:
            if row[0] != "#":
                row = row.strip().split("\t")
                if len(row) >= 6:
                    samples.append(row)
        f.close()

    else:
        print(f"{bcolors.FAIL}Error: {sample_file} is missing or empty.{bcolors.ENDC}\n")
        sys.exit()

    return samples

def init(sample_file, threads, memory, run_coverage):
    # Read and check file
    samples = read_file(sample_file)

    # Loop samples
    for sample in samples:
        try:
            indiv.init(read1=sample[0], read2=sample[1], long=sample[2], prefix=sample[3], genus=sample[4], threads=threads, memory=memory, run_coverage=run_coverage)
        except:
            print(f"{bcolors.WARNING}Warning: {sample[0]} computation failed.{bcolors.ENDC}\n")
            continue

