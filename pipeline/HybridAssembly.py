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

# Internal
from src import indiv, multi
#from src import qc
#from src import correction
#from src import assembly
#from src import coverage
#from src import validation
#from src import annotation

# External
#import os
import sys
#import json
import time
import argparse
#import subprocess
#import statistics as st # Python v3.4 or above

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

# Options
#--------
def params():
    # Parser
    parser = argparse.ArgumentParser(prog='HybridAssembly', description='Python pipeline to perform hybrid assembly of nanopore and illumina sequencing data. It is based on Benjamin\'s Perry pipeline and uses flye and Unicycler as the nanopore and hybrid assemblers, respectively.')

    # Subparsers
    subparsers = parser.add_subparsers(help="Sub-commands help")
    parser_indiv = subparsers.add_parser('indiv', help='Individual')
    parser_multi = subparsers.add_parser('multi', help='Multiple')
    
    # Arguments Indiv
    parser_indiv.add_argument('-r1', '--illumina-forward', dest='short_forward', action='store', help='R1/Forward illumina reads', required=True)
    parser_indiv.add_argument('-r2', '--illumina-reverse', dest='short_reverse', action='store', help='R2/Reverse illumina reads', required=True)
    parser_indiv.add_argument('-n', '--nanopore', dest='long', action='store', help='Nanopore reads', required=True)
    parser_indiv.add_argument('-p', '--prefix', dest='prefix', action='store', help='Prefix', required=True)
    parser_indiv.add_argument('-t', '--threads', dest='threads', action='store', help='Threads (default: %(default)s)', default=1, type=int)
    #parser_indiv.add_argument('-f', '--force', dest='boolean', action='store', help='Force recomputation', choices=['True','False'], default=False)

    parser_indiv.set_defaults(func="indiv")

    # Arguments Multi
    parser_multi.add_argument('-l', '--list-samples', dest='samples', action='store', help='List of samples (process multiples samples)', required=True)
    parser_multi.add_argument('-t', '--threads', dest='threads', action='store', help='Threads (default: %(default)s)', default=1, type=int)
    #parser_multi.add_argument('-f', '--force', dest='boolean', action='store', help='Force recomputation', choices=['True','False'], default=False)

    parser_multi.set_defaults(func="multi")

    # Args
    args = parser.parse_args()

    return args

# Checks
#-------

# Empty exec
def empty_exec(args):
    if len(vars(args)) == 0:
        print("Please, select a mode. For more information run \"python HybridAssembly.py --help\"")
        sys.exit()

# Execution
#----------
if __name__ == '__main__':
    
    start = time.time()
    start_time = time.strftime("%Y-%m-%d %H:%M:%S")

    # Parameters
    args = params()

    # Check if args are empty
    empty_exec(args)

    # Call mode
    if args.func == "indiv":
        indiv.init(read1=args.short_forward, read2=short_reverse, long=args.long, prefix=args.prefix, threads=args.threads) 
    elif args.func == "multi":
        multi.init(samples=args.samples, threads=args.threads)

