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

# Arguments
#----------
class Arguments(object):

    def __init__(self):
        parser = argparse.ArgumentParser(
            description='Prokaryotic Hybrid Assembly',
            usage='''PHA.py <mode> [<args>]

    The pipeline can run in two modes:
        indiv     Process one individual strains
        multi     Process multiple strains
        ''')
        parser.add_argument('mode', help='Run mode')
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.mode):
            print('Unrecognized mode')
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.mode)()


    def indiv(self):
        # Parser
        parser = argparse.ArgumentParser(description='Process one individual strains')

        # Sub-Arguments
        parser.add_argument('-r1', '--illumina-forward', dest='short_forward', action='store', help='R1/Forward illumina reads', required=True)
        parser.add_argument('-r2', '--illumina-reverse', dest='short_reverse', action='store', help='R2/Reverse illumina reads', required=True)
        parser.add_argument('-n', '--nanopore', dest='long', action='store', help='Nanopore reads', required=True)
        parser.add_argument('-p', '--prefix', dest='prefix', action='store', help='Prefix', required=True)
        parser.add_argument('-t', '--threads', dest='threads', action='store', help='Threads (default: %(default)s)', default=1, type=int)
        parser.add_argument('-f', '--force', dest='boolean', action='store', help='Force recomputation', choices=['True','False'], default=False)

        # Settings
        args = parser.parse_args(sys.argv[2:])
        if args.prefix == None:
            parser.print_help()
        else:
            print('Processing strain: {}'.format(args.prefix))

        return args

    def multi(self):
        # Parser
        parser = argparse.ArgumentParser(description='Process multiple strains')

        # Sub-Arguments
        parser.add_argument('-l', '--list-samples', dest='samples', action='store', help='List of samples (process multiples samples)')
        parser.add_argument('-t', '--threads', dest='threads', action='store', help='Threads (default: %(default)s)', default=1, type=int)
        parser.add_argument('-f', '--force', dest='boolean', action='store', help='Force recomputation', choices=['True','False'], default=False)

        # Settings
        args = parser.parse_args(sys.argv[2:])
        if args.samples == None:
            parser.print_help()
        else:
            print('Processing strains in file: {}'.format(args.samples))

        return args


if __name__ == '__main__':
    args = Arguments()

    start = time.time()
    start_time = time.strftime("%Y-%m-%d %H:%M:%S")
