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
from src import qc
from src import correction
from src import assembly
from src import coverage
from src import validation
#from src import annotation

# External
import glob, os, sys
from subprocess import check_call

# Functions
def init(samples, threads, run_coverage):
    pass
