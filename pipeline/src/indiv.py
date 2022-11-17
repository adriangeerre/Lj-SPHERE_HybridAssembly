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
#from src import validation
#from src import annotation

# External
#import glob, os, sys
import subprocess

# Functions
def init(read1, read2, long, prefix, threads):
    print(read1, read2, long, prefix)

    # Conda path
    conda_path = conda_path()[0]

    # QC
    qc.qc_illumina(illumina_1=read1, illumina_2=read2, out_dir=f"01-QC/{prefix}/Illumina", threads=threads, prefix=prefix, conda_path=conda_path)
    qc.qc_nanopore(nanopore=long, out_dir=f"01-QC/{prefix}", threads=threads, conda_path=conda_path)
    
    # Correction
    read_corr_1 = ".".join(read1.split("/")[-1].split(".")[:-2]) + ".fastq.00.0_0.cor.fastq.gz"
    read_corr_2 = ".".join(read2.split("/")[-1].split(".")[:-2]) + ".fastq.00.0_0.cor.fastq.gz"

    correction.correct_illumina(illumina_1=read1, illumina_2=read2, illumina_corr_1=read_corr_1, illumina_corr_2=read_corr_2, out_dir=f"10-Correction/{prefix}", threads=threads, conda_path=conda_path)
    correction.correct_nanopore(nanopore=long, illumina_corr=f"10-Correction/{prefix}/Illumina/corrected/illumina.corrected.fastq.gz", out_dir=f"10-Correction/{prefix}/Nanopore", threads=threads, conda_path=conda_path)

    # Assembly
    assembly.flye_assembly(nanopore_corr=f"10-Correction/{prefix}/Nanopore/nanopore.corrected.fasta", out_dir=f"20-Assembly/{prefix}", threads=threads, conda_path=conda_path)
    assembly.unicycler(assembly=f"20-Assembly/{prefix}/assembly.fasta", illumina_corr_1=read_corr_1, illumina_corr_2=read_corr_2, nanopore_corr=f"10-Correction/{prefix}/Nanopore/nanopore.corrected.fasta", out_dir=f"30-HybridAssembly/{prefix}", threads=threads, conda_path=conda_path)

    # AFTER CHECK!

    # Coverage
    coverage.align_illumina(assembly="30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), illumina_corr_1=read_corr_1, illumina_corr_2=read_corr_2, out_dir="30-HybridAssembly/{}/Align".format(prefix), threads=threads, conda_path=conda_path)
    coverage.align_nanopore(assembly, nanopore_corr, out_dir="30-HybridAssembly/{}/Align".format(prefix), threads=threads, conda_path=conda_path)
    coverage.coverage(in_dir, out_dir, conda_path=conda_path)
    coverage.plot_coverage(in_dir, out_dir, prefix=prefix, conda_path=conda_path)

def conda_path():
    # Exec terminal command
    cmd = "conda info | grep 'base environment' | awk '{print $4}'"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()

    # Out to string
    out = str(out.strip()).replace("'","")[1:]

    return [out, err]

def check_contig_number(hybrid_assembly, nanopore_draft):
    pass
