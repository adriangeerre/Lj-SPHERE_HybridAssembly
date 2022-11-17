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
import json
import subprocess

# Internal
from src import annotation, assembly, correction, coverage, qc, validation


# Functions
def init(read1, read2, long, prefix, genus, threads):
    print(read1, read2, long, prefix)

    # Conda path
    conda_path = conda_path()[0]

    # QC
    #---

    # Illumina QC
    qc.qc_illumina(illumina_1=read1, illumina_2=read2, out_dir=f"01-QC/{prefix}/Illumina", threads=threads, prefix=prefix, conda_path=conda_path)

    # Nanopore QC
    qc.qc_nanopore(nanopore=long, out_dir=f"01-QC/{prefix}", threads=threads, conda_path=conda_path)
    
    # Correction
    #-----------

    # Define Illumina correct output
    read_corr_1 = ".".join(read1.split("/")[-1].split(".")[:-2]) + ".fastq.00.0_0.cor.fastq.gz"
    read_corr_2 = ".".join(read2.split("/")[-1].split(".")[:-2]) + ".fastq.00.0_0.cor.fastq.gz"

    # Correct Illumina reads
    correction.correct_illumina(illumina_1=read1, illumina_2=read2, illumina_corr_1=read_corr_1, illumina_corr_2=read_corr_2, out_dir=f"10-Correction/{prefix}", threads=threads, conda_path=conda_path)

    # Correct nanopore reads
    correction.correct_nanopore(nanopore=long, illumina_corr=f"10-Correction/{prefix}/Illumina/corrected/illumina.corrected.fastq.gz", out_dir=f"10-Correction/{prefix}/Nanopore", threads=threads, conda_path=conda_path)

    # Assembly
    #---------

    # Nanopore draft (Flye)
    assembly.flye_assembly(nanopore_corr=f"10-Correction/{prefix}/Nanopore/nanopore.corrected.fasta", out_dir=f"20-Assembly/{prefix}", threads=threads, conda_path=conda_path)

    # Hybrid Assembly (Unicycler)
    assembly.unicycler(assembly=f"20-Assembly/{prefix}/assembly.fasta", illumina_corr_1=read_corr_1, illumina_corr_2=read_corr_2, nanopore_corr=f"10-Correction/{prefix}/Nanopore/nanopore.corrected.fasta", out_dir=f"30-HybridAssembly/{prefix}/unicycler", threads=threads, conda_path=conda_path)

    # Check contigs
    #--------------
    d = f"20-Assembly/{prefix}/flye/assembly.fasta"
    h = f"30-HybridAssembly/{prefix}/unicycler/assembly.fasta"
    
    cmdd = f"grep -c '^>' {d}"
    cmdh = f"grep -c '^>' {h}"
    dp = subprocess.Popen(cmdd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    hp = subprocess.Popen(cmdh, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    outd,errd = dp.communicate()
    outh,errh = hp.communicate()

    if int(outd) == int(outh):
        pass

    # Coverage
    #---------

    # Align illumina
    coverage.align_illumina(assembly=f"30-HybridAssembly/{prefix}/unicycler/assembly.fasta", illumina_corr_1=read_corr_1, illumina_corr_2=read_corr_2, out_dir=f"30-HybridAssembly/{prefix}/Align", threads=threads, conda_path=conda_path)

    # Align nanopore
    coverage.align_nanopore(assembly=f"30-HybridAssembly/{prefix}/unicycler/assembly.fasta", nanopore_corr=f"10-Correction/{prefix}/Nanopore/nanopore.corrected.fasta", out_dir=f"30-HybridAssembly/{prefix}/Align", threads=threads, conda_path=conda_path)
    
    # Obtain coverage
    coverage.coverage(in_dir="30-HybridAssembly/{prefix}/Align", out_dir=f"30-HybridAssembly/{prefix}/Coverage", conda_path=conda_path)

    # Plot coverage
    coverage.plot_coverage(in_dir=f"30-HybridAssembly/{prefix}/Align", out_dir=f"30-HybridAssembly/{prefix}/Coverage/CovPlots", prefix=prefix, conda_path=conda_path)

    # Validation
    #-----------
    # Quast
    validation.quast(assembly=f"30-HybridAssembly/{prefix}/unicycler/assembly.fasta", out_dir=f"30-HybridAssembly/{prefix}/Quast", threads=threads, conda_path=conda_path)

    # Assembly validation
    validation.assembly_validation(assembly=f"30-HybridAssembly/{prefix}/unicycler/assembly.fasta", database=genus, out_dir=f"40-Validation/{prefix}/unicycler", threads=threads, conda_path=conda_path)

    # Check validation
    #-----------------

    # Validate BUSCO
    if os.path.isdir(f"40-Validation/{prefix}/unicycler/Busco"):
        try:
            b = [i for i in os.listdir(f"40-Validation/{prefix}/unicycler/Busco") if i[-5:] == ".json"][0]
            bf = open(f"40-Validation/{prefix}/unicycler/Busco/{b}")
            bd = json.load(bf)
            bv = validation.validate_busco(bd)
            bf.close()
        except:
            pass
            #print(f"{bcolors.FAIL}Error: 40-Validation/{folder}/unicycler/Busco summary file is missing or empty.{bcolors.ENDC}\n")

    # Validate CheckM
    if os.path.isdir(f"40-Validation/{prefix}/unicycler/CheckM"):
        c = f"40-Validation/{prefix}/unicycler/CheckM/results.tsv"
        try:
            cv = validation.validate_checkm(c)
        except:
            pass
            #print(f"{bcolors.FAIL}Error: {c} is missing or empty.{bcolors.ENDC}\n")


    # Annotation
    #-----------
    if "bv" in globals() and "cv" in globals():
        if bv == "Pass" and cv == "Pass":
            if os.path.isdir(f"30-HybridAssembly/{prefix}/unicycler"):
                out_dir_yaml = f"30-HybridAssembly/{prefix}/unicycler/{prefix}"
                if not os.path.exists(out_dir_yaml + ".submol.yml") and not os.path.exists(out_dir_yaml + ".input.yml"):
                    annotation.pgap_files_creator(genus = genus, assembly = f"30-HybridAssembly/{prefix}/unicycler/assembly.fasta", out_dir = out_dir_yaml)
                    
                if os.path.exists(out_dir_yaml + '.submol.yml') and os.path.exists(out_dir_yaml + '.input.yml') and genus != "NA":
                    annotation.annotation(input_yaml=f"{out_dir_yaml}.input.yml", out_dir=f"50-Annotation/{prefix}/unicycler", memory=4, threads=1,  conda_path=conda_path)

def conda_path():
    # Exec terminal command
    cmd = "conda info | grep 'base environment' | awk '{print $4}'"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out,err = proc.communicate()

    # Out to string
    out = str(out.strip()).replace("'","")[1:]

    return [out, err]

def check_contig_number(hybrid_assembly, nanopore_draft):
    pass
