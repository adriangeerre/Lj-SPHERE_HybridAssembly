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

import json
# External
import os
import subprocess

# Internal
from src import annotation, assembly, correction, coverage, qc, validation


# Functions
def conda_path():
    # Exec terminal command
    cmd = "conda info | grep 'base environment' | awk '{print $4}'"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out,err = proc.communicate()

    # Out to string
    out = str(out.strip()).replace("'","")[1:]

    return [out, err]

def check_contig_number(hybrid_assembly, nanopore_draft):
    # Count contigs
    cmdd = f"grep -c '^>' {nanopore_draft}"
    cmdh = f"grep -c '^>' {hybrid_assembly}"
    dp = subprocess.Popen(cmdd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    hp = subprocess.Popen(cmdh, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    outd,errd = dp.communicate()
    outh,errh = hp.communicate()

    # Return category
    if int(outd) == int(outh):
        return "Equal"
    elif int(outd) < int(outh):
        return "Below"
    else:
        return "Above" # Exit

def init(read1, read2, long, prefix, genus, threads, run_coverage):

    # Conda path
    cpath = conda_path()[0]

    # QC
    #---

    # Illumina QC
    qc.qc_illumina(illumina_1=read1, illumina_2=read2, out_dir=f"01-QC/{prefix}/Illumina", threads=threads, conda_path=cpath)

    # Nanopore QC
    qc.qc_nanopore(nanopore=long, out_dir=f"01-QC/{prefix}", threads=threads, conda_path=cpath)
    
    # Correction
    #-----------

    # Define Illumina correct output
    read_corr_1 = ".".join(read1.split("/")[-1].split(".")[:-2]) + ".fastq.00.0_0.cor.fastq.gz"
    read_corr_2 = ".".join(read2.split("/")[-1].split(".")[:-2]) + ".fastq.00.0_0.cor.fastq.gz"

    # Correct Illumina reads
    correction.correct_illumina(illumina_1=read1, illumina_2=read2, illumina_corr_1=read_corr_1, illumina_corr_2=read_corr_2, out_dir=f"10-Correction/{prefix}", threads=threads, conda_path=cpath)

    # Correct nanopore reads
    correction.correct_nanopore(nanopore=long, illumina_corr=f"10-Correction/{prefix}/Illumina/corrected/illumina.corrected.fastq.gz", out_dir=f"10-Correction/{prefix}/Nanopore", threads=threads, conda_path=cpath)

    # Assembly
    #---------

    # Nanopore draft (Flye)
    assembly.flye_assembly(nanopore_corr=f"10-Correction/{prefix}/Nanopore/nanopore.corrected.fasta", out_dir=f"20-Assembly/{prefix}", threads=threads, conda_path=cpath)

    # Hybrid Assembly (Unicycler)
    assembly.unicycler(assembly=f"20-Assembly/{prefix}/assembly.fasta", illumina_corr_1=read_corr_1, illumina_corr_2=read_corr_2, nanopore_corr=f"10-Correction/{prefix}/Nanopore/nanopore.corrected.fasta", out_dir=f"30-HybridAssembly/{prefix}/unicycler", threads=threads, conda_path=cpath)

    # Check contigs
    #--------------
    d = f"20-Assembly/{prefix}/flye/assembly.fasta"
    h = f"30-HybridAssembly/{prefix}/unicycler/assembly.fasta"
    
    if os.path.exists(h) and os.path.exists(d):
        category = check_contig_number(h, d)
    elif not os.path.exists(h) and os.path.exists(d):
        category = "Below"

    # Define path
    #------------

    if category == "Equal":
        assembly_folder = "30-HybridAssembly"
        software = "unicycler"
    elif category == "Below":
        assembly_folder = "20-Assembly"
        software = "flye"

    # Coverage
    #---------
    
    if run_coverage == True: # Define title for Draft or HA
        # Align illumina
        coverage.align_illumina(assembly=f"{assembly_folder}/{prefix}/{software}/assembly.fasta", illumina_corr_1=read_corr_1, illumina_corr_2=read_corr_2, out_dir=f"{assembly_folder}/{prefix}/Align", threads=threads, conda_path=cpath)

        # Align nanopore
        coverage.align_nanopore(assembly=f"{assembly_folder}/{prefix}/{software}/assembly.fasta", nanopore_corr=f"10-Correction/{prefix}/Nanopore/nanopore.corrected.fasta", out_dir=f"{assembly_folder}/{prefix}/Align", threads=threads, conda_path=cpath)
        
        # Obtain coverage
        coverage.coverage(in_dir="{assembly_folder}/{prefix}/Align", out_dir=f"{assembly_folder}/{prefix}/Coverage", conda_path=cpath)

        # Plot coverage
        coverage.plot_coverage(in_dir=f"{assembly_folder}/{prefix}/Align", out_dir=f"{assembly_folder}/{prefix}/Coverage/CovPlots", prefix=prefix, conda_path=cpath)

    # Validation
    #-----------
    # Quast
    validation.quast(assembly=f"{assembly_folder}/{prefix}/{software}/assembly.fasta", out_dir=f"{assembly_folder}/{prefix}/Quast", threads=threads, conda_path=cpath)

    # Assembly validation
    validation.assembly_validation(assembly=f"{assembly_folder}/{prefix}/{software}/assembly.fasta", database=genus, out_dir=f"40-Validation/{prefix}/{software}", threads=threads, conda_path=cpath)

    # Check validation
    #-----------------

    # Validate BUSCO
    if os.path.isdir(f"40-Validation/{prefix}/{software}/Busco"):
        try:
            b = [i for i in os.listdir(f"40-Validation/{prefix}/{software}/Busco") if i[-5:] == ".json"][0]
            bf = open(f"40-Validation/{prefix}/{software}Busco/{b}")
            bd = json.load(bf)
            bv = validation.validate_busco(bd)
            bf.close()
        except:
            pass
            #print(f"{bcolors.FAIL}Error: 40-Validation/{folder}/unicycler/Busco summary file is missing or empty.{bcolors.ENDC}\n")

    # Validate CheckM
    if os.path.isdir(f"40-Validation/{prefix}/{software}/CheckM"):
        c = f"40-Validation/{prefix}/{software}/CheckM/results.tsv"
        try:
            cv = validation.validate_checkm(c)
        except:
            pass
            #print(f"{bcolors.FAIL}Error: {c} is missing or empty.{bcolors.ENDC}\n")


    # Annotation
    #-----------
    if "bv" in globals() and "cv" in globals():
        if bv == "Pass" and cv == "Pass":
            if os.path.isdir(f"{assembly_folder}/{prefix}/{software}"):
                out_dir_yaml = f"{assembly_folder}/{prefix}/{software}/{prefix}"
                if not os.path.exists(out_dir_yaml + ".submol.yml") and not os.path.exists(out_dir_yaml + ".input.yml"):
                    annotation.pgap_files_creator(genus = genus, assembly = f"{assembly_folder}/{prefix}/{software}/assembly.fasta", out_dir = out_dir_yaml)
                    
                if os.path.exists(out_dir_yaml + '.submol.yml') and os.path.exists(out_dir_yaml + '.input.yml') and genus != "NA":
                    annotation.annotation(input_yaml=f"{out_dir_yaml}.input.yml", out_dir=f"50-Annotation/{prefix}/{software}", memory=4, threads=1,  conda_path=cpath)

    if os.path.isdir(f"50-Annotation/{prefix}/{software}/annotation"):
        a = f"50-Annotation/{prefix}/{software}/annotation/annot.gff"
        try:
            av = validation.validate_pgap(a)
            if av == "Pass" and category == "Equal":
                # Create folder
                if not os.path.isdir("60-Genomes/Complete"): os.makedirs("60-Genomes/Complete")
                
                # Move assembly to Complete Genome
                os.system(f"cp {assembly_folder}/{prefix}/unicycler/assembly.fasta 60-Genomes/Complete/{prefix}.assembly")
                
                # Remove assembly from Improved Genome
                if os.path.exists(f"60-Genomes/Improved/{prefix}.assembly"): os.remove(f"60-Genomes/Improved/{prefix}.assembly.fasta")
                
            elif av == "Pass" and category == "Below":
                # Create folder
                if not os.path.isdir("60-Genomes/Improved"): os.makedirs("60-Genomes/Improved")
                
                # Move assembly to Improved Genome
                os.system(f"cp {assembly_folder}/{prefix}/flye/assembly.fasta 60-Genomes/Improved/{prefix}.assembly.fasta")
                
        except:
            if os.path.isdir(f"50-Annotation/{prefix}/annotation"):
                pass
				# print(f"{bcolors.FAIL}Error: {a} is missing or empty.{bcolors.ENDC}")

