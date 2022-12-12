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
import logging

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

def init(read1, read2, long, prefix, genus, threads, memory, run_coverage):

    # Conda path
    cpath = conda_path()[0]

    # Logging
    if os.path.isdir("01-Logs") == False: os.makedirs("01-Logs")
    logging.basicConfig(filename=f'01-Logs/{prefix}.log', format='%(levelname)s:%(message)s', encoding='utf-8', level=logging.INFO)

    # QC
    #---

    logging.info('QC:')

    # Illumina QC
    r1qc=read1.split("/")[-1].replace(".fastq.gz","_fastqc.zip")
    r2qc=read2.split("/")[-1].replace(".fastq.gz","_fastqc.zip")
    if not os.path.exists(f"01-QC/{prefix}/Illumina/{r1qc}") or not os.path.exists(f"01-QC/{prefix}/Illumina/{r2qc}"):
        qc.qc_illumina(illumina_1=read1, illumina_2=read2, out_dir=f"01-QC/{prefix}/Illumina", threads=threads, conda_path=cpath, logfile=f"01-Logs/{prefix}.log")

    # Nanopore QC
    if not os.path.exists(f"01-QC/{prefix}/Nanopore/{prefix}_NanoStats.txt"):
        qc.qc_nanopore(nanopore=long, out_dir=f"01-QC/{prefix}/Nanopore", threads=threads, conda_path=cpath, logfile=f"01-Logs/{prefix}.log")
    
    # Correction
    #-----------

    #logging.info('### Illumina Correction ###')

    # Define Illumina correct output
    read_corr_1 = ".".join(read1.split("/")[-1].split(".")[:-2]) + ".fastq.00.0_0.cor.fastq.gz"
    read_corr_2 = ".".join(read2.split("/")[-1].split(".")[:-2]) + ".fastq.00.0_0.cor.fastq.gz"

    # Correct Illumina reads
    if not os.path.exists(f"10-Correction/{prefix}/Illumina/corrected/{read_corr_1}") or not os.path.exists(f"10-Correction/{prefix}/Illumina/corrected/{read_corr_2}"):
        correction.correct_illumina(illumina_1=read1, illumina_2=read2, illumina_corr_1=read_corr_1, illumina_corr_2=read_corr_2, out_dir=f"10-Correction/{prefix}/Illumina", threads=threads, memory=memory, conda_path=cpath, logfile=f"01-Logs/{prefix}.log")

    logging.info('Nanopore Correction:')

    # Correct nanopore reads
    if not os.path.exists(f"10-Correction/{prefix}/Nanopore/nanopore.corrected.fasta"):
        correction.correct_nanopore(nanopore=long, illumina_corr=f"10-Correction/{prefix}/Illumina/corrected/illumina.corrected.fastq.gz", out_dir=f"10-Correction/{prefix}/Nanopore", threads=threads, conda_path=cpath, logfile=f"01-Logs/{prefix}.log")

    # Assembly
    #---------

    logging.info('Nanopore Draft Assembly:')

    # Nanopore draft (Flye)
    if not os.path.exists(f"20-Assembly/{prefix}/flye/assembly.fasta"):
        assembly.flye_assembly(nanopore_corr=f"10-Correction/{prefix}/Nanopore/nanopore.corrected.fasta", out_dir=f"20-Assembly/{prefix}/flye", threads=threads, conda_path=cpath, logfile=f"01-Logs/{prefix}.log")

    logging.info('Hybrid Assembly:')

    # Hybrid Assembly (Unicycler)
    if not os.path.exists(f"30-HybridAssembly/{prefix}/unicycler/assembly.fasta"):
        assembly.unicycler(assembly=f"20-Assembly/{prefix}/flye/assembly.fasta", illumina_corr_1=f"10-Correction/{prefix}/Illumina/corrected/{read_corr_1}", illumina_corr_2=f"10-Correction/{prefix}/Illumina/corrected/{read_corr_2}", nanopore_corr=f"10-Correction/{prefix}/Nanopore/nanopore.corrected.fasta", out_dir=f"30-HybridAssembly/{prefix}/unicycler", threads=threads, conda_path=cpath, logfile=f"01-Logs/{prefix}.log")

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
        coverage.align_illumina(assembly=f"{assembly_folder}/{prefix}/{software}/assembly.fasta", illumina_corr_1=read_corr_1, illumina_corr_2=read_corr_2, out_dir=f"{assembly_folder}/{prefix}/Align", threads=threads, conda_path=cpath, logfile=f"01-Logs/{prefix}.log")

        # Align nanopore
        coverage.align_nanopore(assembly=f"{assembly_folder}/{prefix}/{software}/assembly.fasta", nanopore_corr=f"10-Correction/{prefix}/Nanopore/nanopore.corrected.fasta", out_dir=f"{assembly_folder}/{prefix}/Align", threads=threads, conda_path=cpath, logfile=f"01-Logs/{prefix}.log")
        
        # Obtain coverage
        coverage.coverage(in_dir="{assembly_folder}/{prefix}/Align", out_dir=f"{assembly_folder}/{prefix}/Coverage", conda_path=cpath, logfile=f"01-Logs/{prefix}.log")

        # Plot coverage
        coverage.plot_coverage(in_dir=f"{assembly_folder}/{prefix}/Align", out_dir=f"{assembly_folder}/{prefix}/Coverage/CovPlots", prefix=prefix, conda_path=cpath, logfile=f"01-Logs/{prefix}.log")

    # Validation
    #-----------
    # Quast
    if not os.path.exists(f"{assembly_folder}/{prefix}/Quast/report.txt"):
        validation.quast(assembly=f"{assembly_folder}/{prefix}/{software}/assembly.fasta", out_dir=f"{assembly_folder}/{prefix}/Quast", threads=threads, conda_path=cpath, logfile=f"01-Logs/{prefix}.log")

    # Assembly validation
    validation.assembly_validation(assembly=f"{assembly_folder}/{prefix}/{software}/assembly.fasta", database=genus, out_dir=f"40-Validation/{prefix}/{software}", threads=threads, conda_path=cpath, logfile=f"01-Logs/{prefix}.log")

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
                    annotation.annotation(input_yaml=f"{out_dir_yaml}.input.yml", out_dir=f"50-Annotation/{prefix}/{software}", memory=4, threads=1,  conda_path=cpath, logfile=f"01-Logs/{prefix}.log")

            if os.path.isdir(f"50-Annotation/{prefix}/{software}/annotation"):
                a = f"50-Annotation/{prefix}/{software}/annotation/annot.gff"
                try:
                    av = validation.validate_pgap(a)                
                except:
                    if os.path.isdir(f"50-Annotation/{prefix}/annotation"):
                        pass
                        # print(f"{bcolors.FAIL}Error: {a} is missing or empty.{bcolors.ENDC}")
        elif bv == "Fail" or cv == "Fail":
            pass

            # 16S genus identification
            #new_genus = taxonomy.taxonomy(assembly=f"{assembly_folder}/{prefix}/{software}/assembly.fasta")

            if genus != new_genus and os.path.isdir(f"{assembly_folder}/{prefix}/{software}"):
                out_dir_yaml = f"{assembly_folder}/{prefix}/{software}/{prefix}"
                if not os.path.exists(out_dir_yaml + ".submol.yml") and not os.path.exists(out_dir_yaml + ".input.yml"):
                    annotation.pgap_files_creator(genus = new_genus, assembly = f"{assembly_folder}/{prefix}/{software}/assembly.fasta", out_dir = out_dir_yaml)
                    
                if os.path.exists(out_dir_yaml + '.submol.yml') and os.path.exists(out_dir_yaml + '.input.yml') and genus != "NA":
                    annotation.annotation(input_yaml=f"{out_dir_yaml}.input.yml", out_dir=f"50-Annotation/{prefix}/{software}", memory=4, threads=1,  conda_path=cpath, logfile=f"01-Logs/{prefix}.log")

                if os.path.isdir(f"50-Annotation/{prefix}/{software}/annotation"):
                    a = f"50-Annotation/{prefix}/{software}/annotation/annot.gff"
                    try:
                        av = validation.validate_pgap(a)
                    except:
                        if os.path.isdir(f"50-Annotation/{prefix}/annotation"):
                            pass
                            # print(f"{bcolors.FAIL}Error: {a} is missing or empty.{bcolors.ENDC}")
        else:
            # Assembly did not pass
            pass

        # Organize genomes
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