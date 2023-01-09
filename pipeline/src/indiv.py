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

# Log Header
#--------
def exec_info(read1, read2, long, prefix, order, genus, threads, memory, run_coverage):
    # Define longest argument value
    maxlen = max([len(str(i)) for i in [read1, read2, long, prefix, order, genus, threads, memory, run_coverage]])

    maxlen += 14
    message = (
    f"\n{'-'*maxlen}\n"
    f"\nSummary\n"
    f"{' '*2}Short reads:\n{' '*4}{read1}\n{' '*4}{read2}\n"
    f"{' '*2}Long reads: {long}\n"
    f"{' '*2}Prefix: {prefix}\n"
    f"{' '*2}Order: {order}\n"
    f"{' '*2}Genus: {genus}\n"
    f"{' '*2}Threads: {threads}\n"
    f"{' '*2}Memory: {memory}\n"
    f"{' '*2}Run Coverage: {run_coverage}\n"
    f"\n{'-'*maxlen}\n\n"
    )
    return(message)

# Functions
def coverage(assembly_folder, prefix, software, reads_1, reads_2, threads, cpath):
    # Align illumina
    coverage.align_illumina(assembly=f"{assembly_folder}/{prefix}/{software}/assembly.fasta", illumina_corr_1=reads_1, illumina_corr_2=reads_2, out_dir=f"{assembly_folder}/{prefix}/Align", threads=threads, conda_path=cpath, logfile=f".logs/{prefix}.log")

    # Align nanopore
    coverage.align_nanopore(assembly=f"{assembly_folder}/{prefix}/{software}/assembly.fasta", nanopore_corr=f"10-Correction/{prefix}/Nanopore/nanopore.corrected.fasta", out_dir=f"{assembly_folder}/{prefix}/Align", threads=threads, conda_path=cpath, logfile=f".logs/{prefix}.log")
    
    # Obtain coverage
    coverage.coverage(in_dir="{assembly_folder}/{prefix}/Align", out_dir=f"{assembly_folder}/{prefix}/Coverage", conda_path=cpath, logfile=f".logs/{prefix}.log")

    # Plot coverage
    coverage.plot_coverage(in_dir=f"{assembly_folder}/{prefix}/Align", out_dir=f"{assembly_folder}/{prefix}/Coverage/CovPlots", prefix=prefix, conda_path=cpath, logfile=f".logs/{prefix}.log")

def conda_path():
    # Exec terminal command
    cmd = "conda info | grep 'base environment' | awk '{print $4}'"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out,err = proc.communicate()

    # Out to string
    out = str(out.strip()).replace("'","")[1:]

    return [out, err]

def init(read1, read2, long, prefix, order, genus, threads, memory, run_coverage):

    # Conda path
    cpath = conda_path()[0]

    # Logging
    # -------
    if os.path.isdir(".logs") == False: os.makedirs(".logs") # Create folder
    if os.path.exists(f".logs/{prefix}.log"): os.remove(f".logs/{prefix}.log") # Re-new log file
    LOG_FILE=f".logs/{prefix}.log"
    LOG_FORMAT = "%(asctime)s %(levelname)-8s %(message)s"
    LOG_LEVEL = "INFO"
    LOG_DATEFMT = "%d/%m/%Y %H:%M:%S"
    logging.basicConfig(filename=LOG_FILE, format=LOG_FORMAT, level=getattr(logging, LOG_LEVEL),  datefmt=LOG_DATEFMT)
    logger = logging.getLogger()

    # Logging header
    f = open(LOG_FILE, "a")
    f.write(exec_info(read1=read1, read2=read2, long=long, prefix=prefix, order=order, genus=genus, threads=threads, memory=memory, run_coverage=run_coverage))
    f.close()

    # QC
    #---
    logger.info('**** QC ****')

    # Illumina QC
    r1qc=read1.split("/")[-1].replace(".fastq.gz","_fastqc.zip")
    r2qc=read2.split("/")[-1].replace(".fastq.gz","_fastqc.zip")
    if not os.path.exists(f"01-QC/{prefix}/Illumina/{r1qc}") or not os.path.exists(f"01-QC/{prefix}/Illumina/{r2qc}"):
        qc.qc_illumina(illumina_1=read1, illumina_2=read2, out_dir=f"01-QC/{prefix}/Illumina", threads=threads, conda_path=cpath, logfile=f".logs/{prefix}.log")

    # Nanopore QC
    if not os.path.exists(f"01-QC/{prefix}/Nanopore/{prefix}_NanoStats.txt"):
        qc.qc_nanopore(nanopore=long, out_dir=f"01-QC/{prefix}/Nanopore", threads=threads, conda_path=cpath, logfile=f".logs/{prefix}.log")
    
    # Correction
    #-----------
    logger.info('**** Illumina Correction ****')

    # Define Illumina correct output
    read_corr_1 = ".".join(read1.split("/")[-1].split(".")[:-2]) + ".fastq.00.0_0.cor.fastq.gz"
    read_corr_2 = ".".join(read2.split("/")[-1].split(".")[:-2]) + ".fastq.00.0_0.cor.fastq.gz"

    # Correct Illumina reads
    if not os.path.exists(f"10-Correction/{prefix}/Illumina/corrected/{read_corr_1}") or not os.path.exists(f"10-Correction/{prefix}/Illumina/corrected/{read_corr_2}"):
        correction.correct_illumina(illumina_1=read1, illumina_2=read2, illumina_corr_1=read_corr_1, illumina_corr_2=read_corr_2, out_dir=f"10-Correction/{prefix}/Illumina", threads=threads, memory=memory, conda_path=cpath, logfile=f".logs/{prefix}.log")

    logger.info('**** Nanopore Correction ****')

    # Correct nanopore reads
    if not os.path.exists(f"10-Correction/{prefix}/Nanopore/nanopore.corrected.fasta"):
        correction.correct_nanopore(nanopore=long, illumina_corr=f"10-Correction/{prefix}/Illumina/corrected/illumina.corrected.fastq.gz", out_dir=f"10-Correction/{prefix}/Nanopore", threads=threads, conda_path=cpath, logfile=f".logs/{prefix}.log")

    # Assembly
    #---------
    logger.info('**** Nanopore Draft Assembly ****')

    # Nanopore draft (Flye)
    if not os.path.exists(f"20-Assembly/{prefix}/flye/assembly.fasta"):
        assembly.flye_assembly(nanopore_corr=f"10-Correction/{prefix}/Nanopore/nanopore.corrected.fasta", out_dir=f"20-Assembly/{prefix}/flye", threads=threads, conda_path=cpath, logfile=f".logs/{prefix}.log")

    logger.info('**** Hybrid Assembly ****')

    # Hybrid Assembly (Unicycler)
    if not os.path.exists(f"30-HybridAssembly/{prefix}/unicycler/assembly.fasta"):
        assembly.unicycler(assembly=f"20-Assembly/{prefix}/flye/assembly.fasta", illumina_corr_1=f"10-Correction/{prefix}/Illumina/corrected/{read_corr_1}", illumina_corr_2=f"10-Correction/{prefix}/Illumina/corrected/{read_corr_2}", nanopore_corr=f"10-Correction/{prefix}/Nanopore/nanopore.corrected.fasta", out_dir=f"30-HybridAssembly/{prefix}/unicycler", threads=threads, conda_path=cpath, logfile=f".logs/{prefix}.log")

    # Validation
    #-----------
    logger.info('**** Hybrid Assembly summary: Quast ****')

    # Folders
    assembly_folder="30-HybridAssembly"
    software="unicycler"

    # Quast
    if not os.path.exists(f"{assembly_folder}/{prefix}/Quast/report.txt"):
        validation.quast(assembly=f"{assembly_folder}/{prefix}/{software}/assembly.fasta", out_dir=f"{assembly_folder}/{prefix}/Quast", threads=threads, conda_path=cpath, logfile=f".logs/{prefix}.log")

    logger.info('**** Hybrid Assembly validation: Busco and CheckM ****')

    # Assembly validation
    busco_dict = {'Actinomycetales': 'actinobacteria_class_odb10', 'Flavobacteriales': 'flavobacteriales_odb10', 'Bacillales': 'bacillales_odb10', 'Burkholderiales': 'burkholderiales_odb10', 'Caulobacterales': 'alphaproteobacteria_odb10', 'Rhizobiales': 'rhizobiales_odb10', 'Sphingomonadales': 'sphingomonadales_odb10', 'Pseudomonadales': 'pseudomonadales_odb10', 'Xanthomonadales': 'xanthomonadales_odb10'}
    
    try:
        # Define Busco database
        database=busco_dict[order]
        if not os.path.exists(f"40-Validation/{prefix}/{software}/Busco/short_summary.specific.{database}.Busco.txt") or not os.path.exists(f"40-Validation/{prefix}/{software}/CheckM/results.tsv") :
            validation.assembly_validation(assembly=f"{assembly_folder}/{prefix}/{software}/assembly.fasta", database=database, out_dir=f"40-Validation/{prefix}/{software}", threads=threads, conda_path=cpath, logfile=f".logs/{prefix}.log")
    except:
        logger.error(f"Database for {order} not considered for Busco execution. You can modify the object \"busco_dict\" to include the required database. Please, check \"busco --list-datasets\" to know which Order should be selected.")
        logger.warning("If CheckM failed with \"FileNotFoundError: [Errno 2] No such file or directory: '~/.checkm/hmms/phylo.hmm'\", use 'checkm data setRoot <checkm_data_dir>' to specify the location of CheckM database files.")
        pass

    # Check validation
    #-----------------

    # Validate BUSCO
    if os.path.isdir(f"40-Validation/{prefix}/{software}/Busco"):
        try:
            bha = [i for i in os.listdir(f"40-Validation/{prefix}/{software}/Busco") if i[-5:] == ".json"][0]
            bfha = open(f"40-Validation/{prefix}/{software}Busco/{bha}")
            bdha = json.load(bfha)
            bvha = validation.validate_busco(bdha)
            bfha.close()
            logger.info(f'Busco validation: {bvha}')
        except:
            logger.error(f"Folder \"40-Validation/{prefix}/{software}/Busco\" is missing or empty.")
            pass

    # Validate CheckM
    if os.path.isdir(f"40-Validation/{prefix}/{software}/CheckM"):
        c = f"40-Validation/{prefix}/{software}/CheckM/results.tsv"
        try:
            cvha = validation.validate_checkm(c)
            logger.info(f'CheckM validation: {cvha}')
        except:
            logger.error(f"File \"{c}\" is missing or empty.")
            pass

    # Annotation
    #-----------
    if "bvha" in globals() and "cvha" in globals():
        if bvha == "Pass" and cvha == "Pass":
            logger.info('**** Hybrid Assembly annotation: PGAP ****')

            # Annotation of Hybrid Assembly
            if os.path.isdir(f"{assembly_folder}/{prefix}/{software}"):
                out_dir_yaml = f"{assembly_folder}/{prefix}/{software}/{prefix}"
                if not os.path.exists(out_dir_yaml + ".submol.yml") and not os.path.exists(out_dir_yaml + ".input.yml"):
                    annotation.pgap_files_creator(genus = genus, assembly = f"{assembly_folder}/{prefix}/{software}/assembly.fasta", out_dir = out_dir_yaml)
                    
                if os.path.exists(out_dir_yaml + '.submol.yml') and os.path.exists(out_dir_yaml + '.input.yml') and genus != "NA":
                    annotation.annotation(input_yaml=f"{out_dir_yaml}.input.yml", out_dir=f"50-Annotation/{prefix}/{software}", memory=4, threads=1,  conda_path=cpath, logfile=f".logs/{prefix}.log")

            # Validate Hybrid Assembly annotation
            if os.path.isdir(f"50-Annotation/{prefix}/{software}/annotation"):
                a = f"50-Annotation/{prefix}/{software}/annotation/annot.gff"
                try:
                    avha = validation.validate_pgap(a)                
                except:
                    if os.path.isdir(f"50-Annotation/{prefix}/annotation"):
                        logger.error(f"File \"{a}\" is missing or empty.")
                        pass

            # Complete Genome
            if "avha" in globals():
                if avha == "Pass":
                    # Create folder
                    if not os.path.isdir("60-Genomes/Complete"): os.makedirs("60-Genomes/Complete")
                    
                    # Move assembly to Complete Genome
                    os.system(f"cp {assembly_folder}/{prefix}/unicycler/assembly.fasta 60-Genomes/Complete/{prefix}.assembly")
                    
                    # Remove assembly from Improved Genome
                    if os.path.exists(f"60-Genomes/Improved/{prefix}.assembly"): os.remove(f"60-Genomes/Improved/{prefix}.assembly.fasta")

                    # Coverage
                    if run_coverage == True:
                        # Logger
                        logger.info('**** Running coverage HA ****')

                        # Call
                        coverage(assembly_folder=assembly_folder, prefix=prefix, software=software, reads_1=read_corr_1, reads2=read_corr_2, threads=threads, cpath=cpath)

                else:
                    logger.warning(f"Sample {prefix} failed the hybrid assembly annotation validation.")
            else:
                logger.error(f"Hybrid assembly annotation validation for sample {prefix} was not completed. Please, check the logs.")
        
        # HA assembly or annotation failed (all variables already checked in globals)
        if bvha == "Fail" or cvha == "Fail" or avha == "Fail":
            # Validation failed (Draft)
            logger.warning(f"Sample {prefix} failed the assembly or annotation validation at the hybrid assembly level. Instead, the nanopore draft will be validated in search to recover an assembly for {prefix}.")

            # Use Draft
            logger.info('**** Draft Assembly validation: Busco and CheckM ****')

            # Folders
            assembly_folder="20-Assembly"
            software="flye"

            try:
                if not os.path.exists(f"40-Validation/{prefix}/{software}/Busco/short_summary.specific.{database}.Busco.txt") or not os.path.exists(f"40-Validation/{prefix}/{software}/CheckM/results.tsv") :
                    validation.assembly_validation(assembly=f"{assembly_folder}/{prefix}/{software}/assembly.fasta", database=database, out_dir=f"40-Validation/{prefix}/{software}", threads=threads, conda_path=cpath, logfile=f".logs/{prefix}.log")
            except:
                logger.error(f"Database for {order} not considered for Busco execution. You can modify the object \"busco_dict\" to include the required database. Please, check \"busco --list-datasets\" to know which Order should be selected.")
                pass

            # Validate Draft BUSCO
            if os.path.isdir(f"40-Validation/{prefix}/{software}/Busco"):
                try:
                    bd = [i for i in os.listdir(f"40-Validation/{prefix}/{software}/Busco") if i[-5:] == ".json"][0]
                    bfd = open(f"40-Validation/{prefix}/{software}Busco/{bd}")
                    bdd = json.load(bfd)
                    bvd = validation.validate_busco(bdd)
                    bfd.close()
                    logger.info(f'Busco Draft validation: {bvd}')
                except:
                    logger.error(f"Folder \"40-Validation/{prefix}/{software}/Busco\" is missing or empty.")
                    pass

            # Validate Draft CheckM
            if os.path.isdir(f"40-Validation/{prefix}/{software}/CheckM"):
                c = f"40-Validation/{prefix}/{software}/CheckM/results.tsv"
                try:
                    cvd = validation.validate_checkm(c)
                    logger.info(f'CheckM Draft validation: {cvd}')
                except:
                    logger.error(f"File \"{c}\" is missing or empty.")
                    pass

            # Annotation
            if "bvd" in globals() and "cvd" in globals():
                if bvd == "Pass" and cvd == "Pass":
                    logger.info('**** Draft Assembly annotation: PGAP ****')

                    # Annotation of Hybrid Assembly
                    if os.path.isdir(f"{assembly_folder}/{prefix}/{software}"):
                        out_dir_yaml = f"{assembly_folder}/{prefix}/{software}/{prefix}"
                        if not os.path.exists(out_dir_yaml + ".submol.yml") and not os.path.exists(out_dir_yaml + ".input.yml"):
                            annotation.pgap_files_creator(genus = genus, assembly = f"{assembly_folder}/{prefix}/{software}/assembly.fasta", out_dir = out_dir_yaml)
                            
                        if os.path.exists(out_dir_yaml + '.submol.yml') and os.path.exists(out_dir_yaml + '.input.yml') and genus != "NA":
                            annotation.annotation(input_yaml=f"{out_dir_yaml}.input.yml", out_dir=f"50-Annotation/{prefix}/{software}", memory=4, threads=1,  conda_path=cpath, logfile=f".logs/{prefix}.log")

                    # Validate Hybrid Assembly annotation
                    if os.path.isdir(f"50-Annotation/{prefix}/{software}/annotation"):
                        a = f"50-Annotation/{prefix}/{software}/annotation/annot.gff"
                        try:
                            avd = validation.validate_pgap(a)                
                        except:
                            if os.path.isdir(f"50-Annotation/{prefix}/annotation"):
                                logger.error(f"File \"{a}\" is missing or empty.")
                                pass
                    
                    # Improved Genome
                    if "avd" in globals():
                        if avd == "Pass":
                            # Create folder
                            if not os.path.isdir("60-Genomes/Improved"): os.makedirs("60-Genomes/Improved")
                            
                            # Move assembly to Improved Genome
                            os.system(f"cp {assembly_folder}/{prefix}/flye/assembly.fasta 60-Genomes/Improved/{prefix}.assembly.fasta")

                            # Coverage
                            if run_coverage == True:
                                # Logger
                                logger.info('**** Running coverage Draft ****')

                                # Call
                                coverage(assembly_folder=assembly_folder, prefix=prefix, software=software, reads_1=read_corr_1, reads2=read_corr_2, threads=threads, cpath=cpath)

                        else:
                            logger.warning(f"Sample {prefix} failed the draft annotation validation. Sample will be drop and no genome will be reported.")
                    else:
                        logger.error(f"Draft annotation validation for sample {prefix} was not completed. Please, check the logs.")
                else:
                    # Validation failed (Draft)
                    logger.warning(f"Sample {prefix} failed the assembly validation at the nanopore draft level (after failing the hybrid assembly validation). Sample will be drop and no genome will be reported.")
            else:
                # No validation variables (Draft)
                logger.error(f"Assembly validation for sample {prefix} was not completed. Please, check the logs.")
    else:
        # No validation variables (HA)
        logger.error(f"Validation for sample {prefix} was not completed. Please, check the logs.")



