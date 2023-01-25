# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  	: Adrián Gómez Repollés
# Email			: adrian.gomez@mbg.au.dk
# Created Date	: 24/01/2023
# version 		: '1.2'
# ---------------------------------------------------------------------------
""" Pipeline (gwf) to perform hybrid assembly of nanopore and illumina sequen
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
from src import annotation

# External
import os
import sys
import json
# import statistics as st # Python v3.4 or above
from gwf import Workflow

gwf = Workflow(defaults={"account": "CCRP_Data"})


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

# Database
#---------

file_path = "LjSphere_taxonomy.csv"

# Create dict from csv
def create_dict(path, col):
	d = {}
	f = open(path, "r")
	for l in f:
		l = l.split("\t")
		d[l[0]] = l[col].strip()
	return d

LjTaxa = create_dict(file_path, 4)
LjGenus = create_dict(file_path, 6)

# Busco databases
busco_dict = {'Actinomycetales': 'actinobacteria_class_odb10', 'Flavobacteriales': 'flavobacteriales_odb10', 'Bacillales': 'bacillales_odb10', 'Burkholderiales': 'burkholderiales_odb10', 'Caulobacterales': 'alphaproteobacteria_odb10', 'Rhizobiales': 'rhizobiales_odb10', 'Sphingomonadales': 'sphingomonadales_odb10', 'Pseudomonadales': 'pseudomonadales_odb10', 'Xanthomonadales': 'xanthomonadales_odb10', 'NA': 'NA'}

# Execution
#----------

## Strains/Isolates (nanopore illumina_1 illumina_2 folder_name)
# Example: 00-Data/TA1/SRR11719982_pass.fastq.gz 00-Data/TA1/SRR3927460_pass_1.fastq.gz 00-Data/TA1/SRR3927460_pass_2.fastq.gz TA1 False
# Open file:
file = "strains.tsv"

try:
	f = open(file, 'r')
except:
	print(f"{bcolors.FAIL}Error: {file} summary file not found.{bcolors.ENDC}")
	sys.exit()

# Log folder
if not os.path.isdir(".logs"): os.makedirs(".logs")

# First Loop
for row in f:

	if row[0] == "#":
		continue

	# Divide line
	row = row.strip().split("\t")

	# Variables
	folder = row[3]
	run_coverage = row[4]
	illumina_corr_1 = row[1].split("/")[-1].replace(".gz","") + ".00.0_0.cor.fastq.gz"
	illumina_corr_2 = row[2].split("/")[-1].replace(".gz","") + ".00.0_0.cor.fastq.gz"

	# Log file
	if os.path.exists(f".logs/{folder}.log"): os.remove(f".logs/{folder}.log")
	log = open(f".logs/{folder}.log", "w")
	log.write(f"Sample: {folder}\n")

	# 01 QC
	gwf.target_from_template('{}_01_qc_illumina'.format(folder), qc.qc_illumina(illumina_1=row[1], illumina_2=row[2], out_dir="01-QC/{}/Illumina".format(folder), threads=4, memory=8, folder=folder))
	gwf.target_from_template('{}_01_qc_nanopore'.format(folder), qc.qc_nanopore(nanopore=row[0], out_dir="01-QC/{}/Nanopore".format(folder), threads=4, memory=8, folder=folder))

	# 10 Correction
	gwf.target_from_template('{}_10_correct_illumina'.format(folder), correction.correct_illumina(illumina_1=row[1], illumina_2=row[2], illumina_corr_1=illumina_corr_1, illumina_corr_2=illumina_corr_2, out_dir="10-Correction/{}/Illumina".format(folder), threads=4, memory=24, folder=folder))
	gwf.target_from_template('{}_10_correct_nanopore'.format(folder), correction.correct_nanopore(nanopore=row[0], illumina_corr="10-Correction/{}/Illumina/corrected/illumina.corrected.fastq.gz".format(folder), out_dir="10-Correction/{}/Nanopore".format(folder), threads=4, memory=24, folder=folder))

	# 20 Assembly
	gwf.target_from_template('{}_20_assembly_flye'.format(folder), assembly.flye_assembly(nanopore_corr="10-Correction/{}/Nanopore/nanopore.corrected.fasta".format(folder), out_dir="20-Assembly/{}/flye".format(folder), threads=8, memory=64, folder=folder))

	# 30 Unicycler
	gwf.target_from_template("{}_30_hybrid_assembly_unicycler".format(folder), assembly.unicycler(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), illumina_corr_1="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_1), illumina_corr_2="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_2), nanopore_corr="10-Correction/{}/Nanopore/nanopore.corrected.fasta".format(folder), out_dir="30-HybridAssembly/{}/unicycler".format(folder), threads=8, memory=64, folder=folder))

	# 30 Quast
	gwf.target_from_template("{}_30_quast_hybridassembly".format(folder), validation.quast(assembly="30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), out_dir="30-HybridAssembly/{}/Quast".format(folder), threads=4, memory=32))

	# 30 Coverage
	if run_coverage == "True":
		# Log file
		log.write(f"Performing coverage for {folder}\n")

		# Alignment Illumina
		gwf.target_from_template("{}_30_align_illumina".format(folder), coverage.align_illumina(assembly="30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), illumina_corr_1="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_1), illumina_corr_2="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_2), out_dir="30-HybridAssembly/{}/Align".format(folder), threads=4, memory=16))

		# Alignment Nanopore
		gwf.target_from_template("{}_30_align_nanopore".format(folder), coverage.align_nanopore(assembly="30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), nanopore_corr="10-Correction/{}/Nanopore/nanopore.corrected.fasta".format(folder), out_dir="30-HybridAssembly/{}/Align".format(folder), threads=4, memory=16))

		# Coverage
		gwf.target_from_template("{}_30_coverage".format(folder), coverage.coverage(in_dir="30-HybridAssembly/{}/Align".format(folder), out_dir="30-HybridAssembly/{}/Coverage".format(folder), memory=8))

		# Coverage Plot
		gwf.target_from_template("{}_30_plot_coverage".format(folder), coverage.plot_coverage(in_dir="30-HybridAssembly/{}/Coverage".format(folder), out_dir="30-HybridAssembly/{}/Coverage/CovPlots".format(folder), memory=8, folder=folder))

	# 40 Validation
	database = busco_dict[LjTaxa[folder]]
	if database != "NA":
		gwf.target_from_template("{}_40_validation_hybridassembly".format(folder), validation.assembly_validation(assembly="30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), database=database, out_dir="40-Validation/{}/unicycler".format(folder), threads=4, memory=24))
	else:
		log.write(f"Validation not performed for {folder}. Taxonomy is equal to \"NA\"\n")

# ---------------------------------------------------------------
	if not os.path.exists(f"60-Genomes/Complete/{folder}.assembly.fasta"):
		# Dict
		val_ha = {"bv": "NA", "cv": "NA", "av": "NA"}

		# Validate BUSCO
		if os.path.isdir("40-Validation/{}/unicycler/Busco".format(folder)):
			try:
				b = [i for i in os.listdir("40-Validation/{}/unicycler/Busco".format(folder)) if i[-5:] == ".json"][0]
				bf = open("40-Validation/{}/unicycler/Busco/{}".format(folder,b))
				bd = json.load(bf)
				val_ha["bv"] = validation.validate_busco(bd)
				bf.close()
			except:
				log.write(f"{bcolors.FAIL}Error: 40-Validation/{folder}/unicycler/Busco summary file is missing or empty.{bcolors.ENDC}\n")

		# Validate CheckM
		if os.path.isdir("40-Validation/{}/unicycler/CheckM".format(folder)):
			c = "40-Validation/{}/unicycler/CheckM/results.tsv".format(folder)
			try:
				val_ha["cv"] = validation.validate_checkm(c)
			except:
				log.write(f"{bcolors.FAIL}Error: {c} is missing or empty.{bcolors.ENDC}\n")

	    # ---------------------------------------------------------------

		# Check validation output
		if val_ha["bv"] == "Pass" and val_ha["cv"] == "Pass":
			# 50 Annotation
			if os.path.isdir("30-HybridAssembly/{}/unicycler".format(folder)):
				out_dir_yaml = "30-HybridAssembly/{}/unicycler/{}".format(folder, folder)
				genus = LjGenus[folder]
				if not os.path.exists(out_dir_yaml + ".submol.yml") and not os.path.exists(out_dir_yaml + ".input.yml"):
					annotation.pgap_files_creator(genus = genus, assembly = "30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), out_dir = out_dir_yaml)
				
				if os.path.exists(out_dir_yaml + '.submol.yml') and os.path.exists(out_dir_yaml + '.input.yml') and genus != "NA":
					gwf.target_from_template("{}_50_annotation_hybridassembly".format(folder), annotation.annotation(assembly="30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), input_yaml = "{}.input.yml".format(out_dir_yaml), out_dir="50-Annotation/{}/unicycler".format(folder), threads=1, memory=4))
				else:
					log.write(f"PGAP input files might be missing or the defined genus is \"NA\" for sample {folder}\n")

			# Check annotation output
			if os.path.isdir("50-Annotation/{}/unicycler/annotation".format(folder)):
				a = "50-Annotation/{}/unicycler/annotation/annot.gff".format(folder)
				try:
					val_ha["av"] = validation.validate_pgap(a)
				except:
					if os.path.isdir("50-Annotation/{}/annotation".format(folder)):
						log.write(f"{bcolors.FAIL}Error: {a} is missing or empty.{bcolors.ENDC}\n")

			# Complete Genome
			if val_ha["av"] == "Pass":
				# Create folder
				if not os.path.isdir("60-Genomes/Complete"): os.makedirs("60-Genomes/Complete")

				# Move assembly to Complete Genome
				os.system("cp 30-HybridAssembly/{}/unicycler/assembly.fasta 60-Genomes/Complete/{}.assembly.fasta".format(folder, folder))

				# Remove assembly from Improved Genome
				if os.path.exists("60-Genomes/Improved/{}.assembly.fasta".format(folder)): os.remove("60-Genomes/Improved/{}.assembly.fasta".format(folder))

				# Log file
				log.write(f"Genome for {folder} was marked as Completed\n")
	else:
		log.write(f"Genome for {folder} was already done (Completed)\n")
		continue

# ---------------------------------------------------------------

	# If anything fails, the draft is taken into account
	if "Failed" in val_ha.values():

		# 20 Quast
		gwf.target_from_template("{}_20_quast_assembly".format(folder), validation.quast(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), out_dir="20-Assembly/{}/Quast".format(folder), threads=4, memory=32))

		# 20 Coverage
		if run_coverage == "True":
			# Log file
			log.write(f"Performing coverage for {folder}\n")

			# Alignment Illumina
			gwf.target_from_template("{}_20_align_illumina".format(folder), coverage.align_illumina(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), illumina_corr_1="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_1), illumina_corr_2="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_2), out_dir="20-Assembly/{}/Align".format(folder), threads=4, memory=16))

			# Alignment Nanopore
			gwf.target_from_template("{}_20_align_nanopore".format(folder), coverage.align_nanopore(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), nanopore_corr="10-Correction/{}/Nanopore/nanopore.corrected.fasta".format(folder), out_dir="20-Assembly/{}/Align".format(folder), threads=4, memory=16))

			# Coverage
			gwf.target_from_template("{}_20_coverage".format(folder), coverage.coverage(in_dir="20-Assembly/{}/Align".format(folder), out_dir="20-Assembly/{}/Coverage".format(folder), memory=8))

			# Coverage Plot
			gwf.target_from_template("{}_20_plot_coverage".format(folder), coverage.plot_coverage(in_dir="20-Assembly/{}/Coverage".format(folder), out_dir="20-Assembly/{}/Coverage/CovPlots".format(folder), memory=8, folder=folder))

		# 40 Validation
		database = busco_dict[LjTaxa[folder]]
		if database != "NA":
			gwf.target_from_template("{}_40_validation_assembly".format(folder), validation.assembly_validation(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), database=database, out_dir="40-Validation/{}/flye".format(folder), threads=4, memory=24))
		else:
			log.write(f"Validation not performed for {folder}. Taxonomy is equal to \"NA\"\n")

# ---------------------------------------------------------------
	if not os.path.exists(f"60-Genomes/Complete/{folder}.assembly.fasta"):
		# Dict
		val_draft = {"bv": "NA", "cv": "NA", "av": "NA"}

		# Validate BUSCO
		if os.path.isdir("40-Validation/{}/flye/Busco".format(folder)):
			try:
				b = [i for i in os.listdir("40-Validation/{}/flye/Busco".format(folder)) if i[-5:] == ".json"][0]
				bf = open("40-Validation/{}/flye/Busco/{}".format(folder,b))
				bd = json.load(bf)
				val_draft["bv"] = validation.validate_busco(bd)
				bf.close()
			except:
				log.write(f"{bcolors.FAIL}Error: 40-Validation/{folder}/flye/Busco summary file is missing or empty.{bcolors.ENDC}\n")

		# Validate CheckM
		if os.path.isdir("40-Validation/{}/flye/CheckM".format(folder)):
			c = "40-Validation/{}/flye/CheckM/results.tsv".format(folder)
			try:
				val_draft["cv"] = validation.validate_checkm(c)
			except:
				log.write(f"{bcolors.FAIL}Error: {c} is missing or empty.{bcolors.ENDC}\n")

		# ---------------------------------------------------------------

		# Check validation output

		if val_draft["bv"] == "Pass" and val_draft["cv"] == "Pass":
			# 50 Annotation
			if os.path.isdir("20-Assembly/{}/flye".format(folder)):
				out_dir_yaml = "20-Assembly/{}/flye/{}".format(folder, folder)
				genus = LjGenus[folder]
				if not os.path.exists(out_dir_yaml + ".submol.yml") and not os.path.exists(out_dir_yaml + ".input.yml"):
					annotation.pgap_files_creator(genus = genus, assembly = "20-Assembly/{}/flye/assembly.fasta".format(folder), out_dir = out_dir_yaml)
				
				if os.path.exists(out_dir_yaml + '.submol.yml') and os.path.exists(out_dir_yaml + '.input.yml') and genus != "NA":
					gwf.target_from_template("{}_50_annotation_assembly".format(folder), annotation.annotation(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), input_yaml = "{}.input.yml".format(out_dir_yaml), out_dir="50-Annotation/{}/flye/".format(folder), threads=1, memory=4))
				else:
					log.write(f"PGAP input files might be missing or the defined genus is \"NA\" for sample {folder}\n")

			# Check annotation output
			if os.path.isdir("50-Annotation/{}/flye/annotation".format(folder)):
				a = "50-Annotation/{}/flye/annotation/annot.gff".format(folder)
				try:
					val_draft["av"] = validation.validate_pgap(a)
				except:
					log.write(f"{bcolors.FAIL}Error: {a} is missing or empty.{bcolors.ENDC}\n")

			# Improved Genome
			if val_draft["av"] == "Pass" and not os.path.exists(f"60-Genomes/Improved/{folder}.assembly.fasta"):
				# Create folder
				if not os.path.isdir("60-Genomes/Improved"): os.makedirs("60-Genomes/Improved")

				# Move assembly to Improved Genome
				os.system("cp 20-Assembly/{}/flye/assembly.fasta 60-Genomes/Improved/{}.assembly.fasta".format(folder, folder))

				# Log file
				log.write(f"Genome for {folder} was marked as Improved\n")
	else:
		log.write(f"Genome for {folder} was already done (Completed)\n")
		continue

# ---------------------------------------------------------------

	if val_draft["bv"] == "Failed" or val_draft["cv"] == "Failed":
		log.write(f"{bcolors.FAIL}Error: draft assembly validation failed for {folder}.{bcolors.ENDC}\n")

	# Close log file
	log.close()

# SummaryTable
gwf.target('SummaryTableCompleteGenomes', inputs=[file], outputs=['SummaryTableCompleteGenomes.tsv']) << """
echo -n "Hello " > greeting.txt
cat name.txt >> greeting.txt
"""