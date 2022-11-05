# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  	: Adrián Gómez Repollés
# Email			: adrian.gomez@mbg.au.dk
# Created Date	: 25/03/2022
# version 		: '1.0'
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
import statistics as st # Python v3.4 or above
from gwf import Workflow, AnonymousTarget

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
# Example: 00-Data/TA1/SRR11719982_pass.fastq.gz 00-Data/TA1/SRR3927460_pass_1.fastq.gz 00-Data/TA1/SRR3927460_pass_2.fastq.gz TA1
# Open file:
file = "strains.tsv"

try:
	f = open(file, 'r')
except:
	print(f"{bcolors.FAIL}Error: {file} summary file not found.{bcolors.ENDC}")
	sys.exit()

# First Loop
for row in f:

	if row[0] == "#":
		continue

	# Divide line
	row = row.strip().split("\t")

	# Variables
	folder = row[3]
	illumina_corr_1 = row[1].split("/")[-1].replace(".gz","") + ".00.0_0.cor.fastq.gz"
	illumina_corr_2 = row[2].split("/")[-1].replace(".gz","") + ".00.0_0.cor.fastq.gz"

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

# ContigTable
gwf.target('ContigTable', inputs=[file], outputs=['ContigTable.tsv']) << """
# Remove file if existed
if [ -f ContigTable.tsv ]; then
		rm -rf ContigTable.tsv
	fi

# Create ContigTable
for iso in $(cat strains.tsv | cut -f 4); do
	# Get number of contigs
	cor1=$(cat strains.tsv | grep -w ${iso} | cut -f 2)
	cor2=$(cat strains.tsv | grep -w ${iso} | cut -f 2)
	drcnt=$(grep "^>" 20-Assembly/${iso}/flye/assembly.fasta | wc -l)
	hacnt=$(grep "^>" 30-HybridAssembly/${iso}/unicycler/assembly.fasta | wc -l)
	# Categorize contig comparison
	if [ ${drcnt} -lt ${hacnt} ]; then
		comp="Below"
	elif [ ${drcnt} -gt ${hacnt} ]; then
		comp="Above"
	else
		comp="Equal"
	fi
	# Create file
	printf "${iso}\t${cor1}\t${cor2}\t${drcnt}\t${hacnt}\t${comp}\n" >> ContigTable.tsv
done
"""

# Second Loop
try:
	f = open("ContigTable.tsv", 'r')
except:
	print(f"{bcolors.WARNING}WARNING: summary file with assembly contigs not found, run \"ContigTable\" after the hybrid assemblies.{bcolors.ENDC}\n")
	#sys.exit()

for row in f:

	if row[0] == "#":
		continue

	# Divide line
	row = row.strip().split("\t")

	# Variables
	folder = row[0]
	comp = row[5].strip()
	illumina_corr_1 = row[1].split("/")[-1].replace(".gz","") + ".00.0_0.cor.fastq.gz"
	illumina_corr_2 = row[2].split("/")[-1].replace(".gz","") + ".00.0_0.cor.fastq.gz"

	if comp == "Equal":
		# 30 Quast
		gwf.target_from_template("{}_30_quast_hybridassembly".format(folder), validation.quast(assembly="30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), out_dir="30-HybridAssembly/{}/Quast".format(folder), threads=4, memory=32))

		# 30 Alignment
		gwf.target_from_template("{}_30_align_illumina".format(folder), coverage.align_illumina(assembly="30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), illumina_corr_1="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_1), illumina_corr_2="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_2), out_dir="30-HybridAssembly/{}/Align".format(folder), threads=4, memory=16))
		gwf.target_from_template("{}_30_align_nanopore".format(folder), coverage.align_nanopore(assembly="30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), nanopore_corr="10-Correction/{}/Nanopore/nanopore.corrected.fasta".format(folder), out_dir="30-HybridAssembly/{}/Align".format(folder), threads=4, memory=16))

		# 30 Coverage
		gwf.target_from_template("{}_30_coverage".format(folder), coverage.coverage(in_dir="30-HybridAssembly/{}/Align".format(folder), out_dir="30-HybridAssembly/{}/Coverage".format(folder), memory=8))

		# 30 Coverage Plot
		gwf.target_from_template("{}_30_plot_coverage".format(folder), coverage.plot_coverage(in_dir="30-HybridAssembly/{}/Coverage".format(folder), out_dir="30-HybridAssembly/{}/Coverage/CovPlots".format(folder), memory=8, folder=folder))

		# 40 Validation
		database = busco_dict[LjTaxa[folder]]
		if database != "NA":
			gwf.target_from_template("{}_40_validation_hybridassembly".format(folder), validation.assembly_validation(assembly="30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), database=database, out_dir="40-Validation/{}/unicycler".format(folder), threads=4, memory=24))

			# Validate BUSCO
			try:
				b = [i for i in os.listdir("40-Validation/{}/unicycler/Busco".format(folder)) if i[-5:] == ".json"][0]
				bf = open("40-Validation/{}/unicycler/Busco/{}".format(folder,b))
				bd = json.load(bf)
				bv = validation.validate_busco(bd)
				bf.close()
			except:
				if os.path.isdir("40-Validation/{}/unicycler/Busco".format(folder)):
					print(f"{bcolors.FAIL}Error: 40-Validation/{folder}/unicycler/Busco summary file is missing or empty.{bcolors.ENDC}\n")
					continue

			# Validate CheckM
			c = "40-Validation/{}/unicycler/CheckM/results.tsv"
			try:
				cv = validation.validate_checkm(c)
			except:
				if os.path.isdir("40-Validation/{}/unicycler/CheckM".format(folder)):
					print(f"{bcolors.FAIL}Error: {c} is missing or empty.{bcolors.ENDC}\n")
					continue

		# Check validation output
		if "bv" in globals() and "cv" in globals():
			if bv == "Pass" and cv == "Pass":
				# 50 Annotation
				out_dir_yaml = "30-HybridAssembly/{}/unicycler/{}".format(folder, folder)
				genus = LjGenus[folder]
				if not os.path.exists(out_dir_yaml + ".submol.yml") and not os.path.exists(out_dir_yaml + ".input.yml"):
					annotation.pgap_files_creator(genus = genus, assembly = "30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), out_dir = out_dir_yaml)
					
				if os.path.exists(out_dir_yaml + '.submol.yml') and os.path.exists(out_dir_yaml + '.input.yml') and genus != "NA":
					gwf.target_from_template("{}_50_annotation_hybridassembly".format(folder), annotation.annotation(assembly="30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), input_yaml = "{}.input.yml".format(out_dir_yaml), out_dir="50-Annotation/{}/annotation".format(folder), threads=1, memory=4))

				# Check annotation output
				a = "50-Annotation/{}/annotation/annot.gff"
				try:
					av = validation.validate_pgap(a)
					if av == "Pass":
						# Move assembly to Complete Genome
						continue
				except:
					print(f"{bcolors.FAIL}Error: {c} is missing or empty.{bcolors.ENDC}")
					continue

			elif bv == "Failed" and cv == "Failed":
				# New 16S identification! (Here is where SyFi could be used!)
				#	1. Extract 16S sequences
				#	2. Remove duplicates
				#	3. Blast online
				#	4. Compare pre and new taxonomy (update pre to new)
				#	4.1. If equal: Stop
				#	4.2. If different: Re-run Busco and Annotation
				continue

	elif comp == "Below":
		# 20 Quast
		gwf.target_from_template("{}_20_quast_assembly".format(folder), validation.quast(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), out_dir="20-Assembly/{}/flye/Quast".format(folder), threads=4, memory=32))

		# 20 Alignment
		gwf.target_from_template("{}_20_align_illumina".format(folder), coverage.align_illumina(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), illumina_corr_1="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_1), illumina_corr_2="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_2), out_dir="20-Assembly/{}/Align".format(folder), threads=4, memory=16))
		gwf.target_from_template("{}_20_align_nanopore".format(folder), coverage.align_nanopore(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), nanopore_corr="10-Correction/{}/Nanopore/nanopore.corrected.fasta".format(folder), out_dir="20-Assembly/{}/Align".format(folder), threads=4, memory=16))

		# 20 Coverage
		gwf.target_from_template("{}_20_coverage".format(folder), coverage.coverage(in_dir="20-Assembly/{}/Align".format(folder), out_dir="20-Assembly/{}/Coverage".format(folder), memory=8))

		# 20 Coverage Plot
		gwf.target_from_template("{}_20_plot_coverage".format(folder), coverage.plot_coverage(in_dir="20-Assembly/{}/Coverage".format(folder), out_dir="20-Assembly/{}/Coverage/CovPlots".format(folder), memory=8, folder=folder))

		# 40 Validation
		database = busco_dict[LjTaxa[folder]]
		if database != "NA":
			gwf.target_from_template("{}_40_validation_assembly".format(folder), validation.assembly_validation(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), database=database, out_dir="40-Validation/{}/flye".format(folder), threads=4, memory=24))

			# Validate BUSCO
			try:
				b = [i for i in os.listdir("40-Validation/{}/flye/Busco".format(folder)) if i[-5:] == ".json"][0]
				bf = open("40-Validation/{}/flye/Busco/{}".format(folder,b))
				bd = json.load(bf)
				bv = validation.validate_busco(bd)
				bf.close()
			except:
				if os.path.isdir("40-Validation/{}/flye/Busco".format(folder)):
					print(f"{bcolors.FAIL}Error: 40-Validation/{folder}/flye/Busco summary file is missing or empty.{bcolors.ENDC}\n")
					continue

			# Validate CheckM
			c = "40-Validation/{}/flye/CheckM/results.tsv"
			try:
				cv = validation.validate_checkm(c)
			except:
				if os.path.isdir("40-Validation/{}/flye/CheckM".format(folder)):
					print(f"{bcolors.FAIL}Error: {c} is missing or empty.{bcolors.ENDC}\n")
					continue

		# Check validation output
		if "bv" in globals() and "cv" in globals():
			if bv == "Pass" and cv == "Pass":
				# 50 Annotation
				out_dir_yaml = "20-Assembly/{}/flye/{}".format(folder, folder)
				genus = LjGenus[folder]
				if not os.path.exists(out_dir_yaml + ".submol.yml") and not os.path.exists(out_dir_yaml + ".input.yml"):
					annotation.pgap_files_creator(genus = genus, assembly = "20-Assembly/{}/flye/assembly.fasta".format(folder), out_dir = out_dir_yaml)
				
				if os.path.exists(out_dir_yaml + '.submol.yml') and os.path.exists(out_dir_yaml + '.input.yml') and genus != "NA":
					gwf.target_from_template("{}_50_annotation_assembly".format(folder), annotation.annotation(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), input_yaml = "{}.input.yml".format(out_dir_yaml), out_dir="50-Annotation/{}/annotation".format(folder), threads=1, memory=4))

				# Check annotation output
				a = "50-Annotation/{}/annotation/annot.gff"
				try:
					av = validation.validate_pgap(a)
					if av == "Pass":
						# Move assembly to Improved Genome
						continue
				except:
					print(f"{bcolors.FAIL}Error: {c} is missing or empty.{bcolors.ENDC}")
					continue

			elif bv == "Failed" and cv == "Failed":
				continue

# SummaryTable
gwf.target('SummaryTableCompleteGenomes', inputs=[file], outputs=['SummaryTableCompleteGenomes.tsv']) << """
echo -n "Hello " > greeting.txt
cat name.txt >> greeting.txt
"""