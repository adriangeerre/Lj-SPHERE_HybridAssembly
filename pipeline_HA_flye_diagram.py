# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  	: Adrián Gómez Repollés
# Email			: adrian.gomez@mbg.au.dk
# Created Date	: 25/03/2022
# version 		: '1.0'
# ---------------------------------------------------------------------------
""" Pipeline (gwf) to perform hybrid assembly of nanopore and illumina sequen
-cing data. It is based on Benjamin Perry pipeline and uses flye as the main
assembler. """ 
# ---------------------------------------------------------------------------
# Imports 
# ---------------------------------------------------------------------------

from genericpath import isdir
from html import entities
import os
import sys
import yaml
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

# Functions
#----------

# Quality Control
def qc_illumina(illumina_1, illumina_2, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("01-QC/" + folder + "/Illumina") == False: os.makedirs("01-QC/" + folder + "/Illumina")

	# Names
	fqc_illumina_1 = illumina_1.split("/")[-1].split(".")[0] + "_fastqc"
	fqc_illumina_2 = illumina_2.split("/")[-1].split(".")[0] + "_fastqc"

	# GWF
	inputs = ["{}".format(illumina_1), "{}".format(illumina_2)]
	outputs = ["{}/{}.zip".format(out_dir, fqc_illumina_1), "{}/{}.html".format(out_dir, fqc_illumina_1), "{}/{}.zip".format(out_dir, fqc_illumina_2), "{}/{}.html".format(out_dir, fqc_illumina_2)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '2:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# FastQC
	conda activate ha-flye
	fastqc -t {threads} -o {out_dir} {illumina_1} {illumina_2}
	'''.format(illumina_1=illumina_1, illumina_2=illumina_2, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def qc_nanopore(nanopore, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("01-QC/" + folder) == False: os.makedirs("01-QC/" + folder)

	# GWF
	inputs = ["{}".format(nanopore)]
	outputs = ["{}/{}-NanoPlot-report.html".format(out_dir, nanopore.split("/")[-1].replace(".fastq.gz",""), )]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '4:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# NanoPlot
	conda activate NanoQC
	NanoPlot -o {out_dir} -p {prefix} --info_in_report --N50 --title {title} --fastq {nanopore} --threads {threads}
	'''.format(nanopore=nanopore, out_dir=out_dir, prefix="{}-".format(nanopore.split("/")[-1].replace(".fastq.gz","")), title=nanopore.split("/")[-1], threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)	

# Correction
def correct_illumina(illumina_1, illumina_2, illumina_corr_1, illumina_corr_2, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("10-Correction/" + folder) == False: os.makedirs("10-Correction/" + folder)

	# GWF
	inputs = ["{}".format(illumina_1), "{}".format(illumina_2)]
	outputs = ["{}/corrected/illumina.corrected.fastq.gz".format(out_dir), "{}/corrected/{}".format(out_dir, illumina_corr_1), "{}/corrected/{}".format(out_dir, illumina_corr_2)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '24:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# SPAdes
	conda activate SPAdes
	spades.py -1 {illumina_1} -2 {illumina_2} -o {out_dir} --only-error-correction -t {threads} -m {memory}
	cat {out_dir}/corrected/{illumina_corr_1} {out_dir}/corrected/{illumina_corr_2} > {out_dir}/corrected/illumina.corrected.fastq.gz
	'''.format(illumina_1=illumina_1, illumina_2=illumina_2, illumina_corr_1=illumina_corr_1, illumina_corr_2=illumina_corr_2, out_dir=out_dir, threads=threads, memory=memory)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)	

def correct_nanopore(nanopore, illumina_corr, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("10-Correction/" + folder + "/Nanopore") == False: os.makedirs("10-Correction/" + folder + "/Nanopore")

	# GWF
	inputs = ["{}".format(nanopore), "{}".format(illumina_corr)]
	outputs = ["{}/nanopore.corrected.fasta".format(out_dir)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '24:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# LoRDEC
	conda activate ha-flye
	lordec-correct -i {nanopore} -2 {illumina_corr} -k 19 -s 4 -T {threads} -p -o {out_dir}/nanopore.kmer19.fasta
	lordec-correct -i {out_dir}/nanopore.kmer19.fasta -2 {illumina_corr} -k 31 -s 3 -T {threads} -p -o {out_dir}/nanopore.kmer31.fasta
	lordec-correct -i {out_dir}/nanopore.kmer31.fasta -2 {illumina_corr} -k 41 -s 3 -T {threads} -p -o {out_dir}/nanopore.corrected.fasta
	'''.format(nanopore=nanopore, illumina_corr=illumina_corr, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Assembler
def flye_assembly(nanopore_corr, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("20-Assembly/" + folder) == False: os.makedirs("20-Assembly/" + folder)

	# GWF
	inputs = ["{}".format(nanopore_corr)]
	outputs = ["{}/assembly.fasta".format(out_dir)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '12:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Flye
	conda activate ha-flye
	flye --nano-corr {nanopore_corr} --plasmids --out-dir {out_dir} --threads {threads}
	'''.format(nanopore_corr=nanopore_corr, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Unicycler
def unicycler(assembly, illumina_corr_1, illumina_corr_2, nanopore_corr, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("30-HybridAssembly/" + folder) == False: os.makedirs("30-HybridAssembly/" + folder)

	# GWF
	inputs = ["{}".format(assembly), "{}".format(illumina_corr_1), "{}".format(illumina_corr_2), "{}".format(nanopore_corr)]
	outputs = ["{}/assembly.fasta".format(out_dir)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '12:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Unicycler
	conda activate Unicycler
	unicycler -1 {illumina_corr_1} -2 {illumina_corr_2} --existing_long_read_assembly {assembly} -l {nanopore_corr} --threads {threads} --keep 2 --verbosity 2 -o {out_dir}
	'''.format(assembly=assembly, illumina_corr_1=illumina_corr_1, illumina_corr_2=illumina_corr_2, nanopore_corr=nanopore_corr, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Quast
def quast(assembly, out_dir, threads, memory):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# GWF
	inputs = ["{}".format(assembly)]
	outputs = ["{}/report.tsv".format(out_dir)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '12:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Quast
	conda activate Quast
	quast -o {out_dir} -t {threads} {assembly}
	'''.format(assembly=assembly, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Alignment
def align_illumina(assembly, illumina_corr_1, illumina_corr_2, out_dir, threads, memory):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# GWF
	inputs = ["{}".format(illumina_corr_1), "{}".format(illumina_corr_2), "{}".format(assembly)]
	outputs = ["{}/Illumina.sort.bam".format(out_dir), "{}/Illumina.sort.bai".format(out_dir)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '4:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Bowtie2 + Samtools
	conda activate ha-flye
	mkdir -p {in_dir}/index_contigs
	bowtie2-build --threads {threads} {assembly} {in_dir}/index_contigs/index
	bowtie2 -x {in_dir}/index_contigs/index -1 {illumina_corr_1} -2 {illumina_corr_2} -S {out_dir}/Illumina.sam --threads {threads} 2> {out_dir}/bowtie2.log

	samtools view -bS {out_dir}/Illumina.sam -@ {threads} > {out_dir}/Illumina.bam
	samtools sort -o {out_dir}/Illumina.sort.bam -O bam {out_dir}/Illumina.bam -@ {threads}
	samtools index -b {out_dir}/Illumina.sort.bam {out_dir}/Illumina.sort.bai -@ {threads}
	rm {out_dir}/Illumina.sam {out_dir}/Illumina.bam
	'''.format(assembly=assembly, illumina_corr_1=illumina_corr_1, illumina_corr_2=illumina_corr_2, in_dir='/'.join(assembly.split("/")[:-1]), out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def align_nanopore(assembly, nanopore_corr, out_dir, threads, memory):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# GWF
	inputs = ["{}".format(nanopore_corr), "{}".format(assembly)]
	outputs = ["{}/Nanopore.sort.bam".format(out_dir), "{}/Nanopore.sort.bai".format(out_dir)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '4:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Minimap2 + Samtools
	conda activate ha-flye
	minimap2 -a -o {out_dir}/Nanopore.sam -t {threads} -x map-ont {assembly} {nanopore_corr}

	samtools view -bS {out_dir}/Nanopore.sam -@ {threads} > {out_dir}/Nanopore.bam
	samtools sort -o {out_dir}/Nanopore.sort.bam -O bam {out_dir}/Nanopore.bam -@ {threads}
	samtools index -b {out_dir}/Nanopore.sort.bam {out_dir}/Nanopore.sort.bai -@ {threads}
	rm {out_dir}/Nanopore.sam {out_dir}/Nanopore.bam
	'''.format(assembly=assembly, nanopore_corr=nanopore_corr, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Coverage
def coverage(in_dir, out_dir, memory):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# GWF
	inputs = ["{}/Illumina.sort.bam".format(in_dir), "{}/Nanopore.sort.bam".format(in_dir)]
	outputs = ["{}/Illumina.cov".format(out_dir), "{}/Nanopore.cov".format(out_dir)]
	options = {'cores': 1,'memory': '{}g'.format(memory), 'queue': 'short', 'walltime': '2:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Bedtools
	conda activate Bedtools
	bedtools genomecov -d -ibam {in_dir}/Illumina.sort.bam > {out_dir}/Illumina.cov
	bedtools genomecov -d -ibam {in_dir}/Nanopore.sort.bam > {out_dir}/Nanopore.cov
	'''.format(in_dir=in_dir, out_dir=out_dir)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Plot Coverage
def plot_coverage(in_dir, out_dir, memory, folder):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# Number of contigs
	if os.path.exists("30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder)):
		num_contigs = len([1 for line in open("30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder)) if line.startswith(">")])
	else:
		num_contigs = 1

	# GWF
	inputs = ["{}/Illumina.cov".format(in_dir), "{}/Nanopore.cov".format(in_dir)]
	outputs = ["{}/{}.pdf".format(out_dir, num) for num in range(1,num_contigs+1) if os.path.exists("30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder))]
	options = {'cores': 1,'memory': '{}g'.format(memory), 'queue': 'short', 'walltime': '2:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# R
	conda activate Renv
	Rscript coverage.R -i {in_dir} -o {out_dir}
	'''.format(in_dir=in_dir, out_dir=out_dir)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Annotation
def annotation(assembly, input_yaml, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# GWF
	inputs = ["{}".format(assembly), "{}".format(input_yaml)]
	outputs = ["{}/annot.{}".format(out_dir, ext) for ext in ["faa","fna","gbk","gff","sqn"]]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'short', 'walltime': '12:00:00'}

	spec='''
	# Cache and tmp folders
	export SINGULARITY_CACHEDIR=/scratch/$SLURM_JOBID
	export SINGULARITY_TMPDIR=/scratch/$SLURM_JOBID

	# PGAP (already in path)
	python /home/agomez/programas/PGAP/pgap.py -d -n --no-internet --ignore-all-errors --docker singularity -o {out_dir} --memory {memory} --container-path ~/programas/SingularityImages/pgap_2022-08-11.build6275.sif {input_yaml}
	'''.format(input_yaml=input_yaml, out_dir=out_dir, folder=folder, memory=memory)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def assembly_validation(assembly, database, out_dir, threads, memory):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# GWF
	inputs = ["{}".format(assembly)]
	outputs = ["{}/Busco/short_summary.specific.{}.Busco.txt".format(out_dir,database), "{}/CheckM/results.tsv".format(out_dir)]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '2:00:00'}
	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# BUSCO
	conda activate Busco
	busco -m genome -i {assembly} -o {out_dir}/Busco -l {database} -f

	# CheckM
	conda activate CheckM
	checkm lineage_wf -x fasta -t {threads} {assembly_folder} {assembly} {out_dir}/CheckM --reduced_tree
	checkm qa {out_dir}/CheckM/lineage.ms -f {out_dir}/CheckM/results.tsv --tab_table -t {threads}
	'''.format(assembly=assembly, assembly_folder=assembly.strip("assembly.fasta"), database=database, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# def diamond(out_dir, diamond_db, threads, memory, folder):
# 	# Folder structure
# 	if os.path.isdir("80-Validation") == False: os.mkdir("80-Validation")
# 	if os.path.isdir("80-Validation/" + folder) == False: os.mkdir("80-Validation/" + folder)
# 	if os.path.isdir("80-Validation/" + folder + "/Diamond") == False: os.mkdir("80-Validation/" + folder + "/Diamond")

# 	# GWF
# 	inputs = ["70-Prokka/{}/prokka_{}.faa".format(folder, folder)]
# 	outputs = ["{}/{}.diamond.tsv".format(out_dir, folder)]
# 	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '6:00:00'}
# 	spec='''
# 	# Source conda to work with environments
# 	source ~/programas/minconda3.9/etc/profile.d/conda.sh

# 	# Diamond (db pre-computed)
# 	conda activate Diamond
# 	diamond blastp --query 70-Prokka/{folder}/prokka_{folder}.faa --db {diamond_db} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen sseq qseq --threads {threads} --out {out_dir}/{folder}.diamond.tsv --un {out_dir}/{folder}.diamond.unalign.fasta --unfmt fasta --fast --tmpdir /scratch/$SLURM_JOBID/ --parallel-tmpdir /scratch/$SLURM_JOBID/
# 	'''.format(out_dir=out_dir, diamond_db=diamond_db, folder=folder, threads=threads)

# 	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def fastani(hybrid_assembly, illumina_genomes, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("90-FastANI") == False: os.makedirs("90-FastANI")
	if os.path.isdir("90-FastANI/" + folder) == False: os.makedirs("90-FastANI/" + folder)

	# GWF
	inputs = ["{}".format(hybrid_assembly), "{}".format(illumina_genomes)]
	outputs = ["{}/{}.fastani.tsv".format(out_dir, folder)]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '00:30:00'}
	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# FastANI
	conda activate FastANI
	fastANI -q {hybrid_assembly} --rl {illumina_genomes} -o {out_dir}/{folder}.fastani.tsv -t {threads}
	'''.format(hybrid_assembly=hybrid_assembly, illumina_genomes=illumina_genomes, out_dir=out_dir, folder=folder, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Functions
#----------

def pgap_files_creator(genus, assembly, out_dir):
	# Submol
	dct_submol = {'topology':"'circular'", 'location':"'chromosome'", 'organism':{'genus_species':"'{}'".format(genus)}}
	submol = os.path.abspath(".") + "/" + out_dir + '.submol.yml'
	text = yaml.dump(dct_submol, sort_keys=False)
	text = text.replace("'''", "'") # PGAP only works if submol has single quotes for the values but not the keys.
	with open(submol, 'w') as yaml_file:
		yaml_file.write(text)

	# Input
	assembly = os.path.abspath(".") + "/" + assembly
	dct_input = {'fasta':{'class':'File', 'location':'{}'.format(assembly)}, 'submol':{'class':'File', 'location':'{}'.format(submol)}}
	input = os.path.abspath(".") + "/" + out_dir + '.input.yml'
	with open(input, 'w') as yaml_file:
		yaml.dump(dct_input, yaml_file)

# IMPLEMENT NO GENE OVERLAP! (AGAT)
def genbank_to_dict(gff):
	# Dict
	features = {}
	# Read file
	f = open(gff, "r")
	# Loop entries
	tlength = 0
	entries = {"gene": [], "CDS": [], "pseudogene": [], "rRNA_5S": [], "rRNA_16S": [], "rRNA_23S": [], "tRNA": []}
	for i in f:
		# Avoid header info
		if i[0] == "#":
			continue	
		i = i.split("\t")
		# Obtain total length
		if i[2] == "region":
			tlength += int(i[4])
		elif i[2] == "rRNA":
			p = i[-1].strip().split(";")
			p = [j for j in p if j[0:7] == "product"]
			p = p[0].split("=")[1].split(" ")[0]
			if p != "":
				entries["rRNA_{}".format(p)].append(abs(int(i[4])-int(i[3])))
		else:
			# Count defined features
			if i[2] in entries.keys():
				entries[i[2]].append(abs(int(i[4])-int(i[3])))
	# Include total length
	entries["tlength"] = [tlength]
	return entries

# Validations
def validate_busco(bd):
	# Dict
	r = {"c": 0, "s": 0, "d": 0, "f": 0, "m": 0}
	# Values
	if (bd["C"] / int(bd["dataset_total_buscos"])) * 100 > 90: r["c"] = 1
	if (bd["S"] / int(bd["dataset_total_buscos"])) * 100 > 90: r["s"] = 1
	if (bd["D"] / int(bd["dataset_total_buscos"])) * 100 < 5: r["d"] = 1
	if (bd["F"] / int(bd["dataset_total_buscos"])) * 100 < 5: r["f"] = 1
	if (bd["M"] / int(bd["dataset_total_buscos"])) * 100 < 5: r["m"] = 1
	# Check
	if sum(list(r.values())) == 5:
		return "Pass"
	else:
		return "Failed"

def validate_checkm(cm):
	# Dict
	r = {"Comp": 0, "Cont": 0}
	# Read file
	f = open(cm, "r")
	l = f.readlines()[1].strip().split()
	f.close()
	# Values
	if l[-3] > 90: r["Comp"] = 1
	if l[-3] < 5: r["Comp"] = 1
	# Check
	if sum(list(r.values())) == 2:
		return "Pass"
	else:
		return "Failed"

def validate_pgap(gff):
	# Dict
	r = {"GeneRatio": 0, "PseudoRatio": 0, "tRNA": 0, "rRNA": 0}
	ribo = {"rRNA_5S": 0, "rRNA_16S": 0, "rRNA_23S": 0}
	# Summary
	d = genbank_to_dict(gff)
	cd = {i: len(j) for i,j in d.items()} # Total genes = genes + pseudogene
	# Pre-knowledge (E. coli based!)
	l5S = 120
	l16S = 1542
	l23S = 2904
	# Ribosomal Value (at least one and median not deviating more than 10%)
	if cd["rRNA_5S"] >= 1 and abs(st.median(d["rRNA_5S"]) - l5S) <= l5S * 0.1: ribo["5S"] = 1
	if cd["rRNA_16S"] >= 1 and abs(st.median(d["rRNA_16S"]) - l16S) <= l16S * 0.1: ribo["16S"] = 1
	if cd["rRNA_23S"] >= 1 and abs(st.median(d["rRNA_23S"]) - l23S) <= l23S * 0.1: ribo["23S"] = 1
	# Values
	tgenes = cd["gene"] + cd["pseudogene"]
	gene_ratio = (tgenes * 1000) / d[tlength][0]
	if gene_ratio > 0.9: r["GeneRatio"] = 1
	if cd["pseudogene"] / tgenes < 0.2: r["PseudoRatio"] = 1
	if cd["tRNA"] > 20: r["tRNA"] = 1
	if sum(list(ribo.values())) == 3: r["rRNA"] = 1
	# Check
	if sum(list(r.values())) == 4:
		return "Pass"
	else:
		return "Failed"

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

# DiamondDB
# diamond_db = "/home/agomez/CCRP_Data/AGR/Data/UniProt/BacterialDB_diamond/uniprot_bacterial"

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
	gwf.target_from_template('{}_01_qc_illumina'.format(folder), qc_illumina(illumina_1=row[1], illumina_2=row[2], out_dir="01-QC/{}/Illumina".format(folder), threads=4, memory=8, folder=folder))
	gwf.target_from_template('{}_01_qc_nanopore'.format(folder), qc_nanopore(nanopore=row[0], out_dir="01-QC/{}/Nanopore".format(folder), threads=4, memory=8, folder=folder))

	# 10 Correction
	gwf.target_from_template('{}_10_correct_illumina'.format(folder), correct_illumina(illumina_1=row[1], illumina_2=row[2], illumina_corr_1=illumina_corr_1, illumina_corr_2=illumina_corr_2, out_dir="10-Correction/{}/Illumina".format(folder), threads=4, memory=24, folder=folder))
	gwf.target_from_template('{}_10_correct_nanopore'.format(folder), correct_nanopore(nanopore=row[0], illumina_corr="10-Correction/{}/Illumina/corrected/illumina.corrected.fastq.gz".format(folder), out_dir="10-Correction/{}/Nanopore".format(folder), threads=4, memory=24, folder=folder))

	# 20 Assembly
	gwf.target_from_template('{}_20_assembly_flye'.format(folder), flye_assembly(nanopore_corr="10-Correction/{}/Nanopore/nanopore.corrected.fasta".format(folder), out_dir="20-Assembly/{}/flye".format(folder), threads=8, memory=64, folder=folder))

	# 30 Unicycler
	gwf.target_from_template("{}_30_hybrid_assembly_unicycler".format(folder), unicycler(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), illumina_corr_1="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_1), illumina_corr_2="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_2), nanopore_corr="10-Correction/{}/Nanopore/nanopore.corrected.fasta".format(folder), out_dir="30-HybridAssembly/{}/unicycler".format(folder), threads=8, memory=64, folder=folder))

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
		gwf.target_from_template("{}_30_quast_hybridassembly".format(folder), quast(assembly="30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), out_dir="30-HybridAssembly/{}/Quast".format(folder), threads=4, memory=32))

		# 30 Alignment
		gwf.target_from_template("{}_30_align_illumina".format(folder), align_illumina(assembly="30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), illumina_corr_1="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_1), illumina_corr_2="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_2), out_dir="30-HybridAssembly/{}/Align".format(folder), threads=4, memory=16))
		gwf.target_from_template("{}_30_align_nanopore".format(folder), align_nanopore(assembly="30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), nanopore_corr="10-Correction/{}/Nanopore/nanopore.corrected.fasta".format(folder), out_dir="30-HybridAssembly/{}/Align".format(folder), threads=4, memory=16))

		# 30 Coverage
		gwf.target_from_template("{}_30_coverage".format(folder), coverage(in_dir="30-HybridAssembly/{}/Align".format(folder), out_dir="30-HybridAssembly/{}/Coverage".format(folder), memory=8))

		# 30 Coverage Plot
		gwf.target_from_template("{}_30_plot_coverage".format(folder), plot_coverage(in_dir="30-HybridAssembly/{}/Coverage".format(folder), out_dir="30-HybridAssembly/{}/Coverage/CovPlots".format(folder), memory=8, folder=folder))

		# 40 Validation
		database = busco_dict[LjTaxa[folder]]
		if database != "NA":
			gwf.target_from_template("{}_40_validation_hybridassembly".format(folder), assembly_validation(assembly="30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), database=database, out_dir="40-Validation/{}/unicycler/".format(folder), threads=4, memory=24))

			# Validate BUSCO
			try:
				b = [i for i in os.listdir("40-Validation/{}/unicycler/Busco").format(folder) if i[-5:] == ".json"][0]
				bf = open("40-Validation/{}/unicycler/Busco/{}".format(folder,b))
				bd = json.load(bf)
				bv = validate_busco(bd)
				bf.close()
			except:
				if os.path.isdir("40-Validation/{}/unicycler/Busco".format(folder)):
					print(f"{bcolors.FAIL}Error: 40-Validation/{folder}/unicycler/Busco summary file is missing or empty.{bcolors.ENDC}\n")
					continue

			# Validate CheckM
			c = "40-Validation/{}/unicycler/CheckM/results.tsv"
			try:
				cv = validate_checkm(c)
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
					pgap_files_creator(genus = genus, assembly = "30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), out_dir = out_dir_yaml)
					
				if os.path.exists(out_dir_yaml + '.submol.yml') and os.path.exists(out_dir_yaml + '.input.yml') and genus != "NA":
					gwf.target_from_template("{}_50_annotation_hybridassembly".format(folder), annotation(assembly="30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder), input_yaml = "{}.input.yml".format(out_dir_yaml), out_dir="50-Annotation/{}/annotation".format(folder), threads=1, memory=4))

				# Check annotation output
				a = "50-Annotation/{}/annotation/annot.gff"
				try:
					av = validate_pgap(a)
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
		gwf.target_from_template("{}_20_quast_assembly".format(folder), quast(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), out_dir="20-Assembly/{}/flye/Quast".format(folder), threads=4, memory=32))

		# 20 Alignment
		gwf.target_from_template("{}_20_align_illumina".format(folder), align_illumina(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), illumina_corr_1="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_1), illumina_corr_2="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_2), out_dir="20-Assembly/{}/Align".format(folder), threads=4, memory=16))
		gwf.target_from_template("{}_20_align_nanopore".format(folder), align_nanopore(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), nanopore_corr="10-Correction/{}/Nanopore/nanopore.corrected.fasta".format(folder), out_dir="20-Assembly/{}/Align".format(folder), threads=4, memory=16))

		# 20 Coverage
		gwf.target_from_template("{}_20_coverage".format(folder), coverage(in_dir="20-Assembly/{}/Align".format(folder), out_dir="20-Assembly/{}/Coverage".format(folder), memory=8))

		# 20 Coverage Plot
		gwf.target_from_template("{}_20_plot_coverage".format(folder), plot_coverage(in_dir="20-Assembly/{}/Coverage".format(folder), out_dir="20-Assembly/{}/Coverage/CovPlots".format(folder), memory=8, folder=folder))

		# 40 Validation
		database = busco_dict[LjTaxa[folder]]
		if database != "NA":
			gwf.target_from_template("{}_40_validation_assembly".format(folder), assembly_validation(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), database=database, out_dir="40-Validation/{}/flye".format(folder), threads=4, memory=24))

			# Validate BUSCO
			try:
				b = [i for i in os.listdir("40-Validation/{}/flye/Busco").format(folder) if i[-5:] == ".json"][0]
				bf = open("40-Validation/{}/flye/Busco/{}".format(folder,b))
				bd = json.load(bf)
				bv = validate_busco(bd)
				bf.close()
			except:
				if os.path.isdir("40-Validation/{}/flye/Busco".format(folder)):
					print(f"{bcolors.FAIL}Error: 40-Validation/{folder}/flye/Busco summary file is missing or empty.{bcolors.ENDC}\n")
					continue

			# Validate CheckM
			c = "40-Validation/{}/flye/CheckM/results.tsv"
			try:
				cv = validate_checkm(c)
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
					pgap_files_creator(genus = genus, assembly = "20-Assembly/{}/flye/assembly.fasta".format(folder), out_dir = out_dir_yaml)
				
				if os.path.exists(out_dir_yaml + '.submol.yml') and os.path.exists(out_dir_yaml + '.input.yml') and genus != "NA":
					gwf.target_from_template("{}_50_annotation_assembly".format(folder), annotation(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), input_yaml = "{}.input.yml".format(out_dir_yaml), out_dir="50-Annotation/{}/annotation".format(folder), threads=1, memory=4, folder=folder))

				# Check annotation output
				a = "50-Annotation/{}/annotation/annot.gff"
				try:
					av = validate_pgap(a)
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