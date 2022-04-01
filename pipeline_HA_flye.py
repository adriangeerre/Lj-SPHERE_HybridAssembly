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

import os
import sys
from datetime import datetime
from gwf import Workflow, AnonymousTarget

gwf = Workflow(defaults={"account": "CCRP_Data"})

# Functions
#----------

# Quality Control
def qc_illumina(illumina_1, illumina_2, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("01-QC") == False: os.mkdir("01-QC")
	if os.path.isdir("01-QC/" + folder) == False: os.mkdir("01-QC/" + folder)
	if os.path.isdir("01-QC/" + folder + "/Illumina") == False: os.mkdir("01-QC/" + folder + "/Illumina")

	# Names
	fqc_illumina_1 = illumina_1.split("/")[-1].split(".")[0] + "_fastqc"
	fqc_illumina_2 = illumina_2.split("/")[-1].split(".")[0] + "_fastqc"

	# GWF
	inputs = ["{}".format(illumina_1), "{}".format(illumina_2)]
	outputs = ["{}/{}.zip".format(out_dir, fqc_illumina_1), "{}/{}.html".format(out_dir, fqc_illumina_1), "{}/{}.zip".format(out_dir, fqc_illumina_2), "{}/{}.html".format(out_dir, fqc_illumina_2)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '2:00:00'}

	spec='''
	fastqc -t {threads} -o {out_dir} {illumina_1} {illumina_2}
	'''.format(illumina_1=illumina_1, illumina_2=illumina_2, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def qc_nanopore(nanopore, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("01-QC") == False: os.mkdir("01-QC")
	if os.path.isdir("01-QC/" + folder) == False: os.mkdir("01-QC/" + folder)

	# GWF
	inputs = ["{}".format(nanopore)]
	outputs = ["{}/{}-NanoPlot-report.html".format(out_dir, nanopore.split("/")[-1].replace(".fastq.gz",""), )]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '4:00:00'}

	spec='''
	NanoPlot -o {out_dir} -p {prefix} --info_in_report --N50 --title {title} --fastq {nanopore} --threads {threads}
	'''.format(nanopore=nanopore, out_dir=out_dir, prefix="{}-".format(nanopore.split("/")[-1].replace(".fastq.gz","")), title=nanopore.split("/")[-1], threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)	

# Correction
def correct_illumina(illumina_1, illumina_2, illumina_corr_1, illumina_corr_2, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("10-Correction") == False: os.mkdir("10-Correction")
	if os.path.isdir("10-Correction/" + folder) == False: os.mkdir("10-Correction/" + folder)

	# GWF
	inputs = ["{}".format(illumina_1), "{}".format(illumina_2)]
	outputs = ["{}/corrected/illumina.corrected.fastq.gz".format(out_dir), "{}/corrected/{}".format(out_dir, illumina_corr_1), "{}/corrected/{}".format(out_dir, illumina_corr_2)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '24:00:00'}

	spec='''
	spades.py -1 {illumina_1} -2 {illumina_2} -o {out_dir} --only-error-correction -t {threads} -m {memory}
	cat {out_dir}/corrected/{illumina_corr_1} {out_dir}/corrected/{illumina_corr_2} > {out_dir}/corrected/illumina.corrected.fastq.gz
	'''.format(illumina_1=illumina_1, illumina_2=illumina_2, illumina_corr_1=illumina_corr_1, illumina_corr_2=illumina_corr_2, out_dir=out_dir, threads=threads, memory=memory)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)	

def correct_nanopore(nanopore, illumina_corr, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("10-Correction") == False: os.mkdir("10-Correction")
	if os.path.isdir("10-Correction/" + folder) == False: os.mkdir("10-Correction/" + folder)
	if os.path.isdir("10-Correction/" + folder + "/Nanopore") == False: os.mkdir("10-Correction/" + folder + "/Nanopore")

	# GWF
	inputs = ["{}".format(nanopore), "{}".format(illumina_corr)]
	outputs = ["{}/nanopore.corrected.fasta".format(out_dir)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '24:00:00'}

	spec='''
	lordec-correct -i {nanopore} -2 {illumina_corr} -k 19 -s 4 -T {threads} -p -o {out_dir}/nanopore.kmer19.fasta
	lordec-correct -i {out_dir}/nanopore.kmer19.fasta -2 {illumina_corr} -k 31 -s 3 -T {threads} -p -o {out_dir}/nanopore.kmer31.fasta
	lordec-correct -i {out_dir}/nanopore.kmer31.fasta -2 {illumina_corr} -k 41 -s 3 -T {threads} -p -o {out_dir}/nanopore.corrected.fasta
	'''.format(nanopore=nanopore, illumina_corr=illumina_corr, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Assembler
def flye_assembly(nanopore_corr, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("20-Assembly") == False: os.mkdir("20-Assembly")
	if os.path.isdir("20-Assembly/" + folder) == False: os.mkdir("20-Assembly/" + folder)

	# GWF
	inputs = ["{}".format(nanopore_corr)]
	outputs = ["{}/assembly.fasta".format(out_dir)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '12:00:00'}

	spec='''
	flye --nano-corr {nanopore_corr} --plasmids --out-dir {out_dir} --threads {threads}
	'''.format(nanopore_corr=nanopore_corr, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Unicycler
def unicycler(assembly, illumina_corr_1, illumina_corr_2, nanopore_corr, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("30-Unicycler") == False: os.mkdir("30-Unicycler")
	if os.path.isdir("30-Unicycler/" + folder) == False: os.mkdir("30-Unicycler/" + folder)

	# GWF
	inputs = ["{}".format(assembly), "{}".format(illumina_corr_1), "{}".format(illumina_corr_2), "{}".format(nanopore_corr)]
	outputs = ["{}/assembly.fasta".format(out_dir)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '12:00:00'}

	spec='''
	unicycler -1 {illumina_corr_1} -2 {illumina_corr_2} --existing_long_read_assembly {assembly} -l {nanopore_corr} --threads {threads} --keep 2 --verbosity 2 -o {out_dir}
	'''.format(assembly=assembly, illumina_corr_1=illumina_corr_1, illumina_corr_2=illumina_corr_2, nanopore_corr=nanopore_corr, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Quast
def quast(assembly, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("31-Quast") == False: os.mkdir("31-Quast")
	if os.path.isdir("31-Quast/" + folder) == False: os.mkdir("31-Quast/" + folder)

	# GWF
	inputs = ["{}".format(assembly)]
	outputs = ["{}/report.tsv".format(out_dir)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '12:00:00'}

	spec='''
	quast -o {out_dir} -t {threads} {assembly}
	'''.format(assembly=assembly, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Alignment
def align_illumina(assembly, illumina_corr_1, illumina_corr_2, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("40-Align") == False: os.mkdir("40-Align")
	if os.path.isdir("40-Align/" + folder) == False: os.mkdir("40-Align/" + folder)

	# GWF
	inputs = ["{}".format(illumina_corr_1), "{}".format(illumina_corr_2), "{}".format(assembly)]
	outputs = ["{}/Illumina.sort.bam".format(out_dir), "{}/Illumina.sort.bai".format(out_dir)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '4:00:00'}

	spec='''
	mkdir -p {in_dir}/index_contigs
	bowtie2-build --threads {threads} {assembly} {in_dir}/index_contigs/index
	bowtie2 -x {in_dir}/index_contigs/index -1 {illumina_corr_1} -2 {illumina_corr_2} -S {out_dir}/Illumina.sam --threads {threads} 2> {out_dir}/bowtie2.log

	samtools view -bS {out_dir}/Illumina.sam -@ {threads} > {out_dir}/Illumina.bam
	samtools sort -o {out_dir}/Illumina.sort.bam -O bam {out_dir}/Illumina.bam -@ {threads}
	samtools index -b {out_dir}/Illumina.sort.bam {out_dir}/Illumina.sort.bai -@ {threads}
	rm {out_dir}/Illumina.sam {out_dir}/Illumina.bam
	'''.format(assembly=assembly, illumina_corr_1=illumina_corr_1, illumina_corr_2=illumina_corr_2, in_dir='/'.join(assembly.split("/")[:-1]), out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def align_nanopore(assembly, nanopore_corr, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("40-Align") == False: os.mkdir("40-Align")
	if os.path.isdir("40-Align/" + folder) == False: os.mkdir("40-Align/" + folder)

	# GWF
	inputs = ["{}".format(nanopore_corr), "{}".format(assembly)]
	outputs = ["{}/Nanopore.sort.bam".format(out_dir), "{}/Nanopore.sort.bai".format(out_dir)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '4:00:00'}

	spec='''
	minimap2 -a -o {out_dir}/Nanopore.sam -t {threads} -x map-ont {assembly} {nanopore_corr}

	samtools view -bS {out_dir}/Nanopore.sam -@ {threads} > {out_dir}/Nanopore.bam
	samtools sort -o {out_dir}/Nanopore.sort.bam -O bam {out_dir}/Nanopore.bam -@ {threads}
	samtools index -b {out_dir}/Nanopore.sort.bam {out_dir}/Nanopore.sort.bai -@ {threads}
	rm {out_dir}/Nanopore.sam {out_dir}/Nanopore.bam
	'''.format(assembly=assembly, nanopore_corr=nanopore_corr, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Coverage
def coverage(in_dir, out_dir, memory, folder):
	# Folder structure
	if os.path.isdir("50-Coverage") == False: os.mkdir("50-Coverage")
	if os.path.isdir("50-Coverage/" + folder) == False: os.mkdir("50-Coverage/" + folder)

	# GWF
	inputs = ["{}/Illumina.sort.bam".format(in_dir), "{}/Nanopore.sort.bam".format(in_dir)]
	outputs = ["{}/Illumina.cov".format(out_dir), "{}/Nanopore.cov".format(out_dir)]
	options = {'cores': 1,'memory': '{}g'.format(memory), 'queue': 'short', 'walltime': '2:00:00'}

	spec='''
	bedtools genomecov -d -ibam {in_dir}/Illumina.sort.bam > {out_dir}/Illumina.cov
	bedtools genomecov -d -ibam {in_dir}/Nanopore.sort.bam > {out_dir}/Nanopore.cov
	'''.format(in_dir=in_dir, out_dir=out_dir)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Plot Coverage
def plot_coverage(in_dir, out_dir, memory, folder):
	# Folder structure
	if os.path.isdir("60-Plots") == False: os.mkdir("60-Plots")
	if os.path.isdir("60-Plots/" + folder) == False: os.mkdir("60-Plots/" + folder)

	# Number of contigs
	num_contigs = len([1 for line in open("30-Unicycler/{}/flye/assembly.fasta".format(folder)) if line.startswith(">")])

	# GWF
	inputs = ["{}/Illumina.cov".format(in_dir), "{}/Nanopore.cov".format(in_dir)]
	outputs = ["{}/{}.pdf".format(out_dir, num) for num in range(1,num_contigs)]
	options = {'cores': 1,'memory': '{}g'.format(memory), 'queue': 'short', 'walltime': '2:00:00'}

	spec='''
	Rscript coverage.R -i {in_dir} -o {out_dir}
	'''.format(in_dir=in_dir, out_dir=out_dir)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Annotation
def prokka_annot(assembly, out_dir, threads, memory, date, folder):
	# Folder structure
	if os.path.isdir("70-Prokka") == False: os.mkdir("70-Prokka")
	if os.path.isdir("70-Prokka/" + folder) == False: os.mkdir("70-Prokka/" + folder)

	# GWF
	inputs = ["{}".format(assembly)]
	outputs = ["{}/PROKKA_{}.{}".format(out_dir, date, ext) for ext in ["err","faa","ffn","fna","fsa","gbk","gff","log","sqn","tbl","tsv","txt"]]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'short', 'walltime': '2:00:00'}

	spec='''
	prokka --force --cpus {threads} --outdir {out_dir} {assembly}
	'''.format(assembly=assembly, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Execution
#----------

## Strains/Isolates (nanopore illumina_1 illumina_2 folder_name)
# Example: 00-Data/TA1/SRR11719982_pass.fastq.gz 00-Data/TA1/SRR3927460_pass_1.fastq.gz 00-Data/TA1/SRR3927460_pass_2.fastq.gz TA1
# Open file:
file = "strains.tsv"

try:
	f = open(file, 'r')
except:
	print("Error: summary file with genomes not found.")
	sys.exit()

# Loop
for row in f:

	# Divide line
	row = row.strip().split("\t")

	# Variables
	folder = row[3]
	illumina_corr_1 = row[1].split("/")[-1].replace(".gz","") + ".00.0_0.cor.fastq.gz"
	illumina_corr_2 = row[2].split("/")[-1].replace(".gz","") + ".00.0_0.cor.fastq.gz"
	date=datetime.today().strftime('%m%d%Y')

	# 01 QC
	gwf.target_from_template('{}_01_qc_illumina'.format(folder), qc_illumina(illumina_1=row[1], illumina_2=row[2], out_dir="01-QC/{}/Illumina".format(folder), threads=4, memory=8, folder=folder))
	gwf.target_from_template('{}_01_qc_nanopore'.format(folder), qc_nanopore(nanopore=row[0], out_dir="01-QC/{}/Nanopore".format(folder), threads=4, memory=8, folder=folder))

	# 10 Correction
	gwf.target_from_template('{}_10_correct_illumina'.format(folder), correct_illumina(illumina_1=row[1], illumina_2=row[2], illumina_corr_1=illumina_corr_1, illumina_corr_2=illumina_corr_2, out_dir="10-Correction/{}/Illumina".format(folder), threads=4, memory=24, folder=folder))
	gwf.target_from_template('{}_10_correct_nanopore'.format(folder), correct_nanopore(nanopore=row[0], illumina_corr="10-Correction/{}/Illumina/corrected/illumina.corrected.fastq.gz".format(folder), out_dir="10-Correction/{}/Nanopore".format(folder), threads=4, memory=24, folder=folder))

	# 20 Assembly
	gwf.target_from_template('{}_20_assembly_flye'.format(folder), flye_assembly(nanopore_corr="10-Correction/{}/Nanopore/nanopore.corrected.fasta".format(folder), out_dir="20-Assembly/{}/flye".format(folder), threads=8, memory=64, folder=folder))

	# 30 Unicycler
	gwf.target_from_template("{}_30_unicycler".format(folder), unicycler(assembly="20-Assembly/{}/flye/assembly.fasta".format(folder), illumina_corr_1="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_1), illumina_corr_2="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_2), nanopore_corr="10-Correction/{}/Nanopore/nanopore.corrected.fasta".format(folder), out_dir="30-Unicycler/{}/flye".format(folder), threads=8, memory=64, folder=folder))

	# 31 Quast
	gwf.target_from_template("{}_31_quast".format(folder), quast(assembly="30-Unicycler/{}/flye/assembly.fasta".format(folder), out_dir="31-Quast/{}/flye".format(folder), threads=8, memory=64, folder=folder))

	# 40 Alignment
	gwf.target_from_template("{}_40_align_illumina".format(folder), align_illumina(assembly="30-Unicycler/{}/flye/assembly.fasta".format(folder), illumina_corr_1="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_1), illumina_corr_2="10-Correction/{}/Illumina/corrected/{}".format(folder, illumina_corr_2), out_dir="40-Align/{}".format(folder), threads=4, memory=16, folder=folder))
	gwf.target_from_template("{}_40_align_nanopore".format(folder), align_nanopore(assembly="30-Unicycler/{}/flye/assembly.fasta".format(folder), nanopore_corr="10-Correction/{}/Nanopore/nanopore.corrected.fasta".format(folder), out_dir="40-Align/{}".format(folder), threads=4, memory=16, folder=folder))

	# 50 Coverage
	gwf.target_from_template("{}_50_coverage".format(folder), coverage(in_dir="40-Align/{}".format(folder), out_dir="50-Coverage/{}".format(folder), memory=8, folder=folder))

	# 60 Coverage Plot
	gwf.target_from_template("{}_60_plot_coverage".format(folder), plot_coverage(in_dir="50-Coverage/{}".format(folder), out_dir="60-Plots/{}".format(folder), memory=8, folder=folder))

	# 70 Annotation
	gwf.target_from_template("{}_70_annotation".format(folder), prokka_annot(assembly="30-Unicycler/{}/flye/assembly.fasta".format(folder), out_dir="70-Prokka/{}".format(folder), threads=4, memory=16, date=date, folder=folder))