## Imports
import glob, os, sys
from subprocess import check_call

# Fastqc
def qc_illumina(illumina_1, illumina_2, out_dir, threads):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# Execution
	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# FastQC
	conda activate ha-flye
	fastqc -t {threads} -o {out_dir} {illumina_1} {illumina_2}
	'''.format(illumina_1=illumina_1, illumina_2=illumina_2, out_dir=out_dir, threads=threads)

	check_call(spec, shell=True)

# NanoPlot
def qc_nanopore(nanopore, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# NanoPlot
	conda activate NanoQC
	NanoPlot -o {out_dir} -p {prefix} --info_in_report --N50 --title {title} --fastq {nanopore} --threads {threads}
	'''.format(nanopore=nanopore, out_dir=out_dir, prefix="{}-".format(nanopore.split("/")[-1].replace(".fastq.gz","")), title=nanopore.split("/")[-1], threads=threads)

	check_call(spec, shell=True)