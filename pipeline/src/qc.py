## Imports
import os
import subprocess

# Fastqc
def qc_illumina(illumina_1, illumina_2, out_dir, threads, conda_path, logfile):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# Execution
	cmd='''/bin/bash -c "
	# Source conda to work with environments	
	source {conda_path}/etc/profile.d/conda.sh
	
	# FastQC
	conda activate ha-flye
	fastqc -t {threads} -o {out_dir} {illumina_1} {illumina_2}"
	'''.format(illumina_1=illumina_1, illumina_2=illumina_2, out_dir=out_dir, threads=threads, conda_path=conda_path)

	# Exec and log
	f = open(logfile, "w")
	subprocess.check_call(cmd, shell=True, stdout=f, stderr=f)
	f.close()

# NanoPlot
def qc_nanopore(nanopore, out_dir, threads, conda_path, logfile):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	cmd='''/bin/bash -c "
	# Source conda to work with environments
	source {conda_path}/etc/profile.d/conda.sh

	# NanoPlot
	conda activate NanoQC
	NanoPlot -o {out_dir} -p {prefix} --info_in_report --N50 --title {title} --fastq {nanopore} --threads {threads}"
	'''.format(nanopore=nanopore, out_dir=out_dir, prefix=".".join(nanopore.split("/")[-1].split(".")[:-2]) + "_", title=nanopore.split("/")[-1], threads=threads, conda_path=conda_path)

	# Exec and log
	f = open(logfile, "w")
	subprocess.check_call(cmd, shell=True, stdout=f, stderr=f)
	f.close()