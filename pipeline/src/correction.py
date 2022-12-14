## Imports
import os
import subprocess

# SPAdes
def correct_illumina(illumina_1, illumina_2, illumina_corr_1, illumina_corr_2, out_dir, threads, memory, conda_path, logfile):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	cmd='''/bin/bash -c "
	# Source conda to work with environments
	source {conda_path}/etc/profile.d/conda.sh
	
	# SPAdes
	conda activate SPAdes
	spades.py -1 {illumina_1} -2 {illumina_2} -o {out_dir} --only-error-correction -t {threads} -m {memory}
	cat {out_dir}/corrected/{illumina_corr_1} {out_dir}/corrected/{illumina_corr_2} > {out_dir}/corrected/illumina.corrected.fastq.gz"
	'''.format(illumina_1=illumina_1, illumina_2=illumina_2, illumina_corr_1=illumina_corr_1, illumina_corr_2=illumina_corr_2, out_dir=out_dir, threads=threads, memory=memory, conda_path=conda_path)

	# Exec and log
	f = open(logfile, "a")
	subprocess.check_call(cmd, shell=True, stdout=f, stderr=f)
	f.close()

# LoRDEC
def correct_nanopore(nanopore, illumina_corr, out_dir, threads, conda_path, logfile):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	cmd='''/bin/bash -c "
	# Source conda to work with environments
	source {conda_path}/etc/profile.d/conda.sh

	# LoRDEC
	conda activate ha-flye
	lordec-correct -i {nanopore} -2 {illumina_corr} -k 19 -s 4 -T {threads} -p -o {out_dir}/nanopore.kmer19.fasta
	lordec-correct -i {out_dir}/nanopore.kmer19.fasta -2 {illumina_corr} -k 31 -s 3 -T {threads} -p -o {out_dir}/nanopore.kmer31.fasta
	lordec-correct -i {out_dir}/nanopore.kmer31.fasta -2 {illumina_corr} -k 41 -s 3 -T {threads} -p -o {out_dir}/nanopore.corrected.fasta"
	'''.format(nanopore=nanopore, illumina_corr=illumina_corr, out_dir=out_dir, threads=threads, conda_path=conda_path)

	# Exec and log
	f = open(logfile, "a")
	subprocess.check_call(cmd, shell=True, stdout=f, stderr=f)
	f.close()