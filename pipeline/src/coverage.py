## Imports
import os
import subprocess

# Alignment Illumina
def align_illumina(assembly, illumina_corr_1, illumina_corr_2, out_dir, threads, conda_path, logfile):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	cmd='''/bin/bash -c "
	# Source conda to work with environments
	source {conda_path}/etc/profile.d/conda.sh

	# Bowtie2 + Samtools
	conda activate ha-flye
	mkdir -p {in_dir}/index_contigs
	bowtie2-build --threads {threads} {assembly} {in_dir}/index_contigs/index
	bowtie2 -x {in_dir}/index_contigs/index -1 {illumina_corr_1} -2 {illumina_corr_2} -S {out_dir}/Illumina.sam --threads {threads} 2> {out_dir}/bowtie2.log

	samtools view -bS {out_dir}/Illumina.sam -@ {threads} > {out_dir}/Illumina.bam
	samtools sort -o {out_dir}/Illumina.sort.bam -O bam {out_dir}/Illumina.bam -@ {threads}
	samtools index -b {out_dir}/Illumina.sort.bam {out_dir}/Illumina.sort.bai -@ {threads}
	rm {out_dir}/Illumina.sam {out_dir}/Illumina.bam"
	'''.format(assembly=assembly, illumina_corr_1=illumina_corr_1, illumina_corr_2=illumina_corr_2, in_dir='/'.join(assembly.split("/")[:-1]), out_dir=out_dir, threads=threads, conda_path=conda_path)

	# Exec and log
	f = open(logfile, "a")
	subprocess.check_call(cmd, shell=True, stdout=f, stderr=f)
	f.close()

# Alignment Nanopore
def align_nanopore(assembly, nanopore_corr, out_dir, threads, conda_path, logfile):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	cmd='''/bin/bash -c "
	# Source conda to work with environments
	source {conda_path}/etc/profile.d/conda.sh

	# Minimap2 + Samtools
	conda activate ha-flye
	minimap2 -a -o {out_dir}/Nanopore.sam -t {threads} -x map-ont {assembly} {nanopore_corr}

	samtools view -bS {out_dir}/Nanopore.sam -@ {threads} > {out_dir}/Nanopore.bam
	samtools sort -o {out_dir}/Nanopore.sort.bam -O bam {out_dir}/Nanopore.bam -@ {threads}
	samtools index -b {out_dir}/Nanopore.sort.bam {out_dir}/Nanopore.sort.bai -@ {threads}
	rm {out_dir}/Nanopore.sam {out_dir}/Nanopore.bam"
	'''.format(assembly=assembly, nanopore_corr=nanopore_corr, out_dir=out_dir, threads=threads, conda_path=conda_path)

	# Exec and log
	f = open(logfile, "a")
	subprocess.check_call(cmd, shell=True, stdout=f, stderr=f)
	f.close()

# Coverage
def coverage(in_dir, out_dir, conda_path, logfile):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	cmd='''/bin/bash -c "
	# Source conda to work with environments
	source {conda_path}/etc/profile.d/conda.sh

	# Bedtools
	conda activate Bedtools
	bedtools genomecov -d -ibam {in_dir}/Illumina.sort.bam > {out_dir}/Illumina.cov
	bedtools genomecov -d -ibam {in_dir}/Nanopore.sort.bam > {out_dir}/Nanopore.cov"
	'''.format(in_dir=in_dir, out_dir=out_dir, conda_path=conda_path)

	# Exec and log
	f = open(logfile, "a")
	subprocess.check_call(cmd, shell=True, stdout=f, stderr=f)
	f.close()

# Plot Coverage
def plot_coverage(in_dir, out_dir, conda_path, logfile):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	cmd='''/bin/bash -c "
	# Source conda to work with environments
	source {conda_path}/etc/profile.d/conda.sh

	# R
	conda activate Renv
	Rscript src/coverage.R -i {in_dir} -o {out_dir}"
	'''.format(in_dir=in_dir, out_dir=out_dir, conda_path=conda_path)

	# Exec and log
	f = open(logfile, "a")
	subprocess.check_call(cmd, shell=True, stdout=f, stderr=f)
	f.close()
