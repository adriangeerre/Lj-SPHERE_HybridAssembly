## Imports
import os
import subprocess

# Flye
def flye_assembly(nanopore_corr, out_dir, threads, conda_path):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	cmd='''/bin/bash -c "
	# Source conda to work with environments
	source {conda_path}/etc/profile.d/conda.sh

	# Flye
	conda activate ha-flye
	flye --nano-corr {nanopore_corr} --plasmids --out-dir {out_dir} --threads {threads}"
	'''.format(nanopore_corr=nanopore_corr, out_dir=out_dir, threads=threads, conda_path=conda_path)

	subprocess.check_call(cmd, shell=True)

# Unicycler
def unicycler(assembly, illumina_corr_1, illumina_corr_2, nanopore_corr, out_dir, threads, conda_path):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	cmd='''/bin/bash -c "
	# Source conda to work with environments
	source {conda_path}/etc/profile.d/conda.sh

	# Unicycler
	conda activate Unicycler
	unicycler -1 {illumina_corr_1} -2 {illumina_corr_2} --existing_long_read_assembly {assembly} -l {nanopore_corr} --threads {threads} --keep 2 --verbosity 2 -o {out_dir}"
	'''.format(assembly=assembly, illumina_corr_1=illumina_corr_1, illumina_corr_2=illumina_corr_2, nanopore_corr=nanopore_corr, out_dir=out_dir, threads=threads, conda_path=conda_path)

	subprocess.check_call(cmd, shell=True)

