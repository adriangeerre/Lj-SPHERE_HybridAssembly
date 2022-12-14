## Imports
import os
import yaml
import subprocess

## Auxiliary function
# PGAP yaml
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

# Annotation
def annotation(input_yaml, out_dir, threads, memory, conda_path, logfile):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	cmd='''
	# Cache and tmp folders
	export SINGULARITY_CACHEDIR=/scratch/$SLURM_JOBID
	export SINGULARITY_TMPDIR=/scratch/$SLURM_JOBID

	# PGAP (already in path)
	python /home/agomez/programas/PGAP/pgap.py -d -n --no-internet --ignore-all-errors --docker singularity -o {out_dir}/annotation --memory {memory} --container-path ~/programas/SingularityImages/pgap_2022-08-11.build6275.sif {input_yaml}
	'''.format(input_yaml=input_yaml, out_dir=out_dir, memory=memory, conda_path=conda_path)

	# Exec and log
	f = open(logfile, "a")
	subprocess.check_call(cmd, shell=True, stdout=f, stderr=f)
	f.close()

