## Imports
import os
from gwf import AnonymousTarget

## GWF function
# Flye
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
