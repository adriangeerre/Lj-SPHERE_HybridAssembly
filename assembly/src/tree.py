import os
from gwf import AnonymousTarget

## GWF function
# Mashtree
def mashtree(in_dir, breps, threads, memory):
	# GWF
	inputs = ["{}".format(in_dir)]
	outputs = ["LjSC.tree", "LjSC.distmat"]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'short', 'walltime': '12:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

    # Split assembly in contigs
    conda activate MashTree
	mashtree_bootstrap.pl --reps {breps} --numcpus {threads} --outmatrix LjSC.distmat {in_dir}/*.fasta > LjSC.tree

	'''.format(in_dir=in_dir, breps=breps, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)