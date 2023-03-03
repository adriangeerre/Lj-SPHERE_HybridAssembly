## GWF function
# Mashtree
def mashtree(assembly, out_dir, in_dir, threads, memory):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# GWF
	inputs = ["{}".format(assembly)]
	outputs = ["{}/LjSC.tree".format(out_dir), "{}/LjSC.distmat".format(out_dir)]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'short', 'walltime': '8:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

    # Split assembly in contigs
    conda activate MashTree
    mashtree --numcpus {threads} --outmatrix {out_dir}/LjSC.distmat --outtree {out_dir}/LjSC.tree --mindepth 0 {in_dir}/*.fasta

    # Plot tree

	'''.format(assembly=assembly, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)