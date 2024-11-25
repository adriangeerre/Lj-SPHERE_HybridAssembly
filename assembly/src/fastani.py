## GWF function
# Fastani
def fastani(hybrid_assembly, illumina_genomes, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("90-FastANI") == False: os.makedirs("90-FastANI")
	if os.path.isdir("90-FastANI/" + folder) == False: os.makedirs("90-FastANI/" + folder)

	# GWF
	inputs = ["{}".format(hybrid_assembly), "{}".format(illumina_genomes)]
	outputs = ["{}/{}.fastani.tsv".format(out_dir, folder)]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '00:30:00'}
	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# FastANI
	conda activate FastANI
	fastANI -q {hybrid_assembly} --rl {illumina_genomes} -o {out_dir}/{folder}.fastani.tsv -t {threads}
	'''.format(hybrid_assembly=hybrid_assembly, illumina_genomes=illumina_genomes, out_dir=out_dir, folder=folder, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)