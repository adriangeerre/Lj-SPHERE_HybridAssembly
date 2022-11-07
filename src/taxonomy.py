## Imports
import os
from gwf import AnonymousTarget

## GWF function
# Obtain 16S
def taxonomy_16S(assembly, target, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# GWF
	inputs = ["{}".format(assembly)]
	outputs = ["{}/recovered_16S.tsv".format(out_dir)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '24:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Blast
	conda activate Blast
    blastn -subject {assembly} -query {target} -outfmt "6 std sseq" | cut -f 1,13 | tr "\t" "\n" | tr -d "-" > {out_dir}/recovered_16S.fasta

    # Remove duplicates
    conda activate SeqKit
    seqkit rmdup -s {out_dir}/recovered_16S.fasta > {out_dir}/recovered_16S.de-rep.tsv

    # Blast
    conda activate Blast
    blastn -db nt -query {out_dir}/recovered_16S.de-rep.tsv -online -outfmt "6 std" > taxonomy_16S.tsv
	'''.format(assembly=assembly, target=target, out_dir=out_dir, threads=threads, memory=memory)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)	