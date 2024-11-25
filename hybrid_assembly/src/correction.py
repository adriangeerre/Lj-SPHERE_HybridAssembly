## Imports
import os
from gwf import AnonymousTarget

## GWF function
# SPAdes
def correct_illumina(illumina_1, illumina_2, illumina_corr_1, illumina_corr_2, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("10-Correction/" + folder) == False: os.makedirs("10-Correction/" + folder)

	# GWF
	inputs = ["{}".format(illumina_1), "{}".format(illumina_2)]
	outputs = ["{}/corrected/illumina.corrected.fastq.gz".format(out_dir), "{}/corrected/{}".format(out_dir, illumina_corr_1), "{}/corrected/{}".format(out_dir, illumina_corr_2)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '24:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# SPAdes
	conda activate SPAdes
	spades.py -1 {illumina_1} -2 {illumina_2} -o {out_dir} --only-error-correction -t {threads} -m {memory}
	cat {out_dir}/corrected/{illumina_corr_1} {out_dir}/corrected/{illumina_corr_2} > {out_dir}/corrected/illumina.corrected.fastq.gz
	'''.format(illumina_1=illumina_1, illumina_2=illumina_2, illumina_corr_1=illumina_corr_1, illumina_corr_2=illumina_corr_2, out_dir=out_dir, threads=threads, memory=memory)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)	

# LoRDEC
def correct_nanopore(nanopore, illumina_corr, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("10-Correction/" + folder + "/Nanopore") == False: os.makedirs("10-Correction/" + folder + "/Nanopore")

	# GWF
	inputs = ["{}".format(nanopore), "{}".format(illumina_corr)]
	outputs = ["{}/nanopore.corrected.fasta".format(out_dir)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '24:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# LoRDEC
	conda activate ha-flye
	lordec-correct -i {nanopore} -2 {illumina_corr} -k 19 -s 4 -T {threads} -p -o {out_dir}/nanopore.kmer19.fasta
	lordec-correct -i {out_dir}/nanopore.kmer19.fasta -2 {illumina_corr} -k 31 -s 3 -T {threads} -p -o {out_dir}/nanopore.kmer31.fasta
	lordec-correct -i {out_dir}/nanopore.kmer31.fasta -2 {illumina_corr} -k 41 -s 3 -T {threads} -p -o {out_dir}/nanopore.corrected.fasta
	'''.format(nanopore=nanopore, illumina_corr=illumina_corr, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)