## Imports
import os
from gwf import AnonymousTarget

## GWF function
# Fastqc
def qc_illumina(illumina_1, illumina_2, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("01-QC/" + folder + "/Illumina") == False: os.makedirs("01-QC/" + folder + "/Illumina")

	# Names
	fqc_illumina_1 = illumina_1.split("/")[-1].split(".")[0] + "_fastqc"
	fqc_illumina_2 = illumina_2.split("/")[-1].split(".")[0] + "_fastqc"

	# GWF
	inputs = ["{}".format(illumina_1), "{}".format(illumina_2)]
	outputs = ["{}/{}.zip".format(out_dir, fqc_illumina_1), "{}/{}.html".format(out_dir, fqc_illumina_1), "{}/{}.zip".format(out_dir, fqc_illumina_2), "{}/{}.html".format(out_dir, fqc_illumina_2)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '2:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# FastQC
	conda activate ha-flye
	fastqc -t {threads} -o {out_dir} {illumina_1} {illumina_2}
	'''.format(illumina_1=illumina_1, illumina_2=illumina_2, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# NanoPlot
def qc_nanopore(nanopore, out_dir, threads, memory, folder):
	# Folder structure
	if os.path.isdir("01-QC/" + folder) == False: os.makedirs("01-QC/" + folder)

	# GWF
	inputs = ["{}".format(nanopore)]
	outputs = ["{}/{}-NanoPlot-report.html".format(out_dir, nanopore.split("/")[-1].replace(".fastq.gz",""), )]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '4:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# NanoPlot
	conda activate NanoQC
	NanoPlot -o {out_dir} -p {prefix} --info_in_report --N50 --title {title} --fastq {nanopore} --threads {threads}
	'''.format(nanopore=nanopore, out_dir=out_dir, prefix="{}-".format(nanopore.split("/")[-1].replace(".fastq.gz","")), title=nanopore.split("/")[-1], threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)