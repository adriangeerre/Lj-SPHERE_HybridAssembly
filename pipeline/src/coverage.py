## Imports
import os
from Bio import SeqIO
from gwf import AnonymousTarget

# Functions
def gc_by_depth(assembly):
	# GC
	values = {}
	# Read into dictionary
	d = SeqIO.to_dict(SeqIO.parse(assembly, 'fasta'))
	# Compute gc
	for i, j in d.items():
		# Circularity
		circular = max([1 if "circular=" in k else 0 for k in j.description.split(" ")])
		# GC
		g = j.seq.count("G")
		c = j.seq.count("C")
		gc = (g + c)/len(j.seq) * 100
		gsize = len(j.seq)
		values[i] = [gc, gsize, circular]
	# Return
	return(values)

## GWF function
# Alignment Illumina
def align_illumina(assembly, illumina_corr_1, illumina_corr_2, out_dir, threads, memory):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# GWF
	inputs = ["{}".format(illumina_corr_1), "{}".format(illumina_corr_2), "{}".format(assembly)]
	outputs = ["{}/Illumina.sort.bam".format(out_dir), "{}/Illumina.sort.bai".format(out_dir)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '4:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Bowtie2 + Samtools
	conda activate ha-flye
	mkdir -p {in_dir}/index_contigs
	bowtie2-build --threads {threads} {assembly} {in_dir}/index_contigs/index
	bowtie2 -x {in_dir}/index_contigs/index -1 {illumina_corr_1} -2 {illumina_corr_2} -S {out_dir}/Illumina.sam --threads {threads} 2> {out_dir}/bowtie2.log

	samtools view -bS {out_dir}/Illumina.sam -@ {threads} > {out_dir}/Illumina.bam
	samtools sort -o {out_dir}/Illumina.sort.bam -O bam {out_dir}/Illumina.bam -@ {threads}
	samtools index -b {out_dir}/Illumina.sort.bam {out_dir}/Illumina.sort.bai -@ {threads}
	rm {out_dir}/Illumina.sam {out_dir}/Illumina.bam
	'''.format(assembly=assembly, illumina_corr_1=illumina_corr_1, illumina_corr_2=illumina_corr_2, in_dir='/'.join(assembly.split("/")[:-1]), out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Alignment Nanopore
def align_nanopore(assembly, nanopore_corr, out_dir, threads, memory):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# GWF
	inputs = ["{}".format(nanopore_corr), "{}".format(assembly)]
	outputs = ["{}/Nanopore.sort.bam".format(out_dir), "{}/Nanopore.sort.bai".format(out_dir)]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '4:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Minimap2 + Samtools
	conda activate ha-flye
	minimap2 -a -o {out_dir}/Nanopore.sam -t {threads} -x map-ont {assembly} {nanopore_corr}

	samtools view -bS {out_dir}/Nanopore.sam -@ {threads} > {out_dir}/Nanopore.bam
	samtools sort -o {out_dir}/Nanopore.sort.bam -O bam {out_dir}/Nanopore.bam -@ {threads}
	samtools index -b {out_dir}/Nanopore.sort.bam {out_dir}/Nanopore.sort.bai -@ {threads}
	rm {out_dir}/Nanopore.sam {out_dir}/Nanopore.bam
	'''.format(assembly=assembly, nanopore_corr=nanopore_corr, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Coverage
def coverage(in_dir, out_dir, memory):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# GWF
	inputs = ["{}/Illumina.sort.bam".format(in_dir), "{}/Nanopore.sort.bam".format(in_dir)]
	outputs = ["{}/Illumina.cov".format(out_dir), "{}/Nanopore.cov".format(out_dir)]
	options = {'cores': 1,'memory': '{}g'.format(memory), 'queue': 'short', 'walltime': '2:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Bedtools
	conda activate Bedtools
	bedtools genomecov -d -ibam {in_dir}/Illumina.sort.bam > {out_dir}/Illumina.cov
	bedtools genomecov -d -ibam {in_dir}/Nanopore.sort.bam > {out_dir}/Nanopore.cov
	'''.format(in_dir=in_dir, out_dir=out_dir)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Plot Coverage
def plot_coverage(assembly, in_dir, out_dir, memory, folder):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# Contig data
	contig_data = "{}/{}.contig-data.tsv".format(in_dir, folder)
	if os.path.exists(assembly):
		if not os.path.exists(contig_data):
			d = gc_by_depth(assembly)
			w = open(contig_data, "w")
			for i, j in d.items():
				w.write(f"{i}\t{round(j[0],4)}\t{j[1]}\t{j[2]}\n")
			w.close()

	# Define outputs
	outputs = []
	if os.path.exists(contig_data):
		r = open(contig_data, "r")
		contigs = [i.strip().split("\t")[0] for i in r]
		outputs = ["{}/{}.pdf".format(out_dir, contig) for contig in contigs]

	# GWF
	inputs = ["{}/Illumina.cov".format(in_dir), "{}/Nanopore.cov".format(in_dir)]
	options = {'cores': 1,'memory': '{}g'.format(memory), 'queue': 'short', 'walltime': '2:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# R
	conda activate Renv
	Rscript src/coverage.R -i {in_dir} -o {out_dir}
	'''.format(in_dir=in_dir, out_dir=out_dir)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Plot GC vs Coverage
def plot_gc_vs_coverage(assembly, in_dir, out_dir, memory, folder):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# Contig data
	contig_data = "{}/{}.contig-data.tsv".format(in_dir, folder)
	if os.path.exists(assembly):
		if not os.path.exists(contig_data):
			d = gc_by_depth(assembly)
			w = open(contig_data, "w")
			for i, j in d.items():
				w.write(f"{i}\t{round(j[0],4)}\t{j[1]}\t{j[2]}\n")
			w.close()

	# GWF
	inputs = [assembly, "{}/Illumina.cov".format(in_dir), "{}/Nanopore.cov".format(in_dir)]
	outputs = ["{}/{}.Illumina.mean.GCvsCOV.png".format(out_dir, folder), "{}/{}.Illumina.median.GCvsCOV.png".format(out_dir, folder), "{}/{}.Nanopore.mean.GCvsCOV.png".format(out_dir, folder), "{}/{}.Nanopore.median.GCvsCOV.png".format(out_dir, folder)]
	options = {'cores': 1,'memory': '{}g'.format(memory), 'queue': 'short', 'walltime': '2:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# R
	conda activate Renv
	Rscript src/gcvscov.R -a {assembly} -p {folder} -i {in_dir} -o {out_dir}
	'''.format(assembly=assembly, folder=folder, in_dir=in_dir, out_dir=out_dir)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)