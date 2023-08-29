## Imports
import os
import matplotlib.pyplot as plt
from gwf import AnonymousTarget

# Functions
def gc_by_depth(assembly):
	# GC
	values = {}
	# Read into dictionary
	d = SeqIO.to_dict(SeqIO.parse(assembly, 'fasta'))
	# Compute gc
	for i, j in d.items():
		# Depth
		depth = [k[6:][:-1] for k in j.description.split(" ") if "depth=" in k][0]
		# GC
		g = j.seq.count("G")
		c = j.seq.count("C")
		gc = (g + c)/len(j.seq) * 100
		l = len(j.seq)
		values[i] = [gc, float(depth), l]
	# Return
	return(values)

def plot_gc_vs_cov(assembly, outpath, name):
	# Get values
	d = gc_by_depth(assembly)
	# Divide information
	names = list(d.keys())
	gc = [i[0] for i in d.values()]
	cov = [i[1] for i in d.values()]
	lens = [i[2] / 1000 for i in d.values()]
	# Plot
	plt.figure(figsize=(8, 6))
	plt.scatter(gc, cov, s=lens, c='blue', marker='o', edgecolors='black', alpha=0.3)
	for i, name in enumerate(names):
		plt.annotate(name, (gc[i], cov[i]), textcoords="offset points", xytext=(0,10), ha='center')
	plt.title('GC Content vs Coverage')
	plt.xlabel('GC Content')
	plt.ylabel('Coverage')
	plt.savefig(f"{outpath}/{name}.GCvsCOV.png")

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
def plot_coverage(in_dir, out_dir, memory, folder):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# Number of contigs
	if os.path.exists("30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder)):
		num_contigs = len([1 for line in open("30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder)) if line.startswith(">")])
	else:
		num_contigs = 1

	# GWF
	inputs = ["{}/Illumina.cov".format(in_dir), "{}/Nanopore.cov".format(in_dir)]
	outputs = ["{}/{}.pdf".format(out_dir, num) for num in range(1,num_contigs+1) if os.path.exists("30-HybridAssembly/{}/unicycler/assembly.fasta".format(folder))]
	options = {'cores': 1,'memory': '{}g'.format(memory), 'queue': 'short', 'walltime': '2:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# R
	conda activate Renv
	Rscript src/coverage.R -i {in_dir} -o {out_dir}
	'''.format(in_dir=in_dir, out_dir=out_dir)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)