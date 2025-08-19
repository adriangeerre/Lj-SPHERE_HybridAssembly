# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By    : Adrián Gómez Repollés
# Email         : adrian.gomez@mbg.au.dk
# Created Date  : 24/01/2024
# version       : '1.0'
# ---------------------------------------------------------------------------
""" Pipeline to identify plasmid contigs in bacterial genomes. """
# ---------------------------------------------------------------------------
# Imports 
# ---------------------------------------------------------------------------

# Internal

# External
import os
import sys
import glob
from Bio import SeqIO
from gwf import Workflow, AnonymousTarget

gwf = Workflow(defaults={"account": "CCRP_Data"})

## Functions

# geNomad
def geNomad(fasta, outdir, name, db, threads, memory):
	# Folder structure
	if not os.path.isdir(outdir): os.makedirs(outdir)

	# GWF
	inputs = [fasta]
	outputs = [f"{outdir}/{name}_summary/{name}_{ext}" for ext in ["plasmid.fna","plasmid_genes.tsv","plasmid_proteins.faa","plasmid_summary.tsv","virus.fna","virus_genes.tsv","virus_proteins.faa","virus_summary.tsv"]]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '02:00:00'}
	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Identify plasmids
	conda activate genomad
	genomad end-to-end --threads {threads} {fasta} {outdir} {db} 
	'''.format(fasta=fasta, outdir=outdir, db=db, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# PLASMe
def PLASMe(fasta, coverage, identity, outdir, db, threads, memory):
	# Folder structure
	if not os.path.isdir(outdir): os.makedirs(outdir)

	# GWF
	inputs = [fasta]
	outputs = [f"{outdir}/results_report.csv"]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '02:00:00'}
	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Link database
	if [ ! -h DB ]; then ln -s {db}; fi

	# Identify plasmids
	conda activate plasme
	python ~/programas/PLASMe/PLASMe.py -c {coverage} -i {identity} -t {threads} --temp {outdir}/temp {fasta} {outdir}/results
	'''.format(fasta=fasta, coverage=coverage, identity=identity, db=db, outdir=outdir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# viralVerify
def viralVerify(fasta, outdir, name, db, threads, memory):
	# Folder structure
	if not os.path.isdir(outdir): os.makedirs(outdir)

	# GWF
	inputs = [fasta]
	outputs = [f"{outdir}/{name}_result_table.csv"]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '02:00:00'}
	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Identify plasmids
	conda activate viralVerify
	python ~/programas/viralVerify/viralVerify-1.1/bin/viralverify -f {fasta} -o {outdir} --hmm {db} -t {threads} -p
	'''.format(fasta=fasta, outdir=outdir, db=db, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Plasmer
def Plasmer(fasta, outdir, name, db, threads, memory):
	# Folder structure
	if not os.path.isdir(outdir): os.makedirs(outdir)

	# GWF
	inputs = [fasta]
	outputs = [f"{outdir}/results/{name}.{ext}" for ext in ["plasmer.predClass.tsv","plasmer.predPlasmids.fa","plasmer.predPlasmids.taxon","plasmer.predProb.tsv","plasmer.shorterM.fasta"]]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '02:00:00'}
	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Identify plasmids
	conda activate Plasmer
	Plasmer -g {fasta} -p {name} -d {db} -t {threads} -l 0 -o {outdir}
	'''.format(fasta=fasta, outdir=outdir, name=name, db=db, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# PlasmidHunter
def PlasmidHunter(fasta, outdir, threads, memory):
	# Folder structure
	if not os.path.isdir(outdir): os.makedirs(outdir)

	# GWF
	inputs = [fasta]
	outputs = [f"{outdir}/predictions.tsv"]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '02:00:00'}
	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Identify plasmids
	conda activate plasmidhunter
	plasmidhunter -i {fasta} -o {outdir} -c {threads}
	'''.format(fasta=fasta, outdir=outdir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Overlap
def overlap(name, folder, threads, memory):
	# GWF
	inputs = [f"01-results/geNomad/{folder}/{folder}_summary/{folder}_plasmid_summary.tsv", f"01-results/PLASMe/{folder}/results_report.csv", f"01-results/viralVerify/{folder}/{folder}_result_table.csv", f"01-results/Plasmer/{folder}/results/{name}.plasmer.predClass.tsv", f"01-results/PlasmidHunter/{folder}/predictions.tsv"]
	outputs = [f"02-Overlap/{folder}.overlap.tsv"]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '00:10:00'}
	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Identify plasmids
	python src/overlap.py -f {folder} -n {name}
	'''.format(folder=folder, name=name)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Divide contigs (FNA)
def contigs2fna(paths, contigs, outdir):
	# Folder structure
	if not os.path.isdir(outdir): os.makedirs(outdir)

	# Load paths into dictionary
	lines = open(paths, "r").readlines()
	paths = {i.strip().split()[0]: i.strip().split()[2] for i in lines}

	# Loop contigs
	with open(contigs, "r") as lines:
		for line in lines:
			if line.split()[0] == "strain":
				continue
			strain = line.strip().split()[0]
			contig = line.strip().split()[1]
			output = f"{outdir}/{strain}_{contig}.fasta"
			if not os.path.exists(output):
				try:
					[SeqIO.write(i, output, "fasta") for i in SeqIO.parse(paths[strain], "fasta") if i.id == contig]
				except:
					print(f"{strain} to do manually!!")

# ANI
def ani(fasta, input_list, outdir, name, threads, memory):
	# Folder structure
	if not os.path.isdir(outdir): os.makedirs(outdir)

	# GWF
	inputs = [fasta]
	genames = [i.split("/")[1].split(".")[0] for i in open(input_list,"r")]
	outputs = [f"{outdir}/{name}_{ext}.tsv" for ext in genames]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '02:00:00'}
	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Compute one-to-many ANI
	conda activate FastANI
	for genome in $(cat {input_list}); do
		ref=$(basename $genome | cut -d . -f 1)
		fastANI -q {fasta} -r $genome -o {outdir}/{name}_$ref.tsv -t {threads}
	done
	'''.format(fasta=fasta, input_list=input_list, name=name, outdir=outdir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Divide contigs (FAA)
def parse_gff(gff, contig):
	# Gene names
	gnames = []
	# Parse gff
	with open(gff, "r") as lines:
		for line in lines:
			line = line.strip().split("\t")
			if line[0] == contig and line[2] == "gene":
				gnames.append([i[3:].replace("gene-","") for i in line[8].split(";") if i[:3] == "ID="][0])
	# Return
	return gnames

def contigs2faa(contigs, outdir):
	# Folder structure
	if not os.path.isdir(outdir): os.makedirs(outdir)

	# Loop contigs
	with open(contigs, "r") as lines:
		for line in lines:
			if line.split()[0] == "strain":
				continue
			# Variables
			strain = line.strip().split()[0]
			contig = line.strip().split()[1]
			# Files
			gff = f"../02-HybridAssembly/00-HA/50-Annotation/{strain}/unicycler/annotation/annot.gff"
			faa = f"../02-HybridAssembly/00-HA/50-Annotation/{strain}/unicycler/annotation/annot.faa"
			output = f"{outdir}/{strain}_{contig}.faa"
			# Gene names
			try:
				gnames = parse_gff(gff, contig)
				if not os.path.exists(output):
					seqs = [i for i in SeqIO.parse(faa, "fasta") if i.id.split("|")[-1] in gnames]
					SeqIO.write(seqs, output, "fasta")
			except:
				print(f"{strain} to do manually!!")

# Orthofinder
def orthofinder(input_list, indir, outdir, threads, memory):
	# GWF
	inputs = [f"{indir}/{i.strip().split('/')[1].replace('fasta','faa')}" for i in open(input_list,"r")]
	outputs = [f"{outdir}/Results_OvPlas/Orthogroups/Orthogroups.GeneCount.tsv"]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '06:00:00'}
	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Orthofinder
	conda activate OrthoFinder
	orthofinder -f {indir} -t {threads} -a 2 -n OvPlas -o {outdir}
	'''.format(indir=indir, outdir=outdir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

## Execution

### input.tsv:
### <GENOME>	<GENOME>.assembly	<PATH>/<GENOME>.assembly.fasta

# Variables
genomad_db = "~/programas/genomad_db"
plasme_db = "~/programas/PLASMe_db/DB"
viral_hmm = "~/programas/viralVerify/nbc_hmms.hmm"
plasmer_db = "~/programas/Plasmer_db"
file = "input.tsv"

# Input samples
try:
	f = open(file, 'r')
except:
	print(f"Error: {file} file not found.")
	sys.exit()

### First Round
for row in f:

	if row[0] == "#":
		continue

	# Divide line
	row = row.strip().split("\t")
	folder = row[0]
	name = row[1]
	fasta = row[2]

	# geNomad
	gwf.target_from_template(f"{folder}_geNomad", geNomad(fasta=fasta, outdir=f"01-results/geNomad/{folder}", name=folder, db=genomad_db, threads=8, memory=24))

	# PLASMe
	gwf.target_from_template(f"{folder}_PLASMe", PLASMe(fasta=fasta, coverage=0.1, identity=0.1, outdir=f"01-results/PLASMe/{folder}", db=plasme_db, threads=8, memory=24))

	# ViralVerify
	gwf.target_from_template(f"{folder}_viralVerify", viralVerify(fasta=fasta, outdir=f"01-results/viralVerify/{folder}", name=folder, db=viral_hmm, threads=8, memory=24))

	# Plasmer
	gwf.target_from_template(f"{folder}_Plasmer", Plasmer(fasta=fasta, outdir=f"01-results/Plasmer/{folder}", name=name, db=plasmer_db, threads=8, memory=24))

	# PlasmidHunter
	gwf.target_from_template(f"{folder}_PlasmidHunter", PlasmidHunter(fasta=fasta, outdir=f"01-results/PlasmidHunter/{folder}", threads=4, memory=12))

	# Overlap
	gwf.target_from_template(f"{folder}_Overlap", overlap(name=name, folder=folder, threads=1, memory=4))

f.close()

### Second round

# Extract contigs as individual fasta
contigs_file = "full_overlap_contigs.tsv" # Produce manually in Jupyter Notebook (ComputeOverlap.ipynb)
if os.path.exists(contigs_file):
	contigs2fna(file, contigs_file, "03-IndFasta")

# ANI
fastas = glob.glob("03-IndFasta/*.fasta")
if not os.path.exists("ani.list"):
	with open("ani.list", "w") as w:
		[w.write(f"{i}\n") for i in fastas]

for fasta in fastas:
	### ani.list
	## 03-IndFasta/LjRoot34_10.fasta
	## ...
	name = os.path.basename(fasta).split(".")[0]
	gwf.target_from_template(f"FastANI_{name}", ani(fasta=fasta, input_list="ani.list", name=name, outdir=f"04-ANI/{name}", threads=4, memory=12))

### Third round

# Extract contigs as individual protein fasta
if os.path.exists(contigs_file):
	contigs2faa(contigs_file, "05-IndProtFasta")

# Pangenome / ntSync
gwf.target_from_template(f"Plasmid_Pangenome", orthofinder(input_list="ani.list", indir="05-IndProtFasta", outdir=f"06-Pangenome", threads=4, memory=32))