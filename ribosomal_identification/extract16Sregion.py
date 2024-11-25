# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By    : Adrián Gómez Repollés
# Email         : adrian.gomez@mbg.au.dk
# Created Date  : 19/10/2023
# version       : '1.0'
# ---------------------------------------------------------------------------
""" Pipeline to process the deep untargeted RNA-seq perform in NCS Univ."""
# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
 
import os
import sys
import glob
from gwf import Workflow, AnonymousTarget

# conda activate GWF
gwf = Workflow(defaults={"account": "CCRP_Data"})

# Functions
#----------

# Subset Silva database
def subset_silva(database, taxonomy, prefix, region, outdir, threads, memory):
	inputs = [database, taxonomy]
	outputs = [f"{outdir}/{prefix}_{region}.qza", f"{outdir}/{prefix}_{region}.classifier.qza"]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '4:00:00'}

	# Choose primers
	primers = {"V3V4": ("CCTACGGGNGGCWGCAG","GACTACHVGGGTATCTAATCC"), "V5V7": ("AACMGGATTAGATACCCKG","ACGTCATCCCCACCTTCC")}
	primer1 = primers[region][0]
	primer2 = primers[region][1]

	spec = '''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Create folder
	mkdir -p {outdir}

	# Subset database (classifier)
	conda activate qiime2-2022.11
	qiime feature-classifier extract-reads --i-sequences {database} --p-f-primer {primer1} --p-r-primer {primer2} --p-min-length 300 --p-max-length 500 --o-reads {outdir}/{prefix}_{region}.qza --p-n-jobs {threads}

	# Train classifier
	qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads {outdir}/{prefix}_{region}.qza --i-reference-taxonomy {taxonomy} --o-classifier {outdir}/{prefix}_{region}.classifier.qza

	'''.format(database=database, taxonomy=taxonomy, primer1=primer1, primer2=primer2, region=region, prefix=prefix, outdir=outdir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Extract 16S
def extract_16S(fasta, prefix, reference, outdir, threads, memory):
	inputs = [fasta, reference]
	outputs = [f"{outdir}/{prefix}.tsv",f"{outdir}/{prefix}.fasta"]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '1:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Create folder
	mkdir -p {outdir}

	# Blast
	conda activate Blast
	blastn -query {fasta} -subject {reference} -strand both -outfmt "6 std qseq" > {outdir}/{prefix}.tsv

	# Save blast hits as fasta
	cut -f 13 {outdir}/{prefix}.tsv | sed "s/^/>{prefix}\\x0A/g" | sed 's/-//g' > {outdir}/tmp
	
	# Number sequences
	counter=1
	while IFS= read -r line; do
		if [[ $line == ">"* ]]; then
			new_header=">{prefix}.$counter"
			((counter++))
			echo "$new_header"
		else
			echo "$line"
		fi
	done < {outdir}/tmp > {outdir}/{prefix}.pre-clean.fasta
	rm -rf {outdir}/tmp

	# SeqKit
	conda activate SeqKit
	seqkit rmdup -s {outdir}/{prefix}.pre-clean.fasta > {outdir}/{prefix}.fasta
	'''.format(fasta=fasta, prefix=prefix, reference=reference, outdir=outdir)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# 16S taxonomy
def classify_16S(r16S, database, taxonomy, prefix, outdir, threads, memory):
	inputs = [r16S, database, taxonomy]
	outputs = [f"{outdir}/{prefix}_16S.taxonomy.tsv", f"{outdir}/{prefix}.qza"]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '1:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Run Qiime2
	conda activate qiime2-2022.11

	if [ -f {outdir}/{prefix}.qza ]; then
		rm -rf {outdir}/{prefix}.qza
		qiime tools import --type 'FeatureData[Sequence]' --input-path {r16S} --output-path {outdir}/{prefix}.qza
	else
		qiime tools import --type 'FeatureData[Sequence]' --input-path {r16S} --output-path {outdir}/{prefix}.qza
	fi

	qiime feature-classifier classify-consensus-vsearch --i-query {outdir}/{prefix}.qza --i-reference-reads {database} --i-reference-taxonomy {taxonomy} --p-threads {threads} --output-dir {outdir}/Blast_16S --o-classification {outdir}/{prefix}_16S.taxonomy.qza
	qiime tools export --input-path {outdir}/{prefix}_16S.taxonomy.qza --output-path {outdir}
	mv {outdir}/taxonomy.tsv {outdir}/{prefix}_16S.taxonomy.tsv

	'''.format(r16S=r16S, database=database, taxonomy=taxonomy, prefix=prefix, outdir=outdir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Extract region and taxonomy
def extract_region_and_classify(r16S, prefix, region, classifier, outdir, threads, memory):
	# Choose primers
	primers = {"V3V4": ("CCTACGGGNGGCWGCAG","GACTACHVGGGTATCTAATCC"), "V5V7": ("AACMGGATTAGATACCCKG","ACGTCATCCCCACCTTCC")}
	primer1 = primers[region][0]
	primer2 = primers[region][1]

	# Variables
	inputs = [r16S]
	outputs = [f"{outdir}/{prefix}_{region}.fasta", f"{outdir}/{prefix}_{region}.taxonomy.tsv"]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '1:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Load into Qiime2
	conda activate qiime2-2022.11

	# Extract region
	qiime feature-classifier extract-reads --i-sequences {outdir}/{prefix}.qza --p-f-primer {primer1} --p-r-primer {primer2} --p-min-length 250 --p-max-length 800 --o-reads {outdir}/{prefix}_{region}.qza
	qiime tools export --input-path {outdir}/{prefix}_{region}.qza --output-path {outdir}/{prefix}_{region}

	# Move file
	mv {outdir}/{prefix}_{region}/dna-sequences.fasta {outdir}/{prefix}_{region}.pre-clean.fasta
	rmdir {outdir}/{prefix}_{region}

	# Taxonomy
	qiime feature-classifier classify-sklearn --i-classifier {classifier} --i-reads {outdir}/{prefix}_{region}.qza --o-classification {outdir}/{prefix}_{region}.taxonomy.qza
	qiime tools export --input-path {outdir}/{prefix}_{region}.taxonomy.qza --output-path {outdir}/{region}
	mv {outdir}/{region}/taxonomy.tsv {outdir}/{prefix}_{region}.taxonomy.tsv
	rmdir {outdir}/{region}

	# Clean fasta
	conda activate SeqKit
	seqkit rmdup -s {outdir}/{prefix}_{region}.pre-clean.fasta > {outdir}/{prefix}_{region}.fasta

	'''.format(r16S=r16S, prefix=prefix, primer1=primer1, primer2=primer2, region=region, classifier=classifier, outdir=outdir)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Execution
# ---------

# Variables
file = "genomes.tsv" # Genome1	/path/to/genome1
reference = "target.fna"
database = "01-Db/silva-138-99-seqs.qza"
taxonomy = "01-Db/silva-138-99-tax.qza"

# Check that files exist
for i in [file, reference, database, taxonomy]:
	if not os.path.exists(i):
		print(f"ERROR: {i} is missing.")
		sys.exit() 

# Preliminary steps
gwf.target_from_template('Subset_silva_V3V4', subset_silva(database=database, taxonomy=taxonomy, prefix="silva-138-99", region="V3V4", outdir="01-Db", threads=4, memory=8))
gwf.target_from_template('Subset_silva_V5V7', subset_silva(database=database, taxonomy=taxonomy, prefix="silva-138-99", region="V5V7", outdir="01-Db", threads=4, memory=8))

# Loop
f = open(file,"r")

for line in f:
	# Avoid commented lines
	if line[0] == "#":
		continue
	# Split line
	prefix,fasta = line.strip().split("\t")
	# Extract 16S
	gwf.target_from_template(f'Extract_16S_{prefix}', extract_16S(fasta=fasta, prefix=prefix, reference=reference, outdir=f"02-Regions/{prefix}", threads=1, memory=2))
	# Classify 16S
	gwf.target_from_template(f'Classify_16S_{prefix}', classify_16S(r16S=f"02-Regions/{prefix}/{prefix}.fasta", database=database, taxonomy=taxonomy, prefix=prefix, outdir=f"02-Regions/{prefix}", threads=4, memory=8))
	# Extract region and Taxonomy
	for region in ["V3V4","V5V7"]:
		gwf.target_from_template(f'Extract_classify_{region}_{prefix}', extract_region_and_classify(r16S = f"02-Regions/{prefix}/{prefix}.qza", prefix=prefix, region=region, classifier=f"01-Db/silva-138-99_{region}.classifier.qza", outdir=f"02-Regions/{prefix}", threads=4, memory=24))
# Close file
f.close()