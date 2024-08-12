# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  	: Adrián Gómez Repollés
# Email			: adrian.gomez@mbg.au.dk
# Created Date	: 31/07/2024
# version 		: '1.0'
# conda env		: HA
# ---------------------------------------------------------------------------
""" Pipeline to scaffold the resulting LjSC high quality genomes and prepare 
them for data upload.""" 
# ---------------------------------------------------------------------------
# Imports 
# ---------------------------------------------------------------------------

# External
import os
import glob
from gwf import Workflow
from gwf import AnonymousTarget

# GWF configuration
gwf = Workflow(defaults={"account": "CCRP_Data"})

# Target functions
#-----------------

# GTDB-Tk
def gtdbtk(prefix, indir, outdir, extension, threads, memory):
	# Folder structure
	if os.path.isdir(outdir) == False: os.makedirs(outdir)

	# GWF
	inputs = glob.glob(f"{indir}/*.{extension}")
	outputs = [f"{outdir}/{prefix}.bac120.summary.tsv"]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '12:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# GTDB-tk
	conda activate gtdbtk-2.1.1
	gtdbtk classify_wf --genome_dir {indir} --out_dir {outdir} -x {extension} --prefix {prefix} --cpus {threads} --pplacer_cpus {threads} --write_single_copy_genes --keep_intermediates
	'''.format(prefix=prefix, indir=indir, outdir=outdir, extension=extension, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# GTDB-Tk
def quast(infiles, indir, outdir, extension, threads, memory):
	# Folder structure
	if os.path.isdir(outdir) == False: os.makedirs(outdir)

	# GWF
	inputs = infiles
	names = [os.path.basename(i).replace(f".{extension}", "") for i in infiles]
	outputs = [f"{outdir}/{i}/transposed_report.tsv" for i in names]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '12:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# Quast
	conda activate Quast
	for fasta in $(ls {indir} | grep {extension}); do
		name=$(echo $fasta | sed 's/.{extension}//g')
		if [ ! -f {outdir}/$name/report.txt ]; then
			mkdir -p {outdir}/$name
			quast -o {outdir}/$name -t {threads} {indir}/$fasta
		fi
	done
	'''.format(indir=indir, outdir=outdir, extension=extension, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# RagTag: Individual
def ragtag_indiv(ref_fna, que_fna, ref_indir, prefix, extension, outdir, threads, memory):
	# Folder structure
	if os.path.isdir(outdir) == False: os.makedirs(outdir)

	# GWF
	inputs = [que_fna, f"{ref_indir}/{ref_fna}.{extension}"]
	outputs = [f"{outdir}/{prefix}/ragtag.scaffold.agp"]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '12:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# GTDB-tk
	conda activate Scaffold
	ragtag.py scaffold -t {threads} -o {outdir}/{prefix} {ref_indir}/{ref_fna}.{extension} {que_fna} -u
	'''.format(ref_fna=ref_fna, que_fna=que_fna, ref_indir=ref_indir, prefix=prefix, extension=extension, outdir=outdir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# RagTag: Multiple
def ragtag_multi(ref_fnas, que_fna, ref_indir, prefix, extension, outdir, threads, memory):
	# Folder structure
	if os.path.isdir(outdir) == False: os.makedirs(outdir)

	# GWF
	inputs = [que_fna] + [f"{ref_indir}/{i}.{extension}" for i in ref_fnas.split(",")]
	outputs = [f"{outdir}/{prefix}/ragtag.merge.agp"]
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '12:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# GTDB-tk
	conda activate Scaffold
	for fna in $(echo {ref_fnas} | tr "," " "); do
		ragtag.py scaffold -t {threads} -o {outdir}/{prefix}_vs_$fna {ref_indir}/$fna.{extension} {que_fna} -u
	done
	ragtag.py merge {que_fna} {outdir}/{prefix}_*/*.agp -o {outdir}/{prefix} -u
	'''.format(ref_fnas=ref_fnas, que_fna=que_fna, ref_indir=ref_indir, prefix=prefix, extension=extension, outdir=outdir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Execution
#----------

# Variables
extension = "fna"
indir = "00-Data/FNAs"

# Fasta files
faas = glob.glob(f"{indir}/*.{extension}")

# GTDB
outdir = "01-GTDB"
prefix = f"LjSC_{len(faas)}"
gwf.target_from_template(f"GTDBtk", gtdbtk(prefix=prefix, indir=indir, outdir=outdir, extension=extension, threads=16, memory=64))

# Download closest strain
infile = f"{outdir}/{prefix}.bac120.summary.tsv"
gwf.target('NCBI_download', inputs=[infile], outputs=["02-NCBI_genomes", "02-NCBI_genomes/NCBIs.entrez.metadata.tsv"], cores = 1, memory = '8gb', queue = 'normal', walltime = '24:00:00') << f"python src/ncbi_download.py -i {infile} -o 02-NCBI_genomes"

# Organize Genomes in new folder
fnas_one = glob.glob(f"02-NCBI_genomes/*/*/*/data/*/*.fna")
fnas_two = glob.glob(f"02-NCBI_genomes/*/*//data/*/*.fna")
fnas = fnas_one + fnas_two
fnas = [i for i in fnas if "cds_from_genomic.fna" not in i]
outfnas = list(set([f"03-NCBI_FNAs/{os.path.basename(os.path.dirname(i))}.fna" for i in fnas]))
gwf.target('NCBI_organize', inputs = fnas, outputs = outfnas, cores = 1, memory = '4gb', queue = 'normal', walltime = '8:00:00') << f"python src/organize_ncbi_genomes.py -i 02-NCBI_genomes -o 03-NCBI_FNAs"

# Genome quality: NCBI
gwf.target_from_template(f"Genome_quality_NCBI", quast(infiles = outfnas, indir="03-NCBI_FNAs", outdir="04-Quast_NCBI", extension=extension, threads=4, memory=16))

# Genome quality: HA
ha_fnas = glob.glob("00-Data/FNAs/*.fna")
gwf.target_from_template(f"Genome_quality_HA", quast(infiles=ha_fnas, indir=indir, outdir="05-Quast_HA", extension=extension, threads=4, memory=16))

# Compare genome to references
quast_ha = glob.glob("05-Quast_HA/*/transposed_report.tsv")
quast_ncbi = glob.glob("04-Quast_NCBI/*/transposed_report.tsv")
gwf.target('Compare_genomes', inputs = quast_ha + quast_ncbi, outputs = ["genome_comparison.tsv"], cores = 1, memory = '8gb', queue = 'normal', walltime = '4:00:00') << f"python src/compare_genomes.py -i 02-NCBI_genomes -a 05-Quast_HA -n 04-Quast_NCBI"

# Continue after manual curation
# ------------------------------

if os.path.exists("ncbi_genomes_for_scaffold.tsv"):

	# Read file
	f = open("ncbi_genomes_for_scaffold.tsv")
	lines = f.readlines()
	f.close()

	# Loop dictionary
	for line in lines:
		genome,reference = line.strip().split(";")
		# Define run type
		if len(reference.split(",")) == 1:
			# Scaffold individual
			gwf.target_from_template(f"Scaffold_indiv_{genome}", ragtag_indiv(ref_fna=reference, que_fna=f"00-Data/FNAs/{genome}.{extension}", ref_indir="03-NCBI_FNAs", prefix = os.path.basename(genome), extension=extension, outdir="06-Scaffold", threads=4, memory=12))
		else:
			# Scaffold multiple
			gwf.target_from_template(f"Scaffold_multi_{genome}", ragtag_multi(ref_fnas=reference, que_fna=f"00-Data/FNAs/{genome}.{extension}", ref_indir="03-NCBI_FNAs", prefix = os.path.basename(genome), extension=extension, outdir="06-Scaffold", threads=4, memory=12))