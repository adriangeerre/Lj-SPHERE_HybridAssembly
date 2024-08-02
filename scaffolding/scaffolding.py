# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  	: Adrián Gómez Repollés
# Email			: adrian.gomez@mbg.au.dk
# Created Date	: 31/07/2024
# version 		: '1.0'
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

# RagTag
def ragtag(ref_fna, que_fna, outdir, threads, memory):
	# Folder structure
	if os.path.isdir(outdir) == False: os.makedirs(outdir)

	# GWF
	inputs = [ref_fna, que_fna]
	outputs = []
	options = {'cores': '{}'.format(threads),'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '12:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

	# GTDB-tk
	conda activate Scaffold
	ragtag.py scaffold -t {threads} -o {outdir} {ref_fna} {que_fna}
	'''.format(ref_fna=ref_fna, que_fna=que_fna, outdir=outdir, threads=threads)

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
gwf.target('NCBI_download', inputs=[infile], outputs=["02-NCBI_genomes"], cores = 1, memory = '8gb', queue = 'normal', walltime = '12:00:00') << f"python src/ncbi_download.py -i {infile} -o 02-NCBI_genomes"

# Scaffold
#gwf.target_from_template(f"Scaffold_{}", ragtag(ref_fna=reference, que_fna=query, indir=indir, outdir="03-Scaffold", threads=16, memory=64))

# Re-annotate

# Automatic file creation