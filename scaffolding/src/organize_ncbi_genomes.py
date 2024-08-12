## Imports
# --------
import os
import glob
import shutil
import argparse

## Arguments
# ----------

def params():
	# Parser
	parser = argparse.ArgumentParser(prog='ncbi_download.py', description='Download genomes from NCBI using GTDB-Tk output.')

	# Options
	parser.add_argument('-i', '--indir', dest='indir', action='store', help='Directory containing the unzip NCBI data', required=True)
	parser.add_argument('-o', '--outdir', dest='outdir', action='store', help='Output directory', required=True)

	# Args
	args = parser.parse_args()

	return args

## Functions
# ----------

def copy_genome(gpath, outdir):
	# Define new name
	new_name = os.path.basename(os.path.dirname(gpath))
	# Copy file
	if not os.path.exists(f"{outdir}/{new_name}.fna"):
		shutil.copy2(gpath, f"{outdir}/{new_name}.fna")

def organize_genomes(indir, outdir):
	# Find files
	fnas_one = glob.glob(f"{indir}/*/*/*/data/*/*.fna")
	fnas_two = glob.glob(f"{indir}/*/*/data/*/*.fna")
	fnas = fnas_one + fnas_two
	# Filter files
	fnas = [i for i in fnas if "cds_from_genomic.fna" not in i]
	# Copy and rename files
	[copy_genome(i, outdir) for i in fnas]
	
## Execution
# ----------

if __name__ == '__main__':
	# Arguments
	args = params()
	indir = args.indir
	outdir = args.outdir

	# Create folder
	if os.path.isdir(outdir) == False: os.makedirs(outdir)

	# Execute
	organize_genomes(indir, outdir)
