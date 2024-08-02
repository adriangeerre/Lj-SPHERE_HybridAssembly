## Imports
# --------
import os
import sys
import argparse
import subprocess
import pandas as pd

## Arguments
# ----------

def params():
	# Parser
	parser = argparse.ArgumentParser(prog='ncbi_download.py', description='Download genomes from NCBI using GTDB-Tk output.')

	# Options
	parser.add_argument('-i', '--infile', dest='infile', action='store', help='GTDB tsv file', required=True)
	parser.add_argument('-o', '--outdir', dest='outdir', action='store', help='Output directory', required=True)

	# Args
	args = parser.parse_args()

	return args

## Execution
# ----------

# Parse GTDB
def parse_gtdb(infile):
	# Storage of results
	res = {}
	# Read file
	df = pd.read_table(infile)
	# Subset
	sub = df[['user_genome','fastani_reference','closest_placement_reference','other_related_references(genome_id,species_name,radius,ANI,AF)']]
	# Select corresponding reference genome/s
	for row in sub.iterrows():
		# Define genome name
		genome_name = row[1]['user_genome']
		# Extract closest reference
		if isinstance(row[1]['fastani_reference'], str):
			res[genome_name] = row[1]['fastani_reference']
		elif isinstance(row[1]['closest_placement_reference'], str):
			res[genome_name] = row[1]['closest_placement_reference']
		else:
			# Select strains with 90% or more ANI
			references = row[1]['other_related_references(genome_id,species_name,radius,ANI,AF)']
			# Check if also NaN
			if isinstance(references, str):
				references = [i.split(",")[0].strip() for i in references.split(";") if float(i.split(",")[2]) >= 90]
				references = list(set(references))
				res[genome_name] = ";".join(references)				
	# Return results
	return res

# Download metadata
def download_ncbi_metadata(res, outdir, outfile):
	# Loop
	for genome, GCFs in res.items():
		for GCF in GCFs.split(";"):
			# Define execution
			cmd = '''/bin/bash -c "
			# Conda
			source ~/programas/minconda3.9/etc/profile.d/conda.sh
			conda activate Entrez

			# Download dataset
			printf '{genome};{GCF};$(esearch -db assembly -query "{GCF}" | esummary | xtract -pattern DocumentSummary -sep ";" -element SpeciesName,AssemblyStatus,SpeciesTaxid,RefSeq_category)' >> {outdir}/{outfile}
			'''.format(GCF=GCF, genome=genome, outdir=outdir, outfile=outfile)

			# Execute
			try:
				subprocess.check_call(cmd, shell=True)
			except KeyboardInterrupt:
				sys.exit(0)
			except:
				print(f"Error while processing {GCF}")

# Download data
def download_ncbi_datasets(infile, outdir):
	# Parse GTDB results
	res = parse_gtdb(infile)

	# Loop by isolates
	for genome, GCFs in res.items():
		# Create folder
		if not os.path.isdir(f"{outdir}/{genome}"): os.makedirs(f"{outdir}/{genome}")

		# Loop genomes
		for GCF in GCFs.split(";"):
			if not os.path.exists(f"{outdir}/{genome}/{GCF}/data/{GCF}/protein.faa"):
				# Define execution
				cmd = '''/bin/bash -c "
				# Conda
				source ~/programas/minconda3.9/etc/profile.d/conda.sh
				conda activate ncbi_datasets

				# Download dataset
				datasets download genome accession {GCF} --include genome,protein,cds,gff3,gbff --no-progressbar --filename {outdir}/{genome}/{GCF}.zip
				
				# Unzip
				unzip {outdir}/{genome}/{GCF}.zip -d {outdir}/{genome}/
				mv {outdir}/{genome}/ncbi_dataset {outdir}/{genome}/{GCF}
				rm {outdir}/{genome}/README.md"
				'''.format(GCF=GCF, genome=genome, outdir=outdir)

				# Execute
				try:
					subprocess.check_call(cmd, shell=True)
				except KeyboardInterrupt:
					sys.exit(0)
				except:
					print(f"Error while processing {GCF}")

		# Download data
		outfile = "NCBIs.entrez.metadata.tsv"
		if not os.path.exists("{outdir}/{outfile}"):
			download_ncbi_metadata(res, outdir, outfile)
		
## Execution
# ----------

if __name__ == '__main__':
	# Arguments
	args = params()
	infile = args.infile
	outdir = args.outdir

	# Execute
	download_ncbi_datasets(infile, outdir)