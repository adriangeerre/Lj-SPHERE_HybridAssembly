## Imports
# --------
import os
import glob
import argparse
import pandas as pd

## Arguments
# ----------

def params():
	# Parser
	parser = argparse.ArgumentParser(prog='ncbi_download.py', description='Download genomes from NCBI using GTDB-Tk output.')

	# Options
	parser.add_argument('-i', '--indir', dest='indir', action='store', help='Directory with genome organization', required=True)
	parser.add_argument('-a', '--ha_reports', dest='ha_reports', action='store', help='Quast reports for HA genomes', required=True)
	parser.add_argument('-n', '--ncbi_reports', dest='ncbi_reports', action='store', help='Quast reports for NCBI genomes', required=True)

	# Args
	args = parser.parse_args()

	return args

## Functions
# ----------

def define_couples(indir):
	# Define HA genomes
	has = glob.glob(f"{indir}/*")
	has = [i for i in has if not i.endswith(".tsv")]
	# Couples
	couples = {os.path.basename(i): [os.path.basename(j).replace(".zip","") for j in glob.glob(f"{i}/*.zip")] for i in has}
	# Return
	return couples

def compare_couple(query, reference):
	'''
	Comparison: Larger, but less fragmented is assumed to be more completed!
	The comparison is made so we distinguish if our genomes are better than the reference.
	If the reference is on the same quality, we could still use it
	'''
	# Read query
	df_query = pd.read_table(query)
	# Read reference
	df_reference = pd.read_table(reference)
	# New row
	row = [df_query['Assembly'].to_list()[0], df_reference['Assembly'].to_list()[0], int(df_query['Total length']), int(df_reference['Total length']), int(df_query['# contigs']), int(df_reference['# contigs']), int(df_query['N50']), int(df_reference['N50']), int(df_query['L50']), int(df_reference['L50'])]
	# Compare length
	if int(df_query['Total length']) < int(df_reference['Total length']):
		row.append(1)
	else:
		row.append(0)
	# Compare Contigs
	if int(df_query['# contigs']) < int(df_reference['# contigs']):
		row.append(1)
	else:
		row.append(0)
	# Compare N50
	if int(df_query['N50']) > int(df_reference['N50']):
		row.append(1)
	else:
		row.append(0)
	# Compare L50
	if int(df_query['L50']) < int(df_reference['L50']):
		row.append(1)
	else:
		row.append(0)
	# Return result
	return(row)

def compare_genomes(indir, quast_ha, quast_ncbi):
	# Results
	table = []
	# Define couples
	couples = define_couples(indir)
	# Compare genomes
	for i,j in couples.items():
		# Load HA report
		query = f"{quast_ha}/{i}/transposed_report.tsv"
		# Load NCBI reports
		if len(j) > 1:
			for a in j:
				if os.path.exists(f"{quast_ncbi}/{a}/transposed_report.tsv"):
					reference = f"{quast_ncbi}/{a}/transposed_report.tsv"
					# Compare
					res = compare_couple(query, reference)
					table.append(res)
				else:
					res = [i, a, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"]
					table.append(res)
					print(f"Multi: {quast_ncbi}/{a}/transposed_report.tsv is missing!")
					continue
		else:
			if os.path.exists(f"{quast_ncbi}/{j[0]}/transposed_report.tsv"):
				reference = f"{quast_ncbi}/{j[0]}/transposed_report.tsv"
				# Compare
				res = compare_couple(query, reference)
				table.append(res)
			else:
				res = [i, j[0], "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"]
				table.append(res)
				print(f"Indiv: {quast_ncbi}/{j[0]}/transposed_report.tsv is missing!")
				continue
	# Create table
	table = pd.DataFrame(table)
	# Rename columns
	cols = ["query", "reference", "query_length", "reference_length", "query_contigs", "reference_contigs", "query_N50", "reference_N50", "query_L50", "reference_L50", "compared_length", "compared_contig", "compared_N50", "compared_L50"]
	table.columns = cols
	# Save table
	table.to_csv("genome_comparison.tsv", index=False, sep="\t")

## Execution
# ----------

if __name__ == '__main__':
	# Arguments
	args = params()
	indir = args.indir
	ha_reports = args.ha_reports
	ncbi_reports = args.ncbi_reports

	# Execute
	compare_genomes(indir, ha_reports, ncbi_reports)