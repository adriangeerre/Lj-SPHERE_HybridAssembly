## Imports
import os
import statistics as st # Python v3.4 or above
import subprocess

## Auxiliary functions
# Busco
def validate_busco(bd):
	# Dict
	r = {"c": 0, "s": 0, "d": 0, "f": 0, "m": 0}
	# Values
	if (bd["C"] / int(bd["dataset_total_buscos"])) * 100 > 90: r["c"] = 1
	if (bd["S"] / int(bd["dataset_total_buscos"])) * 100 > 90: r["s"] = 1
	if (bd["D"] / int(bd["dataset_total_buscos"])) * 100 < 5: r["d"] = 1
	if (bd["F"] / int(bd["dataset_total_buscos"])) * 100 < 5: r["f"] = 1
	if (bd["M"] / int(bd["dataset_total_buscos"])) * 100 < 5: r["m"] = 1
	# Check
	if sum(list(r.values())) == 5:
		return "Pass"
	else:
		return "Failed"

#  CheckM
def validate_checkm(cm):
	# Dict
	r = {"Comp": 0, "Cont": 0}
	# Read file
	f = open(cm, "r")
	l = f.readlines()[1].strip().split()
	f.close()
	# Values
	if float(l[-3]) > 90: r["Comp"] = 1
	if float(l[-2]) < 5: r["Cont"] = 1
	# Check
	if sum(list(r.values())) == 2:
		return "Pass"
	else:
		return "Failed"

# PGAP
# IMPLEMENT NO GENE OVERLAP! (AGAT)
def gff_to_dict(gff):
	# Dict
	features = {}
	# Read file
	f = open(gff, "r")
	# Loop entries
	tlength = 0
	entries = {"gene": [], "CDS": [], "pseudogene": [], "rRNA_5S": [], "rRNA_16S": [], "rRNA_23S": [], "tRNA": []}
	for i in f:
		# Avoid header info
		if i[0] == "#":
			continue	
		i = i.split("\t")
		# Obtain total length
		if i[2] == "region":
			tlength += int(i[4])
		elif i[2] == "rRNA":
			p = i[-1].strip().split(";")
			p = [j for j in p if j[0:7] == "product"]
			p = p[0].split("=")[1].split(" ")[0]
			if p != "":
				entries["rRNA_{}".format(p)].append(abs(int(i[4])-int(i[3])))
		else:
			# Count defined features
			if i[2] in entries.keys():
				entries[i[2]].append(abs(int(i[4])-int(i[3])))
	# Include total length
	entries["tlength"] = [tlength]
	return entries

def validate_pgap(gff):
	# Dict
	r = {"GeneRatio": 0, "PseudoRatio": 0, "tRNA": 0, "rRNA": 0}
	ribo = {"rRNA_5S": 0, "rRNA_16S": 0, "rRNA_23S": 0}
	# Summary
	d = gff_to_dict(gff)
	cd = {i: len(j) for i,j in d.items()} # Total genes = genes + pseudogene
	# Pre-knowledge (E. coli based!)
	l5S = 120
	l16S = 1542
	l23S = 2904
	# Ribosomal Value (at least one and median not deviating more than 10%)
	if cd["rRNA_5S"] >= 1 and abs(st.median(d["rRNA_5S"]) - l5S) <= l5S * 0.1: ribo["5S"] = 1
	if cd["rRNA_16S"] >= 1 and abs(st.median(d["rRNA_16S"]) - l16S) <= l16S * 0.1: ribo["16S"] = 1
	if cd["rRNA_23S"] >= 1 and abs(st.median(d["rRNA_23S"]) - l23S) <= l23S * 0.1: ribo["23S"] = 1
	# Values
	tgenes = cd["gene"] + cd["pseudogene"]
	gene_ratio = (tgenes * 1000) / d['tlength'][0]
	if gene_ratio > 0.9: r["GeneRatio"] = 1
	if cd["pseudogene"] / tgenes < 0.2: r["PseudoRatio"] = 1
	if cd["tRNA"] >= 20: r["tRNA"] = 1
	if sum(list(ribo.values())) == 3: r["rRNA"] = 1
	# Check
	if sum(list(r.values())) == 4:
		return "Pass"
	else:
		return "Failed"

# Quast
def quast(assembly, out_dir, threads, conda_path, logfile):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	cmd='''
	/bin/bash -c "
	# Source conda to work with environments
	source {conda_path}/etc/profile.d/conda.sh

	# Quast
	conda activate Quast
	quast -o {out_dir} -t {threads} {assembly}"
	'''.format(assembly=assembly, out_dir=out_dir, threads=threads, conda_path=conda_path)

	# Exec and log
	f = open(logfile, "w")
	subprocess.check_call(cmd, shell=True, stdout=f, stderr=f)
	f.close()

# Validation
def assembly_validation(assembly, database, out_dir, threads, conda_path, logfile):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	cmd='''/bin/bash -c "
	# Source conda to work with environments
	source {conda_path}/etc/profile.d/conda.sh

	# BUSCO
	conda activate Busco
	busco -m genome -i {assembly} -o {out_dir}/Busco -l {database} -f

	# CheckM
	conda activate CheckM
	checkm lineage_wf -x fasta -t {threads} {assembly_folder} {out_dir}/CheckM --reduced_tree
	checkm qa {out_dir}/CheckM/lineage.ms {out_dir}/CheckM -f {out_dir}/CheckM/results.tsv --tab_table -t {threads}"
	'''.format(assembly=assembly, assembly_folder=assembly.strip("assembly.fasta"), database=database, out_dir=out_dir, threads=threads, conda_path=conda_path)

	# Exec and log
	f = open(logfile, "w")
	subprocess.check_call(cmd, shell=True, stdout=f, stderr=f)
	f.close()
