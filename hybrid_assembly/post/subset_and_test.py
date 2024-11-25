## Imports
import os
import argparse
import subprocess
import validation
from Bio import SeqIO

## Options
def params():
    # Parser
    parser = argparse.ArgumentParser(prog='rerun_busco.py', description='Script to run and validate busco using a different taxonomic annotation as the one used in the HA workflow.')
   
    # Arguments
    parser.add_argument('--infile', dest='infile', action='store', help='Tab-delimited file containing minimum 4 columns and maximum 5 columns: Assembly, Contigs (comma divided), Action [Keep/Remove], and Busco database. Optional, the extra column Annotation if --validate_annotation is set to True.', required=True)
    parser.add_argument('--validate_annotation', dest='valannot', action='store_true', help='Set to validate the annotation.', default=False)
    parser.add_argument('--threads', dest='threads', action='store', help='Number of threads for CheckM.', default=1)

    # Args
    args = parser.parse_args()

    return args

args = params()

## Functions
# Subset fasta
def subset_fasta(assembly, contigs, action, name):
    # Load assembly
    d = SeqIO.to_dict(SeqIO.parse(assembly),'fasta')
    # Define contigs
    contigs = contigs.split(",")
    # Subset
    if action == "Keep":
        subset = [j for i,j in d.items() if i in contigs]
    elif action == "Remove":
        subset = [j for i,j in d.items() if i not in contigs]
    # Save file
    folder = f"subset_and_test/{name}"
    fasta = f"subset_and_test/{name}/{name}.subset.fasta"
    if not os.path.isdir(folder): os.makedirs(folder)
    SeqIO.write(subset, fasta , "fasta")
    return fasta

# Subset annotation
def subset_annotation(annotation):
    continue

# Busco
def run_busco(fasta, database, name):
    # Run
    outdir = f"subset_and_test/{name}"
    cmd = f"source ~/programas/minconda3.9/etc/profile.d/conda.sh; conda activate Busco; busco -m genome -i {fasta} -o {outdir}/Busco_subset -l {database} -f; conda deactivate"
    try:
        subprocess.run(cmd, shell=True, check=True, text=True)
        return f"{outdir}/Busco_subset"
    except:
        print(f"ERROR: Busco failed to run, something went wrong!")
        continue

# CheckM
def run_checkm(threads, name):
    # Run
    outdir = f"subset_and_test/{name}"
    cmd = f"source ~/programas/minconda3.9/etc/profile.d/conda.sh; conda activate CheckM; checkm lineage_wf -x fasta -t {threads} {outdir} {outdir}/CheckM --reduced_tree; checkm qa {outdir}/CheckM/lineage.ms {outdir}/CheckM -f {outdir}/CheckM/results.tsv --tab_table -t {threads}"
    try:
        subprocess.run(cmd, shell=True, check=True, text=True)
        return f"{outdir}/CheckM/results.tsv"
    except:
        print(f"ERROR: CheckM failed to run, something went wrong!")
        continue

# Validation
def validate_subset(busco, checkm, annotation = None):
    # Validation dictionary
    val_ha = {}
    # Validate BUSCO
    if os.path.exists(busco):
        b = [i for i in os.listdir(busco) if i[-5:] == ".json"][0]
        bf = open(f"{busco}/{b}")
        bd = json.load(bf)
        val_ha["bv"] = validation.validate_busco(bd)
        bf.close()
    # Validate CheckM
    if os.path.exists(checkm):
        val_ha["cv"] = validation.validate_checkm(checkm)
    # Validate Annotation
    if annotation is not None:
        if os.path.exists(annotation):
            val_ha["av"] = validation.validate_pgap(annotation)
    return val_ha

# Workflow
def workflow(assembly, contigs, action, database, annotation = None):
    # Define name
    name = os.path.basename(assembly).split(".")[0]
    # Run steps
    fasta = subset_fasta(assembly, contigs, action, name)
    busco = run_busco(fasta, database, name)
    checkm = run_checkm(args.threads, name)
    if annotation is not None:
        annotation = subset_annotation(annotation)
        result = validate_subset(busco, checkm, annotation)
    else:
        result = validate_subset(busco, checkm)
    return result

## Execution
if not os.path.exists(args.infile):
    print(f"ERROR: {args.infile} does not exits!")
    exit()

f = open(args.infile, 'r')
for iso in f:
    if iso[0] == "#":
        continue
    else:
        iso = iso.strip().split("\t")
        # Check assembly
        if not os.path.exists(iso[0]):
            print(f"ERROR: {iso[0]} does not exits.")
            continue
        # Check action
        if iso[2] not in ['Keep','Remove']:
            print(f"ERROR: {iso[2]} is not \'Keep\' or \'Remove\'.")
            continue
        # Check annotation
        if not os.path.exists(iso[5]):
            print(f"ERROR: {iso[5]} does not exist.")
        # Report
        print(f"Processing: {iso[0]}")
        # Run
        if args.valannot:
            workflow(iso[0], iso[1], iso[2], iso[4], iso[5])
        else:
            workflow(iso[0], iso[1], iso[2], iso[4])
f.close()

