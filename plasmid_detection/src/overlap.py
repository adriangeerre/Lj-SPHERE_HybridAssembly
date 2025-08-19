# Import
import os
import csv
import argparse

## Functions

# Read input
def ReadInput(input, software):
    # geNomad
    if software == "geNomad":
        f = open(input,"r")
        contigs = [line.split()[0] for line in f if line.split()[0] != "seq_name"]
        f.close()
        return contigs
    # PLASMe
    if software == "PLASMe":
        f = open(input,"r")
        contigs = [line.split()[0] for line in f if line.split()[0] != "contig"]
        f.close()
        return contigs
    # viralVerify
    if software == "viralVerify":
        f = open(input,"r")
        contigs = [line.split(",")[0] for line in f if line.split(",")[0] != "Contig name" and line.split(",")[1] == "Plasmid"]
        f.close()
        return contigs
    # Plasmer
    if software == "Plasmer":
        f = open(input,"r")
        contigs = [line.split()[0] for line in f if line.split()[1] == "plasmid"]
        f.close()
        return contigs
    # PlasmidHunter
    if software == "PlasmidHunter":
        f = open(input,"r")
        contigs = [line.split()[0] for line in f if line.split()[1] == "1.0" and line.split()[0] != "Prediction"]
        f.close()
        return contigs

# Overlap
def Overlap(folder, name):
    # Inputs
    genomad = f"01-results/geNomad/{folder}/{folder}_summary/{folder}_plasmid_summary.tsv"
    plasme = f"01-results/PLASMe/{folder}/results_report.csv"
    viralverify = f"01-results/viralVerify/{folder}/{folder}_result_table.csv"
    plasmer = f"01-results/Plasmer/{folder}/results/{name}.plasmer.predClass.tsv"
    plasmidhunter = f"01-results/PlasmidHunter/{folder}/predictions.tsv"

    # Outputs
    genomad = ReadInput(genomad, "geNomad")
    plasme = ReadInput(plasme, "PLASMe")
    viralverify = ReadInput(viralverify, "viralVerify")
    plasmer = ReadInput(plasmer, "Plasmer")
    plasmidhunter = ReadInput(plasmidhunter, "PlasmidHunter")

    # Overlap (0/1)
    plasmid_contigs = list(set(genomad + plasme + viralverify + plasmer + plasmidhunter))
    plasmid_contigs.sort()
    table = []
    table.append(plasmid_contigs)
    for i, lst in enumerate([genomad, plasme, viralverify, plasmer, plasmidhunter], start=1):
        table.append([1 if num in lst else 0 for num in plasmid_contigs])
    table = [list(row) for row in zip(*table)]

    # Columns
    cols = ['contig', 'genomad', 'plasme', 'viralverify', 'plasmer', 'plasmidhunter']
    table.insert(0, cols)

    # Save file
    outdir = "02-Overlap"
    outpath = f"{outdir}/{folder}.overlap.tsv"
    if not os.path.isdir(outdir): os.makedirs(outdir) 
    with open(outpath, 'w', newline='') as tsv:
        w = csv.writer(tsv, delimiter='\t')
        w.writerows(table)

## Execution

# Arguments
parser = argparse.ArgumentParser(description='Plasmid overlap script.')
parser.add_argument('-f','--folder', help='Name of strain', required=True)
parser.add_argument('-n','--name', help='Name of strain including type sufix', required=True)
args = parser.parse_args()

# Run
if __name__ == "__main__":
    try:
        Overlap(folder = args.folder, name = args.name)
    except:
        print("ERROR: overlap failed!")
