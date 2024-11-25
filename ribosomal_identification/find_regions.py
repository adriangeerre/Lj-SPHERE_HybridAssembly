import glob
from Bio.Seq import Seq

folders = glob.glob("01-Regions/*")

for iso in folders:
	# Variables
	name = iso.split("/")[1]
	# Sequences
	r16S = open(f"{iso}/{name}.fasta","r").readlines()[1].strip()
	v3v4 = open(f"{iso}/{name}_V3V4.fasta","r").readlines()[1].strip()
	v5v7 = open(f"{iso}/{name}_V5V7.fasta","r").readlines()[1].strip()
	# Complementary
	cv3v4 = str(Seq(v3v4).reverse_complement())
	cv5v7 = str(Seq(v5v7).reverse_complement())
	# Map sequence
	if v3v4 in r16S and v5v7 in r16S:
		print(f"{name}\tNon-complementary\n")
	elif cv3v4 in r16S and cv5v7 in r16S:
		print(f"{name}\tComplementary\n")
	else:
		print(f"ERROR: {name} v3v4 or v5v7 not found in any direction.")