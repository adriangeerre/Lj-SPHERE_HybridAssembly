import glob
from Bio import SeqIO
from Bio.SeqUtils import GC

fastas = glob.glob("03-IndFasta/*.fasta")

with open("plasmids_len_gc.tsv", "w") as w:
	for fasta in fastas:
		name = fasta.split("/")[1].split(".")[0]
		for contig in SeqIO.parse(fasta, "fasta"):
			w.write(f"{name}\t{len(contig.seq)}\t{GC(contig.seq)}\n")
