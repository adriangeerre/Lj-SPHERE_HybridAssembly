import os
from Bio import SeqIO
from gwf import AnonymousTarget

# Remove empty contigs in annotation
def remove_empty_contigs(assembly, gff, out_dir):
    # Open GFF
    f = open(gff, "r")
    # Define contigs dict
    contigs = {}
    # Count features per contig
    for line in f:
        if line[0] == "#":
            pass
        else:
            line = line.strip().split("\t")
            if line[1] == "Local":
                contigs[line[0]] = {"gene": 0, "CDS": 0, "pseudogene": 0, "rRNA_5S": 0, "rRNA_16S": 0, "rRNA_23S": 0, "tRNA": 0, "tmRNA": 0}
            else:
                if line[2] in ["gene","CDS","pseudogene","tRNA","tmRNA"]:
                    contigs[line[0]][line[2]] += 1
                elif line[2] == "rRNA":
                    p = line[-1].strip().split(";")
                    p = [j for j in p if j[0:7] == "product"]
                    p = p[0].split("=")[1].split(" ")[0]
                    if p != "":
                        contigs[line[0]][f"rRNA_{p}"] += 1

    # Define spurious contigs
    keep = {} # 0/1 = Erase/Keep
    for n,c in contigs.items():
        if sum(c.values()) == 0:
            keep[n] = 0
        elif sum(c.values()) != 0:
            if c['gene'] == 0 and c['CDS'] == 0 and c['pseudogene'] == sum(c.values()):
                keep[n] = 0
            else:
                keep[n] = 1

    # Define name
    name = assembly.split("/")[-1].split(".")[0]

    # Clean assembly
    fasta = SeqIO.to_dict(SeqIO.parse(assembly, "fasta"))
    new = [s for i,s in fasta.items() if keep[i] != 0]
    SeqIO.write(new, f'{out_dir}/{name}.clean.fasta', "fasta")

    return True

## GWF function
# Duplicate contigs
def remove_duplicate_contigs(assembly, prefix, threshold, out_dir, threads, memory):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# GWF
	inputs = ["{}".format(assembly)]
	outputs = ["{}/{}.clean-overlap.fasta".format(out_dir, prefix)]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'short', 'walltime': '4:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

    if [ $(grep "^>" {assembly} | wc -l) -gt 1 ]; then
    
        # Split assembly in contigs
        conda activate SeqKit
        seqkit split -i {assembly} -O {out_dir}/{prefix}

        # MashTree
        conda activate MashTree
        mashtree --numcpus {threads} --outmatrix {out_dir}/{prefix}.distmat --outtree {out_dir}/{prefix}.tree --mindepth 0 {out_dir}/{prefix}/*.fasta

        # Size of contigs
        for fna in $(ls -d {out_dir}/{prefix}/*); do printf "$(basename $fna | sed 's/.fasta//g')\t$(cat $fna | grep -v "^>" | tr -d "\n" | wc -c)\n"; done | sort -rn -k 2 > {out_dir}/contigs_size.tsv

        # Summarize all
        conda activate Renv
        Rscript src/cleaning.R -m {out_dir}/{prefix}.distmat -s {out_dir}/contigs_size.tsv -c {threshold} -o {out_dir}/

        # Clean assembly
        conda activate NGSTools
        cat {out_dir}/remove_contigs.txt | cut -d "." -f 2 | cut -d "_" -f 2 | sed 's/^/>/' | sort | uniq > {out_dir}//tmp.txt
        grep "^>" {assembly} | grep -w -v -f {out_dir}/tmp.txt | sed 's/^>//g' > {out_dir}/keep_contigs.txt
        seqtk subseq {assembly} {out_dir}/keep_contigs.txt > {out_dir}/{prefix}.clean-overlap.fasta
        rm {out_dir}/tmp.txt
    
    else
        cp {assembly} {out_dir}/{prefix}.clean-overlap.fasta
    fi

	'''.format(assembly=assembly, prefix=prefix, threshold=threshold, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
