from Bio import SeqIO

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
def remove_duplicate_contigs(assembly, threshold, out_dir, threads, memory):
	# Folder structure
	if os.path.isdir(out_dir) == False: os.makedirs(out_dir)

	# GWF
	inputs = ["{}".format(assembly)]
	outputs = ["{}/{}.clean.fasta".format(out_dir, assembly)]
	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'short', 'walltime': '4:00:00'}

	spec='''
	# Source conda to work with environments
	source ~/programas/minconda3.9/etc/profile.d/conda.sh

    # Split assembly in contigs
    conda activate SeqKit
    seqkit split -i {assembly} -o {out_dir}

    # FastANI
    conda activate FastANI
    for i in $(ls {out_dir}); do
        mkdir -p {out_dir}/${i}
        for j in $(ls {out_dir}); do
            if [[ ! -f {out_dir}/${i}/${i}_${j}.ani && ! -f {out_dir}/${j}/${j}_${i}.ani ]]
                fastANI -q {out_dir}/${i} -r {out_dir}/${j} -o {out_dir}/${i}/${i}_${j}.ani -t {threads}
        done
    done

    # Summarize all
    cat {out_dir}/*/*.ani | awk '($1 != $2)' | awk '($3 >= {threshold})' > {out_dir}/

	'''.format(assembly=assembly, out_dir=out_dir, threads=threads)

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# cd folder
# cat * | awk '($1 != $2)'

# library(tidyverse)
# df <- read.table("LjRoot178.distmat", header=T, sep="\t")
# row.names(df) <- df$. 
# df <- df[,2:ncol(df)]
# df %>% add_rownames() %>% gather(key, value, -rowname) %>% filter(value <= 0.05) %>% filter(rowname != key) %>% spread(key, value, fill = NA)
# LjRoot34!!!