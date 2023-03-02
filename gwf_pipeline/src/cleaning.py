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
            line = line.split()               
            if line[1] == "Local":
                contigs[line[0]] = {"gene": 0, "CDS": 0, "pseudogene": 0, "rRNA_5S": 0, "rRNA_16S": 0, "rRNA_23S": 0, "tRNA": 0, "tmRNA": 0}
            else:
                if line[2] in ["gene","CDS","pseudogene","rRNA_5S","rRNA_16S","rRNA_23S","tRNA","tmRNA"]:
                    contigs[line[0]][line[2]] += 1
    
    # Define spurious contigs
    keep = {} # 0/1 = Erase/Keep
    for n,c in contigs.items():
        if sum(c.keys()) == 0:
            keep[n] = 0
        elif sum(c.keys()) != 0:
            if c['gene'] == 0 and c['CDS'] == 0 and c['pseudogene'] == sum(c.keys()):
                keep[n] = 0
            else:
                keep[n] = 1

    # Clean assembly
    fasta = SeqIO.to_dict(SeqIO.parse(assembly, "fasta"))
    new = {i:s for i,s in fasta.items() if keep[i] != 0}
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
	
	'''.format()

	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
