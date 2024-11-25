# DiamondDB
# diamond_db = "/home/agomez/CCRP_Data/AGR/Data/UniProt/BacterialDB_diamond/uniprot_bacterial"

## GWF function
# Diamond
# def diamond(out_dir, diamond_db, threads, memory, folder):
# 	# Folder structure
# 	if os.path.isdir("80-Validation") == False: os.mkdir("80-Validation")
# 	if os.path.isdir("80-Validation/" + folder) == False: os.mkdir("80-Validation/" + folder)
# 	if os.path.isdir("80-Validation/" + folder + "/Diamond") == False: os.mkdir("80-Validation/" + folder + "/Diamond")

# 	# GWF
# 	inputs = ["70-Prokka/{}/prokka_{}.faa".format(folder, folder)]
# 	outputs = ["{}/{}.diamond.tsv".format(out_dir, folder)]
# 	options = {'cores': '{}'.format(threads), 'memory': '{}g'.format(memory), 'queue': 'normal', 'walltime': '6:00:00'}
# 	spec='''
# 	# Source conda to work with environments
# 	source ~/programas/minconda3.9/etc/profile.d/conda.sh

# 	# Diamond (db pre-computed)
# 	conda activate Diamond
# 	diamond blastp --query 70-Prokka/{folder}/prokka_{folder}.faa --db {diamond_db} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen sseq qseq --threads {threads} --out {out_dir}/{folder}.diamond.tsv --un {out_dir}/{folder}.diamond.unalign.fasta --unfmt fasta --fast --tmpdir /scratch/$SLURM_JOBID/ --parallel-tmpdir /scratch/$SLURM_JOBID/
# 	'''.format(out_dir=out_dir, diamond_db=diamond_db, folder=folder, threads=threads)

# 	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)