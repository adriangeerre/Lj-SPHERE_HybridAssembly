## Imports
import os
import json
import argparse
import subprocess

## Options
def params():
    # Parser
    parser = argparse.ArgumentParser(prog='rerun_busco.py', description='Script to run and validate busco using a different taxonomic annotation as the one used in the HA workflow.')
   
    # Arguments
    parser.add_argument('--assembly', dest='assembly', action='store', help='Assembly file in fasta format.', required=True)
    parser.add_argument('--outdir', dest='outdir', action='store', help='Output folder', required=True)
    parser.add_argument('--taxonomy', dest='taxa', action='store', help='New taxonomic level to use.', required=True)

    # Args
    args = parser.parse_args()

    return args

args = params()

## Functions

# Busco database
def define_busco_database(taxonomy):
    # Options
    busco_dict = {'Actinomycetales': 'actinobacteria_class_odb10', 'Flavobacteriales': 'flavobacteriales_odb10', 'Bacillales': 'bacillales_odb10', 'Burkholderiales': 'burkholderiales_odb10', 'Caulobacterales': 'alphaproteobacteria_odb10', 'Rhizobiales': 'rhizobiales_odb10', 'Sphingomonadales': 'sphingomonadales_odb10', 'Pseudomonadales': 'pseudomonadales_odb10', 'Xanthomonadales': 'xanthomonadales_odb10', 'Hyphomicrobiales':'alphaproteobacteria_odb10', 'Enterobacterales': 'enterobacterales_odb10'}
    # Check taxonomy
    if taxonomy not in busco_dict.keys():
        print(f"ERROR: Taxonomy \"{taxonomy}\" not defined for Busco!")
        exit()
    else:
        database = busco_dict[taxonomy]
    return database

# Busco
def run_busco(assembly, out_dir, database):
    # Create folder
    if not os.path.isdir(out_dir): os.makedirs(out_dir)
    # 
    cmd = f"source ~/programas/minconda3.9/etc/profile.d/conda.sh; conda activate Busco; busco -m genome -i {assembly} -o {out_dir}/Busco_rerun -l {database} -f; conda deactivate"
    try:
        subprocess.run(cmd, shell=True, check=True, text=True)
    except:
        print(f"ERROR: Busco failed to run, something went wrong!")
        exit()

# Validate Busco
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
    
## Execution
database = define_busco_database(args.taxa)
run_busco(args.assembly, args.outdir, database)
b = [i for i in os.listdir(f"{args.outdir}/Busco_rerun") if i[-5:] == ".json"][0]
bf = open(f"{args.outdir}/Busco_rerun/{b}")
bd = json.load(bf)
result = validate_busco(bd)
print(f"\nThe taxonomy \"{args.taxa}\" resulted in: {result}")
