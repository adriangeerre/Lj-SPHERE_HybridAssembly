# Imports
import os
import glob

# Files
inPath = "00-ControlData/DownData"
files = glob.glob(f"{inPath}/*/*_assembly_report.txt")

# Loop
outFile = "control_contigs.tsv"
with open(outFile, "w") as output:
    for f in files:
        name = "_".join(os.path.basename(f).split("_")[:2])
        with open(f, "r") as report:
            for line in report:
                if line[0] != "#":
                    contig = line.split()[0]
                    access = line.split()[4]
                    ctype = line.split()[3]
                    output.write(f"{name}\t{access}\t{contig}\t{ctype}\n")
