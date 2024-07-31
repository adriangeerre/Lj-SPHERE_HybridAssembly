# LjSC Hybrid Assembly

This repository contains the [GWF](https://gwf.app/) pipeline and auxiliary scripts for the assembly improvement of the _Lotus japonicus_ synthetic community genomes. This pipeline was built to run in an High-Performance Computing enviroment (such as, [GenomeDK](https://genome.au.dk/)) using the Slurm system. For detailed information on how this pipeline works please read: [CITATION!!!]

_Disclaimer:_ The manual curation and scaffolding of the genomes was done after running the pipeline. The code could be found in `pipeline/post` and `scaffolding`.

### Installation

Install conda if not yet installed. Then, Install the required conda environments from the `envs` folder. The environment `HA` is the pipeline's environment and, therefore, the only enviroment to be activated. The other environments are auxiliary and automatically activated during the pipeline run.

Moreover, the NCBI's Prokaryotic Genome Annotation Pipeline or [PGAP](https://github.com/ncbi/pgap) needs to be installed in your working enviroment. Warning: HPC enviroments using Slurm might caused PGAP to fail without a clear reason, in that case I recomment patching this issue with a local computer.

### Required files

Few parts of the pipeline are hard-coded causing some files to be required:

- strains.tsv: This is a tab-delimited file that contains five columns: nanopore_reads, illumina_r1, illumina_r2, folder_name, run_coverage (example: 00-Data/TA1/SRR11719982_pass.fastq.gz 00-Data/TA1/SRR3927460_pass_1.fastq.gz 00-Data/TA1/SRR3927460_pass_2.fastq.gz TA1 False). Each row represent one genome and there should be no repeated rows.
- LjSphere_taxonomy.csv: This is a tab-delimited file containing the name of the isolate and the corresponding taxonomy. This taxonomy is provided to two software: Busco and PGAP, and to the validation steps. 

### Adaptating code before usage

Unfortunately, the code is not yet polished to run automatically in new enviroments and computers. This means that the code or folder names should be adjusted to make it work. More exactly:

- Lines: 32, 52, 67* and 75 on HA.py
- All lines on the src scripts defining the conda path `source ~/programas/minconda3.9/etc/profile.d/conda.sh` to fit your environment's path

*Modify this only if you are trying to run on new data and the taxonomic family or busco database is not present.

### Usage

Before using the pipeline, look into what is and how to use [GWF](https://gwf.app/).

```
conda activate HA
gwf -f HA.py status <YOUR_GENOME>_*
gwf -f HA.py run <YOUR_GENOME>_<STEP>
```

Be aware of the large amount of jobs that are created by the pipeline, I do not recommend to run all at once.

### ToDo

- Remove hard-coded busco database and taxonomy within pipeline script
- Remove hard-coded conda environment in auxiliary scripts
- Provide alternative to PGAP if fails during run or lack of local computer

### Acknowledgments

This pipeline was inspired by Benjamin Perry's hybrid assembly workflow ([Github repository](https://github.com/BenjaminJPerry/HybridAssembly)).
