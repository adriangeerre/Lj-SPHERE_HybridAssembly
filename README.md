# LjSC Hybrid Assembly

This repository contains the [GWF](https://gwf.app/) pipeline and auxiliary scripts for the assembly improvement of the _Lotus japonicus_ synthetic community genomes. This pipeline was built to run in an High-Performance Computing enviroment (such as, [GenomeDK](https://genome.au.dk/)) using the Slurm system. For detailed information on how this pipeline works please read: [CITATION!!!]

_Disclaimer:_ The manual curation and scaffolding of the genomes was done after running the pipeline. The code could be found in `pipeline/post` and `scaffolding`.

### Installation

Install conda if not yet installed. Then, Install the required conda environments from the `envs` folder. The environment `HA` is the pipeline's environment and, therefore, the only enviroment to be activated. The other environments are auxiliary and automatically activated during the pipeline run.

Moreover, the NCBI's Prokaryotic Genome Annotation Pipeline or [PGAP](https://github.com/ncbi/pgap) needs to be installed in your working enviroment. Warning: HPC enviroments using Slurm might caused PGAP to fail without a clear reason, in that case I recomment patching this issue with a local computer.
