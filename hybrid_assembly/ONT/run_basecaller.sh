#!/bin/bash
#SBATCH --account CCRP_Data
#SBATCH --partition gpu
#SBATCH --mem 170G
#SBATCH -c 18
#SBATCH --gres=gpu:1
#SBATCH -t 48:00:00

# Basecall (GPU)
cd $1 # fast5 folder should be inside
guppy_basecaller -i fast5 -s basecall -c dna_r9.4.1_450bps_hac.cfg --num_callers 3 --gpu_runners_per_device 6 --device auto --compress_fastq