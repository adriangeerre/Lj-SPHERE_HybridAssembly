#!/bin/bash
#SBATCH --account CCRP_Data
#SBATCH --partition normal 
#SBATCH --mem 2G
#SBATCH -c 8
#SBATCH -t 12:00:00

## Run demultiplexer (CPU)
cd $1
guppy_barcoder -i basecall/pass/ -s demultiplex --barcode_kits SQK-RBK110-96 --trim_barcodes --num_barcoding_threads 8 --compress_fastq >> demultiplex.log
