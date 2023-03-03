# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  	: Adrián Gómez Repollés
# Email			: adrian.gomez@mbg.au.dk
# Created Date	: 03/03/2023
# version 		: '1.0'
# ---------------------------------------------------------------------------
# """ Script to identify identical contigs given threshold""" 
# ---------------------------------------------------------------------------

# Libraries
# ---------
library(tidyverse)
library(optparse)

# Parameters
# ----------
option_list = list(
  make_option(c("-m", "--distance_matrix"), action="store", default=NA, type='character', help="MashTree distance matrix."),
  make_option(c("-s", "--contigs_size"), action="store", default=NA, type='character', help="Length of contigs."),            
  make_option(c("-c", "--cutoff"), action="store", default=NA, type='character', help="Cutoff to define identical contigs.")
  make_option(c("-o", "--output_path"), action="store", default=95, type='character', help="Path to save output.")
)
opt = parse_args(OptionParser(option_list=option_list))

# Load Data
# ---------

# Read file
df <- read.table(opt$m, header=T, sep="\t")
sizes <- read.table(opt$s, header=T, sep="\t")

# Transform data
# --------------

# Rownames
row.names(df) <- df$. 
df <- df[,2:ncol(df)]

# Subset similar contigs
sdf <- df %>% add_rownames() %>% gather(key, value, -rowname) %>% filter(value <= 1-opt$c) %>% filter(rowname != key) %>% spread(key, value, fill = NA)

# Sort contigs
# ------------



# Save file
# ---------


