# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  	: Adrián Gómez Repollés
# Email			: adrian.gomez@mbg.au.dk
# Created Date	: 03/03/2023
# version 		: '1.0'
# ---------------------------------------------------------------------------
# """Script to identify identical contigs given threshold""" 
# ---------------------------------------------------------------------------

# Libraries
# ---------
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))

# Parameters
# ----------
option_list = list(
  make_option(c("-m", "--distance_matrix"), action="store", default=NA, type='character', help="MashTree distance matrix."),
  make_option(c("-s", "--contigs_size"), action="store", default=NA, type='character', help="Length of contigs."),            
  make_option(c("-c", "--cutoff"), action="store", default=NA, type='numeric', help="Cutoff to define identical contigs (min similarity:0 / max similarity: 100)."),
  make_option(c("-o", "--output_path"), action="store", default=95, type='character', help="Path to save output.")
)
opt = parse_args(OptionParser(option_list=option_list))

# Load Data
# ---------

# Read file
df <- read.table(opt$distance_matrix, header=T, sep="\t")
sizes <- read.table(opt$contigs_size, header=F, sep="\t")

# Transform data
# --------------

# Rownames (distances)
row.names(df) <- df$. 
df <- df[,2:ncol(df)]

# Subset similar contigs
sdf <- df %>% rownames_to_column() %>% gather(key, value, -rowname) %>% filter(value <= 1-(opt$cutoff/100)) %>% filter(rowname != key)
spdf <- sdf %>% spread(key, value, fill = NA)

# Remove duplicates (A/B == B/A)
sdf <- sdf[!duplicated(t(apply(sdf, 1, sort))),]

# Select contig among pair
remove <- c()
for (n in 1:nrow(sdf)) {
    left <- as.character(sdf[n,1])
    right <- as.character(sdf[n,2])
    if (sizes[sizes$V1 == left,]$V2 >= sizes[sizes$V1 == right,]$V2) {
        # Remove left
        remove <- c(remove, left)
    } else {
        # Remove right
        remove <- c(remove, right)
    }
}

# Save file
# ---------

# Contigs and values
write.table(spdf, file=paste(opt$output_path, "/remove_contigs_values.tsv", sep=""), quote=F, row.names=F, col.names=T, sep="\t")

# Contigs to remove
write.table(remove, file=paste(opt$output_path, "/remove_contigs.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
