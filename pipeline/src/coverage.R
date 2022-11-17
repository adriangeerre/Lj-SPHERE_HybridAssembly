# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  	: Adrián Gómez Repollés
# Email			    : adrian.gomez@mbg.au.dk
# Created Date	: 11/03/2022
# version 		  : '1.0'
# ---------------------------------------------------------------------------
# """ Script to calculate coverage stats and plot contigs coverage for 
# Nanopore and Illumina alignments against the hybrid assembly created using
# the same sequencing data. The input data comes from a sorted BAM process
# with Bedtools genomecov to obtain a .cov file.""" 
# ---------------------------------------------------------------------------

# Libraries
# ---------
library(tidyverse)
library(optparse)

# Parameters
# ----------
option_list = list(
  make_option(c("-i", "--input_path"), action="store", default=NA, type='character',
              help="Input path with .cov files."),
  make_option(c("-o", "--output_path"), action="store", default=NA, type='character',
              help="Output path to save plots and stats.")
)
opt = parse_args(OptionParser(option_list=option_list))

# Load Data
# ---------

# Read file
f <- list.files(opt$i, full.names=T)
f <- f[grepl(".cov",f)]
dfs <- lapply(f, function(x) read.delim(x, sep="\t",header=F))

# Rename inside lists
sqs <- str_replace_all(f,".cov","")
names(dfs) <- sqs

# Renames columns
for (i in sqs) {
  colnames(dfs[[i]]) <- c("contig_id","contig_bp","depth")
  dfs[[i]] <- dfs[[i]] %>% mutate(seqs = i)
}

# Transform data
# --------------

# Position
dfs <- lapply(dfs, function(x) x %>% group_by(contig_id) %>% mutate(pos = 1:n()) %>% mutate(bin = pos %/% 1000))

# Arrange
dfs <- lapply(dfs, function(x) x %>% mutate(max_pos = max(pos)) %>% arrange(desc(max_pos),pos) %>% select(-max_pos))

# Mean depth per bin
dfs <- lapply(dfs, function(x) x %>% group_by(contig_id, bin, seqs) %>% summarise(bin_depth = mean(depth)))

# Merge dataframes into one
df <- Reduce(function(...) merge(..., all=T), dfs)


# Plot coverage
# -------------

# Plot and save
for (i in unique(df$contig_id)) {
  p <- df %>% filter(contig_id == i) %>% ggplot() + geom_line(aes(x=bin, y=bin_depth, color=seqs)) +
    scale_color_manual(values = c("firebrick","cyan4")) +
    ylab("Depth") + xlab("Genome position (Kbps)") + labs(title = "Assembly Coverage", color="") +
    theme_bw(base_size = 24)
  ggsave(paste(opt$o, "/", i,".pdf", sep=""), plot = p, device = "pdf", width = 48, height = 32, units = "cm")
}
