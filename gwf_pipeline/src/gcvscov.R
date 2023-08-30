# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  	: Adrián Gómez Repollés
# Email			: adrian.gomez@mbg.au.dk
# Created Date	: 30/08/2023
# version 		: '1.0'
# ---------------------------------------------------------------------------
# """Script to plot GC content versus average coverage per contig""" 
# ---------------------------------------------------------------------------

# Libraries
# ---------
suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))
suppressMessages(library(ggrepel))
suppressMessages(library(optparse))

# Parameters
# ----------
option_list = list(
  make_option(c("-a", "--assembly"), action="store", type='character', help="Path to assembly."),
  make_option(c("-p", "--prefix"), action="store", type='character', help="Name of the isolate."),
  make_option(c("-i", "--in_dir"), action="store", type='character', help="Path to folder including coverage."),
  make_option(c("-o", "--out_dir"), action="store", type='character', help="Path to save output.")
)
opt = parse_args(OptionParser(option_list=option_list))

# Load data
# ---------

# Coverage
ill_cov <- read.table(paste(opt$in_dir, "/Illumina.cov", sep=""), header=F, sep="\t")
ont_cov <- read.table(paste(opt$in_dir, "/Nanopore.cov", sep=""), header=F, sep="\t")

# GC, relative depth, genome size and circularity
contig_data <- read.table(paste(opt$in_dir, "/", opt$prefix, ".contig-data.tsv", sep=""), header=F, sep="\t")

# Transform data
# --------------

# Average coverage per contig
ill_cov <- ill_cov %>% group_by(V1) %>% summarise(ill_mean_cov = mean(V3), ill_meadian_cov = median(V3))
ont_cov <- ont_cov %>% group_by(V1) %>% summarise(ont_mean_cov = mean(V3), ont_meadian_cov = median(V3))

# Merge
cov <- merge(ill_cov, ont_cov, by="V1")
df <- merge(contig_data, cov, by="V1")
colnames(df) <- c("Contig","GC","Gsize","Circularity","Ill_mean_cov","Ill_median_cov","Ont_mean_cov","Ont_median_cov")

# Circularity
df$Circularity <- ifelse(df$Circularity == 0, "No", "Yes")

# Plot
# ----

p1 <- df %>% ggplot() +
    geom_point(aes(x=GC, y=Ill_mean_cov, color=Circularity, size=Gsize / 1000), alpha=0.7) +
    geom_text_repel(aes(x=GC, y=Ill_mean_cov, label = Contig), color="black", alpha=0.8, size=4) +
    scale_color_manual(values = c("blue","orange")) +
    labs(x = "GC%", y = "Coverage", color = "Circular", size="Length (Kbps)") +
    theme_classic()

p2 <- df %>% ggplot() +
    geom_point(aes(x=GC, y=Ill_median_cov, color=Circularity, size=Gsize / 1000), alpha=0.7) +
    geom_text_repel(aes(x=GC, y=Ill_median_cov, label = Contig), color="black", alpha=0.8, size=4) +
    scale_color_manual(values = c("blue","orange")) +
    labs(x = "GC%", y = "Coverage", color = "Circular", size="Length (Kbps)") +
    theme_classic()

p3 <- df %>% ggplot() + 
    geom_point(aes(x=GC, y=Ont_mean_cov, color=Circularity, size=Gsize / 1000), alpha=0.7) +
    geom_text_repel(aes(x=GC, y=Ont_mean_cov, label = Contig), color="black", alpha=0.8, size=4) +
    scale_color_manual(values = c("firebrick","forestgreen")) +
    labs(x = "GC%", y = "Coverage", color = "Circular", size="Length (Kbps)") +
    theme_classic()

p4 <- df %>% ggplot() + 
    geom_point(aes(x=GC, y=Ont_median_cov, color=Circularity, size=Gsize / 1000), alpha=0.7) +
    geom_text_repel(aes(x=GC, y=Ont_median_cov, label = Contig), color="black", alpha=0.8, size=4) +
    scale_color_manual(values = c("firebrick","forestgreen")) +
    labs(x = "GC%", y = "Coverage", color = "Circular", size="Length (Kbps)") +
    theme_classic()

# Save
# ----

png(paste(opt$out_dir, "/",opt$prefix, ".Illumina.mean.GCvsCOV.png", sep=""), width=1600, height=1200, res=300)
p1
dev.off()

png(paste(opt$out_dir, "/",opt$prefix, ".Illumina.median.GCvsCOV.png", sep=""), width=1600, height=1200, res=300)
p2
dev.off()

png(paste(opt$out_dir, "/",opt$prefix, ".Nanopore.mean.GCvsCOV.png", sep=""), width=1600, height=1200, res=300)
p3
dev.off()

png(paste(opt$out_dir, "/",opt$prefix, ".Nanopore.median.GCvsCOV.png", sep=""), width=1600, height=1200, res=300)
p4
dev.off()