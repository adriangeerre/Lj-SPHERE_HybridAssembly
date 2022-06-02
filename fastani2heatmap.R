# Libraries
library(tidyverse)

## Data
# List files
f <- list.files("90-FastANI/", recursive=T, full.names = T)
cols <- c("Nanopore", "Illumina", "ANI", "Bidirectional_mapped_fragments", "Total_query_fragments")

# Read files
dfs <- lapply(f, function(x) read.delim(x, sep="\t", header=F))

## Transform
# Set colnames
dfs <- lapply(dfs, setNames, nm = cols)

# Modify path to strain
dfs <- lapply(dfs, function(x) x %>% mutate(Nanopore = unlist(str_split(Nanopore,"/"))[2], Illumina = str_replace_all(str_replace_all(Illumina, "/home/agomez/CCRP_Data/AGR/Data/Illumina_genomes/assemblies/", ""), ".fna", "")) %>% select(Nanopore, Illumina, ANI))

# Reduce
df <- Reduce(function(...) merge(..., all=T), dfs)

# Match Nanopore and Illumina strains
df <- df %>% filter(Illumina %in% df$Nanopore)

# Matrix
m <- pivot_wider(df, names_from = Illumina, values_from = ANI)
m[is.na(m)] <- 0
m <- as.data.frame(m)
rownames(m) <- m$Nanopore
m <- m %>% select(-Nanopore)
m <- as.matrix(m)

## Plot
# Basic Heatmap
heatmap(m, ylab = "Nanopore", xlab="Illumina", scale = "none", Rowv = NA, Colv = NA)

# Basic Ggplot2 (First add missing comparisons)
all <- unique(df$Nanopore)
for (i in df$Nanopore) {
  x <- (df %>% filter(Nanopore == i))$Illumina
  x <- all[!all %in% x]
  for (j in x) {
    df[nrow(df) + 1,] = c(i,j,70)
  }
}

df$ANI <- as.numeric(df$ANI)
df %>% ggplot() + geom_tile(aes(x=Illumina, y=Nanopore, fill=ANI)) +
  scale_fill_gradient2(
    low = "white",
    mid = "white",
    high = "darkred",
    midpoint = 70,
    breaks = c(0, 70, 80, 90, 100),
    limits = c(0,100),
    labels=format(c(0, "70 (min)", 80, 90, "100 (max)")),
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  ) +
  labs(title = "Kmer Assembly Comparison") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
