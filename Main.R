## ---------------------------
##
## Script name: Main.R
##
## Purpose of script:
##
## Author: Konstantinos Perperidis
##
## Date Created: 2024-03-26
##
## Email: k.a.perperidis@gmail.com
##        std531050@ac.eap.gr
##
## ---------------------------

# TODO remove the following commented block
# install.packages("tidyverse");
# install.packages("dplyr");
# 
# library("dplyr");
# library("tidyr");

library("scales");

# 1 - Read data from file into a variable.
rna.seq.data <- read.csv("rna_seq_data.csv");

# 2.1 - Find Genes wih na-values.
genes.na.values <- colnames(rna.seq.data)[colSums(is.na(rna.seq.data)) > 0];

# 2.2 - Register the percentage of na values per column and remove the ones with pct>0 from the data frame.
pct.col.na <- sapply(rna.seq.data, function(c) (sum(is.na(c) / nrow(rna.seq.data))))
na.genes <- names(pct.col.na[pct.col.na != 0])
rna.seq.data.clean <- rna.seq.data[, !(colnames(rna.seq.data) %in% na.genes)]



