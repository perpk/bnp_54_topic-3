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
 install.packages("dplyr");
# 
library("dplyr");
library("scales");
library("tidyverse");

# 1 - Read data from file into a variable.
rna.seq.data <- read.csv("rna_seq_data.csv");
rna.seq.data <- rna.seq.data %>% remove_rownames %>% column_to_rownames(var = "X")

# 2.
# 2.1 - Find Genes wih na-values.
genes.na.values <- colnames(rna.seq.data)[colSums(is.na(rna.seq.data)) > 0];

# 2.2 - Register the percentage of na values per column and remove the ones with pct>0 from the data frame.
pct.col.na <- sapply(rna.seq.data, function(c) (sum(is.na(c) / nrow(rna.seq.data))));
na.genes <- names(pct.col.na[pct.col.na != 0]);
rna.seq.data.clean <- rna.seq.data[, !(colnames(rna.seq.data) %in% na.genes)];

# 3.
# 3.1 - For each column find the max. expression value.
max.gene.expression <- rna.seq.data.clean %>% summarize_if(is.numeric, max);

# 3.2 - For each column find the min. expression value.
min.gene.expression <- rna.seq.data.clean %>% summarize_if(is.numeric, min);

# 3.3 - For each column find the mean expression value.
mean.gene.expression <- rna.seq.data.clean %>% summarize_if(is.numeric, mean);

# 3.4 - For each column find the median expression value.
median.gene.expression <- rna.seq.data.clean %>% summarize_if(is.numeric, median);

# 3.5 - Add new rows to the dataframe containing all above measures for each column.
rna.seq.data.measured <- rbind(rna.seq.data.clean, Max = max.gene.expression, Min = min.gene.expression, Mean = mean.gene.expression, Median = median.gene.expression);

# 4. - Find how many samples are higher than the median (outliers) for Genes 3 and 15
# 4.1
gene_3.median <- rna.seq.data.measured["Median", "Gene_3"];
gene_3.outliers <- rna.seq.data.clean[rna.seq.data.clean$Gene_3 > gene_3.median,]$Gene_3;

# 4.2
gene_15.median <- rna.seq.data.measured["Median", "Gene_15"];
gene_15.outliers <- rna.seq.data.clean[rna.seq.data.clean$Gene_15 > gene_15.median,]$Gene_15;

# 5. - Create boxplot for each Gene with their corresponding expression values.
boxplot(rna.seq.data.clean, col="cyan")

# 6. - Create bar plots with the Median value of each Gene's expression.
all.median <- rna.seq.data.measured["Median",]
list.median.vals <- unlist(all.median)
barplot(unlist(all.median), 
        horiz = TRUE, 
        las = 1, 
        border = "pink", 
        font = 3, 
        cex.names = 0.8, 
        cex.axis = 0.6, 
        ylim=c(0,20), 
        col=ifelse(unlist(all.median)>=19, "red", "pink"));

# 7. - Segregate data into vectors and Visualize pathological samples in a Pie Chart.

