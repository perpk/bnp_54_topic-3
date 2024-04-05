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

install.packages("dplyr");
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
# 4.1 - Calc. for Gene 3
gene_3.median <- rna.seq.data.measured["Median", "Gene_3"];
gene_3.outliers <- rna.seq.data.clean[rna.seq.data.clean$Gene_3 > gene_3.median,]$Gene_3;

# 4.2 - Calc. for Gene 15
gene_15.median <- rna.seq.data.measured["Median", "Gene_15"];
gene_15.outliers <- rna.seq.data.clean[rna.seq.data.clean$Gene_15 > gene_15.median,]$Gene_15;

print(sprintf("Gene 3  - median = %s ; outliers=[%s]", gene_3.median, paste(gene_3.outliers, collapse=",")));
print(sprintf("Gene 15 - median = %s ; outliers=[%s]", gene_15.median, paste(gene_15.outliers, collapse=",")));

# print(sprintf("The total G-C Content for the given sequence 
#  equals to: %s", percent(calcGCContent(dnaSequence), 
#                          accuracy = 0.2)), quote = FALSE);

# 5. - Create boxplot for each Gene with their corresponding expression values.
boxplot(rna.seq.data.clean, col="cyan", cex.axis=0.8)

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
# 7.1 - Create a vector with the results
results <- c(rep("disease", 12), rep("health", 18))
# 7.2 - Create yet another vector based on the results by providing captions which correspond to the rownames of the cleaned-up DF.
names(results) <-row.names(rna.seq.data.clean);
# 7.3 - Create the PieChart - consolidate the vector as a table which sums up results per category
pie(table(results), col=c("red","green"), border="white", radius=1.06, cex=0.6);
