# Lineages project
# github script 2023
# Ashleigh Cheyne

install.packages(c("tidyverse", "dplyr", "tidyr", "ggplot2", "ggfortify", "ape", "ggtree",
                   "gplots", "d3heatmap", "pheatmap", "dendextend", "magrittr", "corrplot",
                   "mixOmics", "reshape2", "rowr", "limma", "factoextra", "DESeq2",
                   "ggbeeswarm", "beeswarm", "ggrepel", "ggpubr", "pROC", "biomaRt",
                   "GEOquery", "Martrix", "edgeR", "RColorBrewer", "Hmisc", "ComplexHeatmap"))

install.packages(c("rlang", "vctrs", "cli"))

library(viridis)

library(rlang)
library(vctrs)
library(cli)

library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(beeswarm)
library(VariantAnnotation)
library(dummies)
library(ape)
library(Biostrings)
library(ggrtree)
library(gplots)
library(circlize)
library(d3heatmap)
library(pheatmap)
library(dendextend)
library(magrittr)
library(corrplot)
library(mixOmics)
library(reshape2)
library(rowr)
library(limma)
library(factoextra)
library(DESeq2)
library(ggbeeswarm)
library(ggrepel)
library(pROC)
library(preprocessCore)
library(gridExtra)
library(biomaRt)
library(Biobase)
library(GEOquery)
library(ggpubr)
library(Matrix)
library(edgeR)
library(RColorBrewer)
library(Hmisc)
library(ComplexHeatmap)

# colours used throughout
# "#CC6666" "gray67" "#003366"
# lin1, lin2, lin4

##### Loading SNP files #####
# loads all files in current directory ending with .csv
temp = list.files(pattern="*.ann.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

files <- Sys.glob("*.csv")
files <- files[1:15]

lineage_list <- vector(mode="list", length=length(files))
names(lineage_list) <- files
for (filename in files) {
  # read data
  sample <- read.csv(filename)
  lineage_list[[filename]] <- sample %>% dplyr::select(POS, ALT, REF, QUAL) %>%
    mutate(lin=filename)
  
}

new_lin <- bind_rows(lineage_list)

##### Doubling times #####
double <- read.csv("doubling_times.csv", header=TRUE)

double$Order <- c(3, 1, 4, 2, 5, 6, 8, 7, 9, 13, 12, 14, 10, 11, 15)
reorder_sub <- c(3, 1, 4, 2, 5, 6, 8, 7, 9, 10, 11, 14, 12, 15, 13)

# plot doubling times using double
png("Doubling_times_paper.png", height = 750, width = 800)
ggplot(double, aes(x = reorder(Strain, Order), y = Hour, color = as.character(Lineage))) + geom_point(pch = "-", size = 7) +
  geom_errorbar(aes(ymin=Hour-SD, ymax=Hour+SD), width=.35, size = 1) +
  xlab("Strains") +
  ylab("Doubling time (hours)") +
  scale_color_manual("Lineages", values=c("#CC6666", "gray67", "#003366")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size=30, margin=margin(-17, 1, 1, 1)),
        axis.title.y = element_text(size=30, margin=margin(1, 7, 1, 1)),
        axis.text.x = element_text(size=30, angle=-45, hjust=0.1, vjust=0.5),
        axis.text.y = element_text(size=30),
        legend.title = element_text(size=30),
        legend.text = element_text(size=30),
        axis.line = element_line(colour="black")) +
  scale_y_continuous(limits = c(0, 58))
dev.off()

# t-test for doubling times between lineages
t.test(Hour~Lineage, data = double[1:9,]) # lin1vslin2
t.test(Hour~Lineage, data = double[5:15,]) # lin2vslin4
t.test(Hour~Lineage, data = double[c(1:4,10:15),]) # lin1vslin4

# include SD in t-test calculation
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE){
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
}

t.test2(mean(double$Hour[1:4]), mean(double$Hour[5:9]), sd(double$SD[1:4]), sd(double$SD[5:9]), 4, 5)
t.test2(mean(double$Hour[1:4]), mean(double$Hour[10:15]), sd(double$SD[1:4]), sd(double$SD[10:15]), 4, 5)
t.test2(mean(double$Hour[5:9]), mean(double$Hour[10:15]), sd(double$SD[5:9]), sd(double$SD[10:15]), 5, 5)


##### annotate SNP objects #####

# rename annotated SNP files
s345 <- MTB_GT_345_pe.bam.sorted.bam.bcf.variants.vcf.final_variants.ann.csv[1:8]
s119 <- MTB_V119BJ_pe.bam.sorted.bam.bcf.variants.vcf.final_variants.ann.csv[1:8]
s212 <- MTB_V212BJ_pe.bam.sorted.bam.bcf.variants.vcf.final_variants.ann.csv[1:8]
s173 <- MTB_V173IO_pe.bam.sorted.bam.bcf.variants.vcf.final_variants.ann.csv[1:8]
s232 <- MTB_V232IO_pe.bam.sorted.bam.bcf.variants.vcf.final_variants.ann.csv[1:8]
s281 <- MTB_GT_281_pe.sorted.bam.bcf.variants.vcf.final_variants.ann.csv[1:8]
s293 <- MTB_V239EA_pe.bam.sorted.bam.bcf.variants.vcf.final_variants.ann.csv[1:8]
s318 <- MTB_V318EA_pe.bam.sorted.bam.bcf.variants.vcf.final_variants.ann.csv[1:8]
s346 <- MTB_V346IO_pe.bam.sorted.bam.bcf.variants.vcf.final_variants.ann.csv[1:8]
s372 <- MTB_V372IO_pe.bam.sorted.bam.bcf.variants.vcf.final_variants.ann.csv[1:8]
s374 <- MTB_V374BJ_pe.bam.sorted.bam.bcf.variants.vcf.final_variants.ann.csv[1:8]
s440 <- MTB_V440EA_pe.bam.sorted.bam.bcf.variants.vcf.final_variants.ann.csv[1:8]
s639 <- MTB_V639_pe.bam.sorted.bam.bcf.variants.vcf.final_variants.ann.csv[1:8]
s649 <- MTB_GT_649_pe.bam.sorted.bam.bcf.variants.vcf.final_variants.ann.csv[1:8]

# rename strain name
s345$Strain <- "345"
s119$Strain <- "119"
s212$Strain <- "212"
s173$Strain <- "173"
s232$Strain <- "232"
s281$Strain <- "281"
s293$Strain <- "293"
s318$Strain <- "318"
s346$Strain <- "346"
s372$Strain <- "372"
s374$Strain <- "374"
s440$Strain <- "440"
s639$Strain <- "639"
s649$Strain <- "649"

# combine all strains
allsnp_strains <- rbind(s281, s345, s119, s212, s173, s232, s293, s318, s346, s372, s374, s639, s649, s440)

# filter for high quality score
allsnp_strains <- subset(allsnp_strains, allsnp_strains$QUAL > 50)

# clean up SNP annotation based on first annotation
snp_dataset <- allsnp_strains
modsplit=function(x)strsplit(x, "\\,A|\\,G|\\,T|\\,C")
rowmodsplit <- function(x){
  lapply(as.character(x), modsplit)
}
snp_dataset$INFO <- lapply(snp_dataset$INFO, rowmodsplit)
# success[[1]][[1]][[1]][[1]]
# remove layers of the list
snp_dataset$INFO <- sapply(snp_dataset$INFO, `[[`, 1)
snp_dataset$INFO <- sapply(snp_dataset$INFO, `[`, 1)
# split by |
modsplit2=function(x)strsplit(x, "\\|")
rowmodsplit2 <- function(x){
  lapply(as.character(x), modsplit2)
}
snp_dataset$INFO <- lapply(snp_dataset$INFO, rowmodsplit2)
snp_dataset$INFO <- sapply(snp_dataset$INFO, `[[`, 1)
# select the info you want
snp_dataset2 <- snp_dataset
snp_dataset2$INFO <- sapply(snp_dataset2$INFO, `[`, 11)
# snp_dataset2$INFO <- sapply(snp_dataset2$INFO, `[`, 4) # gene name
# select the info you want
snp_dataset_snp <- snp_dataset
snp_dataset_snp$INFO <- sapply(snp_dataset_snp$INFO, `[`, 2)

# keep only single nucleotide variants
snp_dataset_snp <- subset(snp_dataset_snp, snp_dataset_snp$REF %in% c("A", "G", "T", "C"))
snp_dataset_snp <- subset(snp_dataset_snp, snp_dataset_snp$ALT %in% c("A", "G", "T", "C"))
# synonymous and non-synonymous as individual objects
success_syn <- subset(snp_dataset_snp, snp_dataset_snp$INFO == "synonymous_variant")
success_nonsyn <- subset(snp_dataset_snp, snp_dataset_snp$INFO != "synonymous_variant")

# keeping gene info for the non-synonymous snps
snp_dataset_snp2 <- snp_dataset
snp_dataset_snp2$INFO <- sapply(snp_dataset_snp2$INFO, `[`, 4)
snp_dataset_snp2 <- subset(snp_dataset_snp2, snp_dataset_snp2$REF %in% c("A", "G", "T", "C"))
snp_dataset_snp2 <- subset(snp_dataset_snp2, snp_dataset_snp2$ALT %in% c("A", "G", "T", "C"))
snp_dataset_snp2 <- subset(snp_dataset_snp2, snp_dataset_snp2$POS %in% success_nonsyn$POS)

# synonymous snps
df1 <- cbind(success_syn, dummy(success_syn$Strain, sep = "_"))
df_ordered <- df1[order(df1$POS), ]
colnames(df_ordered) <- gsub("success_syn_", "s", colnames(df_ordered))

mine <- df_ordered %>% distinct(POS, .keep_all = TRUE)
mine <- mine$POS
#'%ni%' <- Negate('%in%')

# select SNP from all variants
df_ordered <- subset(df_ordered, df_ordered$REF %in% c("A", "G", "T", "C"))
df_ordered <- subset(df_ordered, df_ordered$ALT %in% c("A", "G", "T", "C"))

empty_table <- list()
for (position in mine){ #for snp position
  zero_one <- list()
  my_table <- df_ordered[df_ordered$POS == position, ] #grab all rows with this position
  lineages_present <- my_table$Strain
  lineages_present <- paste("s", lineages_present, sep="") #grab all the lineages where this snp is present (rename to fit with names)
  lineages_absent <- setdiff(x = colnames(df_ordered)[10:23], y = lineages_present) #grab all lineages where this snp is absent
  zeros <- rep(0, length(lineages_absent)) #for every lineage where its absent add a zero
  ones <- rep(1, length(lineages_present)) #for every lineage where its present add a one
  zeros_ones <- c(zeros, ones)
  names(zeros_ones) <- c(lineages_absent, lineages_present)
  zeros_ones <- zeros_ones[colnames(df_ordered)[10:23]] #reorders the columns to the order of linesgaesin the original table
  empty_table <- rbind(empty_table, zeros_ones) #row = snp, add each row one at a time
}
rownames(empty_table) <- mine

# nonsynonymous snps
df1 <- cbind(success_nonsyn, dummy(success_nonsyn$Strain, sep = "_"))
df_ordered <- df1[order(df1$POS), ]
colnames(df_ordered) <- gsub("success_nonsyn_", "s", colnames(df_ordered))

mine <- df_ordered %>% distinct(POS, .keep_all = TRUE)
mine <- mine$POS

# select SNP from all variants
df_ordered <- subset(df_ordered, df_ordered$REF %in% c("A", "G", "T", "C"))
df_ordered <- subset(df_ordered, df_ordered$ALT %in% c("A", "G", "T", "C"))

empty_table2 <- list()
for (position in mine){ #for snp position
  zero_one <- list()
  my_table <- df_ordered[df_ordered$POS == position, ] #grab all rows with this position
  lineages_present <- my_table$Strain
  lineages_present <- paste("s", lineages_present, sep="") #grab all the lineages where this snp is present, note: grabs according to string which must be found in lineages_present and lineages_absent
  lineages_absent <- setdiff(colnames(df_ordered)[10:23], lineages_present) #grab all lineages where this snp is absent
  zeros <- rep(0, length(lineages_absent)) #for every lineage where its absent add a zero
  ones <- rep(1, length(lineages_present)) #for every lineage where its present add a one
  zeros_ones <- c(zeros, ones)
  names(zeros_ones) <- c(lineages_absent, lineages_present)
  zeros_ones <- zeros_ones[colnames(df_ordered)[10:23]] #reorders the columns to the order of linesgaes in the original table
  empty_table2 <- rbind(empty_table2, zeros_ones) #row = snp, add each row one at a time
}
rownames(empty_table2) <- mine

# all snps
df1 <- cbind(snp_dataset_snp, dummy(snp_dataset_snp$Strain, sep = "_"))
df_ordered <- df1[order(df1$POS), ]
colnames(df_ordered) <- gsub("snp_dataset_snp_", "s", colnames(df_ordered))
#new <- df_ordered[df_ordered == "1", ]
#new2 <- split(df_ordered, df_ordered$POS)

mine <- df_ordered %>% distinct(POS, .keep_all = TRUE)
mine <- mine$POS
#'%ni%' <- Negate('%in%')

# select SNP from all variants
df_ordered <- subset(df_ordered, df_ordered$REF %in% c("A", "G", "T", "C"))
df_ordered <- subset(df_ordered, df_ordered$ALT %in% c("A", "G", "T", "C"))

empty_table3 <- list()
for (position in mine){ #for snp position
  zero_one <- list()
  my_table <- df_ordered[df_ordered$POS == position, ] #grab all rows with this position
  lineages_present <- my_table$Strain
  lineages_present <- paste("s", lineages_present, sep="_") #grab all the lineages where this snp is present
  lineages_absent <- setdiff(colnames(df_ordered)[10:23], lineages_present) #grab all lineages where this snp is absent
  zeros <- rep(0, length(lineages_absent)) #for every lineage where its absent add a zero
  ones <- rep(1, length(lineages_present)) #for every lineage where its present add a one
  zeros_ones <- c(zeros, ones)
  names(zeros_ones) <- c(lineages_absent, lineages_present)
  zeros_ones <- zeros_ones[colnames(df_ordered)[10:23]] #reorders the columns to the order of linesgaesin the original table
  empty_table3 <- rbind(empty_table3, zeros_ones) #row = snp, add each row one at a time
}
rownames(empty_table3) <- mine

##### SNP dataframes #####

# save the files and re-read in to use
# # write.csv(empty_table3, "SNPs_03012020.csv")
# # write.csv(empty_table2, "nonsynonymous_snps03012020.csv")
# # write.csv(empty_table, "synonymous_snps07102020.csv")
# # write.csv(empty_table3, "snps_all.csv")

# read in files
# snps <- read.csv("SNPs_03012020.csv", row=1)
snps_nonsynonymous <- read.csv("nonsynonymous_snps03012020.csv", row=1)
snps_synonymous <- read.csv("synonymous_snps07102020.csv", row = 1)
snps_all2 <- read.csv("snps_all.csv", row=1)

snps_all2 <- as.data.frame(snps_all2)
snps_nonsynonymous <- as.data.frame(snps_nonsynonymous)

# other variant types
# # write.csv(empty_table2, "nonsynonymous_variants.csv")
# variants_nonsynonymous <- read.csv("nonsynonymous_variants.csv", row=1)
# variants_nonsynonymous <- as.data.frame(variants_nonsynonymous)

snps_filtered <- snps_nonsynonymous[ rowSums(snps_nonsynonymous) > 1, ]
snps_filtered <- snps_filtered[ rowSums(snps_filtered) < 13, ]
rownames(snps_filtered) <- make.names(rownames(snps_filtered))

##### SNP analysis #####

# number of SNPs by lineage
snps_synonymous <- as.data.frame(snps_synonymous)
sum2 <- as.data.frame(colSums(snps_synonymous))
sum2$Strains <- rownames(sum2)
sum2$sum <- sum2$`colSums(snps_synonymous)`
sum2$sum <- as.numeric(as.character(sum2$sum))
sum2$lin <- c("lin2", "lin4", "lin2", "lin1", "lin1", "lin4", "lin4", "lin2", "lin1", "lin1", "lin2", "lin4", "lin4", "lin2")
# sum3$lin_colours <- string.to.colors(c("lin1", "lin2", "lin4", "lin2", "lin1", "lin4", "lin4", "lin1", "lin1", "lin1", "lin4", "lin1"))

# stacked barplot showing synonymous vs non-synonymous SNPs
snps_nonsynonymous <- as.data.frame(snps_nonsynonymous)
sum1 <- colSums(snps_nonsynonymous)
sum1 <- as.data.frame(sum1)
sum1$Strains <- rownames(sum1)

sum1$non_syn <- sum1$sum1
sum3 <- sum1[order(match(sum1$Strains,sum2$Strains)),]
sum3$syn <- sum2$sum
sum3$lin <- c("lin2", "lin4", "lin2", "lin1", "lin1", "lin4", "lin4", "lin2", "lin1", "lin1", "lin2", "lin4", "lin4", "lin2")

png("number_nonsynonymoussnps_poster.png", width = 700, height = 575)

ggplot(sum3, aes(x=reorder(Strains, non_syn), y=sum1, fill=lin)) +
  geom_bar(stat="identity") +
  xlab("Strains") +
  ylab("Number of non-synonymous SNPs") +
  # geom_text(aes(label=non_syn), vjust=1.6, color="white", size=3.5) +
  scale_fill_manual("Lineages", values=c("#CC6666", "gray67", "#003366")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size=25, margin=margin(1, 1, 1, 1)),
        axis.title.y = element_text(size=25, margin=margin(1, 7, 1, 1)),
        axis.text.x = element_text(size=21, angle=-45, hjust=0.1, vjust=0.5),
        axis.text.y = element_text(size=21),
        legend.title = element_text(size=25),
        legend.text = element_text(size=22),
        axis.line = element_line(colour="black")) +
  scale_y_continuous(expand = c(0.005, 10), limits = c(0, 1650))

dev.off()

# number of SNPs by lineage
sum01 <- colSums(snps_all2)
sum01 <- as.data.frame(sum01)
sum01$Strains <- rownames(sum01)
sum01$sum01 <- as.numeric(as.character(sum01$sum))
sum_all <- cbind(sum01, sum3)
sum_all <- sum_all[, c(1:2, 5:7)]
colnames(sum_all) <- c("sum_all", "Strains", "sum_non", "sum_syn", "lin")

png("number_allsnps_poster.png", width = 700, height = 575)

sum01$lin <- c("lin2", "lin4", "lin2", "lin1", "lin1", "lin4", "lin4", "lin2", "lin1", "lin1", "lin2", "lin4", "lin4", "lin2")
ggplot(data=sum01, aes(x=reorder(Strains, sum01), y=sum01, fill=lin)) +
  geom_bar(stat="identity") + xlab("Strains") +
  ylab("Number of SNPs") +
  # geom_text(aes(label=sum), vjust=1.6, color="white", size=3.5) +
  scale_fill_manual("Lineages", values=c("#CC6666", "gray67", "#003366")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size=25, margin=margin(1, 1, 1, 1)),
        axis.title.y = element_text(size=25, margin=margin(1, 7, 1, 1)),
        axis.text.x = element_text(size=21, angle=-45, hjust=0.1, vjust=0.5),
        axis.text.y = element_text(size=21),
        legend.title = element_text(size=25),
        legend.text = element_text(size=22),
        axis.line = element_line(colour="black")) +
  scale_y_continuous(expand = c(0.005, 10), limits = c(0, 2250))

dev.off()

# cluster by synonymous and nonsynonymous SNPs
snps_dend <- snps_nonsynonymous
snps_dend <- as.data.frame(t(snps_dend))
# rownames(snps_dend) <- gsub("X", "", rownames(snps_dend))
d <- dist(snps_dend, method = "euclidean") # distance matrix
fit_hc <- hclust(d, method="ward")

plot(fit_hc, hang = -1, 
     main = "Cluster dendrogram", sub = NULL,
     xlab = "Mtb Strain", ylab = "Height", cex=0.85,
     horiz = TRUE) # display dendogram

fit_hc <- as.dendrogram(fit_hc)

colors = c("gray52", "#003366", "#CC5555")
# colors = c("gray47", "#5566AA", "#F17B4B")
fit_hc <- hclust(d, method="ward")
clus3 = cutree(fit_hc, 3)

##### SNP block analysis #####
# SNP blocks to representation score

# select blocks in all the combinations
snp_blocks <- unique(snps_nonsynonymous)

# get gene info
snp_dataset_snp2 <- snp_dataset
# select only SNP variants
snp_dataset_snp2$INFO <- sapply(snp_dataset_snp2$INFO, `[`, 4)
snp_dataset_snp2 <- subset(snp_dataset_snp2, snp_dataset_snp2$REF %in% c("A", "G", "T", "C"))
snp_dataset_snp2 <- subset(snp_dataset_snp2, snp_dataset_snp2$ALT %in% c("A", "G", "T", "C"))

# get gene info
snp_dataset_nonsyn <- snp_dataset
# select only SNP variants
snp_dataset_nonsyn$INFO1 <- sapply(snp_dataset_nonsyn$INFO, `[`, 2)
snp_dataset_nonsyn$INFO2 <- sapply(snp_dataset_nonsyn$INFO, `[`, 3)
snp_dataset_nonsyn$INFO3 <- sapply(snp_dataset_nonsyn$INFO, `[`, 4)
snp_dataset_nonsyn$INFO4 <- sapply(snp_dataset_nonsyn$INFO, `[`, 11)
snp_dataset_nonsyn <- subset(snp_dataset_nonsyn, snp_dataset_nonsyn$REF %in% c("A", "G", "T", "C"))
snp_dataset_nonsyn <- subset(snp_dataset_nonsyn, snp_dataset_nonsyn$ALT %in% c("A", "G", "T", "C"))
snp_dataset_nonsyn <- snp_dataset_nonsyn[!snp_dataset_nonsyn$INFO1 == "synonymous_variant", -7]

# select non-synonymous SNPs only
snp_dataset_snp2 <- subset(snp_dataset_snp2, snp_dataset_snp2$POS %in% snp_dataset_nonsyn$POS)

lppa <- snp_dataset_snp2[snp_dataset_snp2$INFO == "lppA", ]

find_gene <- function(x){
  y <- snp_dataset_snp2[snp_dataset_snp2$INFO == x, ]
  y2 <- length(unique(y$Strain))
  y3 <- 2^y2
  df <- data.frame(x, y2, y3)
  colnames(df) <- c("gene", "n_strains", "n_blocks")
  print(df)
}

number_of_strains <- lapply(unique(snp_dataset_snp2$INFO), find_gene)

n_strains <- tibble(number_of_strains) # generates a tibble

n_strains <- n_strains %>% unnest() # unnests the list to create a dataframe

##### functional categories and snp blocks calculated #####

# load functional categories
func_ann <- read.delim("Biocyc_functionalcat.txt", header=TRUE)
# clean up category annotation
# split by backspaces
categories <- func_ann
modsplit=function(x)strsplit(x, "//")
rowmodsplit <- function(x){
  lapply(as.character(x), modsplit)
}
categories$cat <- lapply(categories$Pathways.of.gene, rowmodsplit)

categories$cat <- sapply(categories$cat, `[[`, 1)

# select blocks in all the combinations
snp_blocks <- unique(snps_nonsynonymous)

###-- get nsSNP blocks #####
# get gene info
snp_dataset_snp2 <- snp_dataset
# select only SNP variants
snp_dataset_snp2$INFO <- sapply(snp_dataset_snp2$INFO, `[`, 4)
snp_dataset_snp2 <- subset(snp_dataset_snp2, snp_dataset_snp2$REF %in% c("A", "G", "T", "C"))
snp_dataset_snp2 <- subset(snp_dataset_snp2, snp_dataset_snp2$ALT %in% c("A", "G", "T", "C"))

# get gene info
snp_dataset_nonsyn <- snp_dataset
# select only SNP variants
snp_dataset_nonsyn$INFO1 <- sapply(snp_dataset_nonsyn$INFO, `[`, 2)
snp_dataset_nonsyn$INFO2 <- sapply(snp_dataset_nonsyn$INFO, `[`, 3)
snp_dataset_nonsyn$INFO3 <- sapply(snp_dataset_nonsyn$INFO, `[`, 4)
snp_dataset_nonsyn$INFO4 <- sapply(snp_dataset_nonsyn$INFO, `[`, 11)
snp_dataset_nonsyn <- subset(snp_dataset_nonsyn, snp_dataset_nonsyn$REF %in% c("A", "G", "T", "C"))
snp_dataset_nonsyn <- subset(snp_dataset_nonsyn, snp_dataset_nonsyn$ALT %in% c("A", "G", "T", "C"))
snp_dataset_nonsyn <- snp_dataset_nonsyn[!snp_dataset_nonsyn$INFO1 == "synonymous_variant", -7]

# select non-synonymous SNPs only
snp_dataset_snp2 <- subset(snp_dataset_snp2, snp_dataset_snp2$POS %in% snp_dataset_nonsyn$POS)


###-- filter by genes in at least 1 snp block with more than 1 snp in #####

cor_snp_block <- cor(t(snps_nonsynonymous))

# select SNP names in blocks
block_names <- rownames(snp_blocks) # gets snp block names
# block_names <- "6112"
listofnames <- list()
for (rowname in block_names){ # cycle through each snp block name in block_names
  block1 <- data.frame(cor_snp_block[rownames(cor_snp_block) == rowname, ]) # subsets correlation matrix by snp block name
  colnames(block1) <- rowname # renames the column
  block1 <- subset(block1, block1[,1] == 1) # selects all perfectly correlated snps with that specific snp block
  block1[, 1] <- rownames(block1) # replaces correlation value with the names of snps which are correlated
  listofnames <- append(listofnames, list(block1)) # appends the list with each dataframe corresponding to each snp block
}

listofnames_t <- tibble(listofnames) # generates a tibble

listofnames_t <- listofnames_t %>% unnest() # unnests the list to create a dataframe

snps_summed <- lapply(listofnames_t, table)
snps_summed <- lapply(snps_summed, as.data.frame)
snps_summed <- do.call(rbind, snps_summed) # make into data frame
snps_summed$snp_block <- rownames(snps_summed) # add snp block names back
snps_summed$snp_block <- gsub("\\..*", "", snps_summed$snp_block) # fix rownames
colnames(genes_summed) <- c("Gene", "Freq", "Snp_block")

snps_wide <- tibble(snps_summed)
snps_wide <- pivot_wider(snps_summed, names_from = snp_block, values_from = Freq)
snps_wide <- snps_wide[, 2:117]

snps_wide <- snps_wide %>%
  replace(is.na(.), 0) %>%
  summarise_all(funs(sum))

#snps_wide is the number of snps in each block

###-- for number of blocks gene is actually in #####
snp_dataset_unique2 <- distinct(snp_dataset_snp2, snp_dataset_snp2$POS, .keep_all=TRUE)
# subsetting the snp_dataset_unique2 data frame to select the gene names corresponding to the snps in each snp block
block_positions <- list()
for (snp_names in listofnames_t){
  gene_blocks <- snp_dataset_unique2[snp_dataset_unique2$POS %in% snp_names, ]
  gene_blocks2 <- as.data.frame(gene_blocks$INFO)
  colnames(gene_blocks2) <- gene_blocks$POS[1]
  block_positions <- append(block_positions, list(gene_blocks2))
}

# convert to data frame for further steps
blocks_of_genes <- tibble(block_positions) # generates a tibble
blocks_of_genes <- blocks_of_genes %>% unnest() # unnests the list to create a dataframe
blocks_of_genes2 <- as.data.frame(blocks_of_genes)

# blocks_of_genes2 is the name of the snps in each block

# get number of blocks the genes are in
n_blocks_gene <- function(x){
  gene_in_block <- blocks_of_genes %>%
    filter_all(any_vars(str_ends(., x)))
  gene_in_block <- gene_in_block[, colSums(is.na(gene_in_block)) != nrow(gene_in_block)]
  in_blocks <- length(colnames(gene_in_block)) # add the blocks
  df <- data.frame(x, in_blocks)
  colnames(df) <- c("gene", "blocks_in")
  print(df)
}

# get the names of the blocks
n_blocks_gene_names <- function(x){
  gene_in_block <- blocks_of_genes %>%
    filter_all(any_vars(str_ends(., x)))
  gene_in_block <- gene_in_block[, colSums(is.na(gene_in_block)) != nrow(gene_in_block)]
  in_blocks <- colnames(gene_in_block) # names of the blocks
  in_blocks2 <- length(colnames(gene_in_block)) # add the blocks
  if(in_blocks2 >= 1){
    df <- data.frame(x, in_blocks, in_blocks2)
    colnames(df) <- c("gene", "blocks_name", "blocks_in")
    print(df)
  }
}

# get total number of SNPs associated with that gene in the whole dataset
all_genes <- pivot_wider(snp_dataset_unique2, names_from = INFO, values_from = POS)

# all_genes is the name of the snp block for each gene for each strain

all_genes_summed <- lapply(all_genes[, 9:1966], table)
all_genes_summed <- lapply(all_genes_summed, sum)
all_genes2 <- as.data.frame(all_genes_summed)
all_genes2 <- as.data.frame(t(all_genes2[1, ]))
colnames(all_genes2) <- "Freq"
all_genes2$Genes <- rownames(all_genes2)

# check
# all_genes2[all_genes2$Genes == "dnaA",]

n_blocks <- lapply(colnames(all_genes[, 9:1966]), n_blocks_gene)
n_blocks_names <- lapply(colnames(all_genes[, 9:1966]), n_blocks_gene_names)

n_blocks_fil <- tibble(n_blocks) # generates a tibble

n_blocks_fil <- n_blocks_fil %>% unnest() # unnests the list to create a dataframe

n_blocks_names <- tibble(n_blocks_names) # generates a tibble
n_blocks_names <- n_blocks_names %>% unnest() # unnests the list to create a dataframe
n_blocks_names <- as.data.frame(n_blocks_names)

# n_blocks <- lapply(all_genes[, 9:2456], n_blocks_gene)

# n_blocks <- tibble(n_blocks) # generates a tibble

# n_blocks <- n_blocks %>% unnest() # unnests the list to create a dataframe
# n_blocks$gene <- n_strains$gene

final_rep <- cbind(n_strains, n_blocks_fil)

final_rep$n_b <- (16384/116) # number of potential blocks (2^14) / actual blocks (116)
final_rep$rep_score <- final_rep$n_b/(final_rep$n_blocks/final_rep$blocks_in)
summary(final_rep) # highest is 141.241, median is 0

final_rep2 <- final_rep
final_rep2$b_over_n <- (final_rep2$blocks_in/116) # number of ones the gene is in / actual blocks (116)
final_rep2$second_score <- final_rep2$b_over_n/(final_rep2$n_blocks/16384)
final_rep2$second_score[final_rep2$second_score == "Inf"] <- 0
summary(final_rep2) # highest is 141.241, median is 0

write.csv(final_rep2, "second_score_thesis.csv")

# filter data
final_rep_top <- final_rep[final_rep$rep_score > 1, ] # remove genes with low counts as this could be by chance
final_rep_top2 <- final_rep_top[final_rep_top$rep_score > 9, ]

rep_score <- final_rep_top[, c(1, 3, 5, 7)]

# add block names
index <- match(rep_score$gene, n_blocks_names$gene)
rep_score$block1_name <- n_blocks_names$blocks_name[index]
duplicated_names <- n_blocks_names[duplicated(as.character(n_blocks_names$gene)), ]
index <- match(rep_score$gene, duplicated_names$gene)
rep_score$block2_name <- duplicated_names$blocks_name[index]
duplicated_names <- duplicated_names[duplicated(as.character(duplicated_names$gene)), ]
index <- match(rep_score$gene, duplicated_names$gene)
rep_score$block3_name <- duplicated_names$blocks_name[index]
duplicated_names <- duplicated_names[duplicated(as.character(duplicated_names$gene)), ]
index <- match(rep_score$gene, duplicated_names$gene)
rep_score$block4_name <- duplicated_names$blocks_name[index]
duplicated_names <- duplicated_names[duplicated(as.character(duplicated_names$gene)), ]
index <- match(rep_score$gene, duplicated_names$gene)
rep_score$block5_name <- duplicated_names$blocks_name[index]
duplicated_names <- duplicated_names[duplicated(as.character(duplicated_names$gene)), ]
index <- match(rep_score$gene, duplicated_names$gene)
rep_score$block6_name <- duplicated_names$blocks_name[index]
duplicated_names <- duplicated_names[duplicated(as.character(duplicated_names$gene)), ]
index <- match(rep_score$gene, duplicated_names$gene)
rep_score$block7_name <- duplicated_names$blocks_name[index]
duplicated_names <- duplicated_names[duplicated(as.character(duplicated_names$gene)), ]
index <- match(rep_score$gene, duplicated_names$gene)
rep_score$block8_name <- duplicated_names$blocks_name[index]
duplicated_names <- duplicated_names[duplicated(as.character(duplicated_names$gene)), ]
index <- match(rep_score$gene, duplicated_names$gene)
rep_score$block9_name <- duplicated_names$blocks_name[index]

write.csv(rep_score[order(rep_score$rep_score, decreasing = TRUE), ], "rep_score_thesis.csv")

index <- match(final_rep$gene, n_blocks_names$gene)
final_rep$block1_name <- n_blocks_names$blocks_name[index]
duplicated_names <- n_blocks_names[duplicated(as.character(n_blocks_names$gene)), ]
index <- match(final_rep$gene, duplicated_names$gene)
final_rep$block2_name <- duplicated_names$blocks_name[index]
duplicated_names <- duplicated_names[duplicated(as.character(duplicated_names$gene)), ]
index <- match(final_rep$gene, duplicated_names$gene)
final_rep$block3_name <- duplicated_names$blocks_name[index]
duplicated_names <- duplicated_names[duplicated(as.character(duplicated_names$gene)), ]
index <- match(final_rep$gene, duplicated_names$gene)
final_rep$block4_name <- duplicated_names$blocks_name[index]
duplicated_names <- duplicated_names[duplicated(as.character(duplicated_names$gene)), ]
index <- match(final_rep$gene, duplicated_names$gene)
final_rep$block5_name <- duplicated_names$blocks_name[index]
duplicated_names <- duplicated_names[duplicated(as.character(duplicated_names$gene)), ]
index <- match(final_rep$gene, duplicated_names$gene)
final_rep$block6_name <- duplicated_names$blocks_name[index]
duplicated_names <- duplicated_names[duplicated(as.character(duplicated_names$gene)), ]
index <- match(final_rep$gene, duplicated_names$gene)
final_rep$block7_name <- duplicated_names$blocks_name[index]
duplicated_names <- duplicated_names[duplicated(as.character(duplicated_names$gene)), ]
index <- match(final_rep$gene, duplicated_names$gene)
final_rep$block8_name <- duplicated_names$blocks_name[index]
duplicated_names <- duplicated_names[duplicated(as.character(duplicated_names$gene)), ]
index <- match(final_rep$gene, duplicated_names$gene)
final_rep$block9_name <- duplicated_names$blocks_name[index]

snp_dataset_unique2 <- distinct(snp_dataset_snp2, snp_dataset_snp2$POS, .keep_all=TRUE)

# subsetting the snp_dataset2 data frame to select the gene names corresponding to the snps in each snp block
block_positions <- list()
for (snp_names in listofnames_t){
  gene_blocks <- snp_dataset_unique2[snp_dataset_unique2$POS %in% snp_names, ]
  gene_blocks2 <- as.data.frame(gene_blocks$INFO)
  colnames(gene_blocks2) <- gene_blocks$POS[1]
  block_positions <- append(block_positions, list(gene_blocks2))
}

# convert to data frame for further steps
blocks_of_genes <- tibble(block_positions) # generates a tibble
blocks_of_genes <- blocks_of_genes %>% unnest() # unnests the list to create a dataframe

# sum the genes
genes_summed <- lapply(blocks_of_genes, table)
genes_summed <- lapply(genes_summed, as.data.frame)
genes_summed <- do.call(rbind, genes_summed) # make into data frame
genes_summed$snp_block <- rownames(genes_summed) # add snp block names back
genes_summed$snp_block <- gsub("\\..*", "", genes_summed$snp_block) # fix rownames
colnames(genes_summed) <- c("Gene", "Freq", "Snp_block")

genes_wide <- tibble(genes_summed)
genes_wide <- pivot_wider(genes_summed, names_from = Snp_block, values_from = Freq)
genes_wide <- genes_wide[, 2:117]

genes_wide <- genes_wide %>%
  replace(is.na(.), 0) %>%
  summarise_all(funs(sum))

numberof_snps <- as.data.frame(t(genes_wide))
colnames(numberof_snps) <- "Number_snps"
numberof_snps$block <- rownames(numberof_snps)
numberof_snps_morethan1 <- numberof_snps[numberof_snps$Number_snps > 1, ]

# add gene annotations
idx <- match(rep_score$gene, categories$Name)
rep_score$protein_bio <- categories$Product[idx]
idx <- match(rep_score$gene, categories$Name)
rep_score$protein_bio <- categories$Pathways.of.gene[idx]
idx <- match(rep_score$gene, categories$Name)
rep_score$pathway <- categories$GO.terms..biological.process.[idx]

# get annotation
gene_names <- read.delim("Mycobacterium_smegmatis_MC2-155_txt_v3.txt")
gene_names_protein <- read.delim("ProteinTable1026_300490.txt")
idx <- match(rep_score$gene, gene_names_protein$Locus)
rep_score$protein <- gene_names_protein$Protein.name[idx]

rep_score_filtered <- rep_score[which(rep_score$block1_name %in% rownames(numberof_snps_morethan1) | rep_score$block2_name %in% rownames(numberof_snps_morethan1) | rep_score$block3_name %in% rownames(numberof_snps_morethan1) | rep_score$block4_name %in% rownames(numberof_snps_morethan1) | rep_score$block5_name %in% rownames(numberof_snps_morethan1) | rep_score$block6_name %in% rownames(numberof_snps_morethan1) | rep_score$block7_name %in% rownames(numberof_snps_morethan1) | rep_score$block8_name %in% rownames(numberof_snps_morethan1) | rep_score$block9_name %in% rownames(numberof_snps_morethan1)), ]
write.csv(rep_score_filtered[order(rep_score_filtered$rep_score, decreasing = TRUE), ], "rep_score_filtered_thesis.csv")

# select only blocks present in more than 1 strain
morethan1_block <- snp_blocks[rownames(snp_blocks) %in% rownames(numberof_snps_morethan1), ]
morethan1_block[length(rownames(morethan1_block))+1, ] <- c("lin2", "lin4", "lin2", "lin1", "lin1", "lin4", "lin4", "lin2", "lin1", "lin1", "lin2", "lin4", "lin4", "lin2")
morethan1_block[length(rownames(morethan1_block))+1, ] <- c("Beijing", "H1", "Beijing", "EA14_VNM", "ZERO", "T1", "T1", "Beijing-Like", "EA14_NVM", "EA14_NVM", "Beijing", "H3", "H3", "Beijing-Like")
morethan1_block <- as.data.frame(t(morethan1_block))
morethan1_block <- morethan1_block[c(4, 9, 10, 5, 1, 3, 11, 8, 14, 2, 6, 7, 12, 13), ]

# get information on those blocks
morethan1_block2 <- lapply(morethan1_block[1:26],as.numeric)
df1 <- sapply(morethan1_block2, sum)
summary(df1)
sum(df1)/26

# this is the same as previous analysis
write.csv(morethan1_block, "block_list_morethan1snpperblock_thesis.csv")

rs57298 <- rbind(subset(rep_score, block1_name == 57298), subset(rep_score, block2_name == 57298))
rs79480 <- rbind(subset(rep_score, block1_name == 79480), subset(rep_score, block2_name == 79480))
rs2867458 <- rbind(subset(rep_score, block8_name == 2867458))
rs1722228 <- rbind(subset(rep_score, block1_name == 1722228), subset(rep_score, block2_name == 1722228))
rs467526 <- rbind(subset(rep_score, block1_name == 467526))
rs1849 <- rbind(subset(rep_score, block1_name == 1849))
rs6112 <- rbind(subset(rep_score, block1_name == 6112))
rs6112_filtered <- subset(rs6112, rs6112$rep_score > 1)
rs15117 <- rbind(subset(rep_score, block2_name == 15117))
rs336560 <- rbind(subset(rep_score, block2_name == 336560), subset(rep_score, block3_name == 336560))
rs25137 <- rbind(subset(rep_score, block1_name == 25137), subset(rep_score, block2_name == 25137), subset(rep_score, block3_name == 25137))

idx <- match(final_rep$gene, categories$Name)
final_rep$protein_bio <- categories$Product[idx]

rs57298 <- rbind(subset(final_rep, block1_name == 57298), subset(final_rep, block2_name == 57298))
rs79480 <- rbind(subset(final_rep, block1_name == 79480), subset(final_rep, block2_name == 79480))
rs1849 <- rbind(subset(final_rep, block1_name == 1849))
rs6112 <- rbind(subset(final_rep, block1_name == 6112))
rs15117 <- rbind(subset(final_rep, block2_name == 15117))
rs25137 <- rbind(subset(final_rep, block1_name == 25137), subset(final_rep, block2_name == 25137), subset(final_rep, block3_name == 25137))

write.csv(rs57298[order(rs57298$rep_score, decreasing = TRUE), ], "rep_score_rs57298block_thesis.csv")
write.csv(rs79480[order(rs79480$rep_score, decreasing = TRUE), ], "rep_score_rs79480block2_thesis.csv")
write.csv(rs2867458[order(rs2867458$rep_score, decreasing = TRUE), ], "rep_score_rs2867458block_thesis.csv")
write.csv(rs1722228[order(rs1722228$rep_score, decreasing = TRUE), ], "rep_score_rs1722228block_thesis.csv")
write.csv(rs467526[order(rs467526$rep_score, decreasing = TRUE), ], "rep_score_rs467526block_thesis.csv")
write.csv(rs1849[order(rs1849$rep_score, decreasing = TRUE), ], "rep_score_rs1849block2_thesis.csv")
write.csv(rs1849[-grep("Rv", rs1849$gene), ], "rep_score_rs1849_norv_genes_thesis.csv")
write.csv(rs6112[order(rs6112$rep_score, decreasing = TRUE), ], "rep_score_rs6112block_thesis.csv")
write.csv(rs15117[order(rs15117$rep_score, decreasing = TRUE), ], "rep_score_rs15117block_thesis.csv")
write.csv(rs25137[order(rs25137$rep_score, decreasing = TRUE), ], "rep_score_rs25137block_thesis.csv")

# look at everything that is correlated with ergothioneine
erg1 <- rbind(subset(final_rep, block3_name == 2122395), subset(final_rep, block1_name == 467526), subset(final_rep, block1_name == 1190093))
erg2 <- rbind(subset(final_rep, block1_name == 2509140), subset(final_rep, block2_name == 3462135), subset(final_rep, block2_name == 21795), subset(final_rep, block1_name == 2939657))
erg3 <- rbind(subset(final_rep, block1_name == 190646), subset(final_rep, block2_name == 368092),
              subset(final_rep, block1_name == 103823), subset(final_rep, block1_name == 103836),
              subset(final_rep, block8_name == 3735749), subset(final_rep, block2_name == 3137237),
              subset(final_rep, block1_name == 595501), subset(final_rep, block1_name == 1481185),
              subset(final_rep, block1_name == 77870))
# subset(final_rep, block1_name == 6112))

narl <- snp_dataset2[snp_dataset2$INFO == "narL", ]
glnq <- snp_dataset2[snp_dataset2$INFO == "glnQ", ]
glna4 <- snp_dataset2[snp_dataset2$INFO == "glnA4", ]

erg_genes <- rbind(narl, glnq, glna4)

# 57298 specific to 2/3 beijing strains, 79480 to 3 beijing and 1 beijing-like strains
# 2867458 all lineage 2 plus 318
# 1722228 most of lineage 2
# 467526 all but lineage 1
geneofinterest <- snp_dataset_snp2[snp_dataset_snp2$INFO == "aceAb", ]
geneofinterest <- snp_dataset_snp2[snp_dataset_snp2$INFO == "serA1", ]
geneofinterest <- snp_dataset_snp2[snp_dataset_snp2$INFO == "PE9", ]
# 57298 contains aceAb which is present in 4 strains 212, 374, 639 and 440, not all beijing strains

# get amino acid info
snp_dataset_amino <- snp_dataset
# select only SNP variants
snp_dataset_amino$INFO <- sapply(snp_dataset_amino$INFO, `[`, 11)
# get amino acid info
snp_dataset_amino <- subset(snp_dataset_amino, snp_dataset_amino$REF %in% c("A", "G", "T", "C"))
snp_dataset_amino <- subset(snp_dataset_amino, snp_dataset_amino$ALT %in% c("A", "G", "T", "C"))
# select non-synonymous SNPs only
snp_dataset_amino <- subset(snp_dataset_amino, snp_dataset_amino$POS %in% success_nonsyn$POS)

# snp_dataset_amino <- snp_dataset_nonsyn
snp_dataset_amino$INFO4 <- snp_dataset_amino$INFO

snp_dataset_amino[snp_dataset_amino$POS %in% geneofinterest$POS, ]
snp_dataset_amino$INFO4 <- gsub("Arg", "R", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("His", "H", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Lys", "K", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Asp", "D", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Glu", "E", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Ser", "S", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Thr", "T", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Asn", "N", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Gln", "Q", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Cys", "C", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Gly", "G", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Pro", "P", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Ala", "A", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Val", "V", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Ile", "I", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Leu", "L", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Met", "M", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Phe", "F", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Tyr", "Y", snp_dataset_amino$INFO4)
snp_dataset_amino$INFO4 <- gsub("Trp", "W", snp_dataset_amino$INFO4)

get_aminoacid <- function(x){
  geneofinterest <- snp_dataset_snp2[snp_dataset_snp2$INFO == x, ]
  df <- snp_dataset_amino[snp_dataset_amino$POS %in% geneofinterest$POS, ]
  # df$gene <- x
  print(df)
}

##### snp block play #####

lin2 <- snp_blocks[, c("s119", "s212", "s374", "s345", "s649")]

##### functional categories cleanup #####

func_ann <- read.delim("Biocyc_functionalcat.txt", header=TRUE)
# clean up category annotation
# split by backspaces
categories <- func_ann
modsplit=function(x)strsplit(x, "//")
rowmodsplit <- function(x){
  lapply(as.character(x), modsplit)
}
categories$cat <- lapply(categories$Pathways.of.gene, rowmodsplit)

categories$cat <- sapply(categories$cat, `[[`, 1)

# select the info you want
snp_dataset5 <- snp_dataset
snp_dataset5$INFO <- sapply(snp_dataset5$INFO, `[`, 4)

# keep only single nucleotide variants
snp_dataset5 <- subset(snp_dataset5, snp_dataset5$REF %in% c("A", "G", "T", "C"))
snp_dataset5 <- subset(snp_dataset5, snp_dataset5$ALT %in% c("A", "G", "T", "C"))
# synonymous and non-synonymous as individual objects
# snp_dataset5_nonsyn <- subset(snp_dataset5, snp_dataset5$POS %in% success_nonsyn$POS)

# snp_dataset3 <- snp_dataset2
snp_dataset5$new <- categories$cat[match(snp_dataset5$INFO, categories$Name)]
# check func_ann[func_ann$INFO == "Rv3922c", ]

allsnps_functions <- snp_dataset5
allsnps_function <- allsnps_functions[allsnps_functions$POS %in% rownames(snps_all2), ]
snps_nonsynonymous_fun <- allsnps_functions[allsnps_functions$POS %in% rownames(snps_nonsynonymous),] # NUMBER PROBLEM

# split by strain
snpsep <- split(snps_nonsynonymous_fun, snps_nonsynonymous_fun$Strain)

snp119 <- as.data.frame(snpsep[1])
snp119 <- unnest(snp119, cols = c(X119.new))
snp119$snp_cat <- gsub(" ", "", snp119$X119.new) # remove spaces
snp119 <- snp119[!duplicated(snp119[c(1,3,4,9,11)]),] # ensure no duplicates of category and snp
snp119_counts <- data.frame(counts = summary(as.factor(snp119$snp_cat)))
snp119_counts$pathways <- rownames(snp119_counts)

snp173 <- as.data.frame(snpsep[2])
snp173 <- unnest(snp173, cols = c(X173.new))
snp173$snp_cat <- gsub(" ", "", snp173$X173.new) # remove spaces
snp173 <- snp173[!duplicated(snp173[c(1,3,4,9,11)]),] # ensure no duplicates of category and snp
snp173_counts <- data.frame(counts = summary(as.factor(snp173$snp_cat)))
snp173_counts$pathways <- rownames(snp173_counts)

snp212 <- as.data.frame(snpsep[3])
snp212 <- unnest(snp212, cols = c(X212.new))
snp212$snp_cat <- gsub(" ", "", snp212$X212.new) # remove spaces
snp212 <- snp212[!duplicated(snp212[c(1,3,4,9,11)]),] # ensure no duplicates of category and snp
snp212_counts <- data.frame(counts = summary(as.factor(snp212$snp_cat)))
snp212_counts$pathways <- rownames(snp212_counts)

snp232 <- as.data.frame(snpsep[4])
snp232 <- unnest(snp232, cols = c(X232.new))
snp232$snp_cat <- gsub(" ", "", snp232$X232.new) # remove spaces
snp232 <- snp232[!duplicated(snp232[c(1,3,4,9,11)]),] # ensure no duplicates of category and snp
snp232_counts <- data.frame(counts = summary(as.factor(snp232$snp_cat)))
snp232_counts$pathways <- rownames(snp232_counts)

snp281 <- as.data.frame(snpsep[5])
snp281 <- unnest(snp281, cols = c(X281.new))
snp281$snp_cat <- gsub(" ", "", snp281$X281.new) # remove spaces
snp281 <- snp281[!duplicated(snp281[c(1,3,4,9,11)]),] # ensure no duplicates of category and snp
snp281_counts <- data.frame(counts = summary(as.factor(snp281$snp_cat)))
snp281_counts$pathways <- rownames(snp281_counts)

snp293 <- as.data.frame(snpsep[6])
snp293 <- unnest(snp293, cols = c(X293.new))
snp293$snp_cat <- gsub(" ", "", snp293$X293.new) # remove spaces
snp293 <- snp293[!duplicated(snp293[c(1,3,4,9,11)]),] # ensure no duplicates of category and snp
snp293_counts <- data.frame(counts = summary(as.factor(snp293$snp_cat)))
snp293_counts$pathways <- rownames(snp293_counts)

snp318 <- as.data.frame(snpsep[7])
snp318 <- unnest(snp318, cols = c(X318.new))
snp318$snp_cat <- gsub(" ", "", snp318$X318.new) # remove spaces
snp318 <- snp318[!duplicated(snp318[c(1,3,4,9,11)]),] # ensure no duplicates of category and snp
snp318_counts <- data.frame(counts = summary(as.factor(snp318$snp_cat)))
snp318_counts$pathways <- rownames(snp318_counts)

snp345 <- as.data.frame(snpsep[8])
snp345 <- unnest(snp345, cols = c(X345.new))
snp345$snp_cat <- gsub(" ", "", snp345$X345.new) # remove spaces
snp345 <- snp345[!duplicated(snp345[c(1,3,4,9,11)]),] # ensure no duplicates of category and snp
snp345_counts <- data.frame(counts = summary(as.factor(snp345$snp_cat)))
snp345_counts$pathways <- rownames(snp345_counts)

snp346 <- as.data.frame(snpsep[9])
snp346 <- unnest(snp346, cols = c(X346.new))
snp346$snp_cat <- gsub(" ", "", snp346$X346.new) # remove spaces
snp346 <- snp346[!duplicated(snp346[c(1,3,4,9,11)]),] # ensure no duplicates of category and snp
snp346_counts <- data.frame(counts = summary(as.factor(snp346$snp_cat)))
snp346_counts$pathways <- rownames(snp346_counts)

snp372 <- as.data.frame(snpsep[10])
snp372 <- unnest(snp372, cols = c(X372.new))
snp372$snp_cat <- gsub(" ", "", snp372$X372.new) # remove spaces
snp372 <- snp372[!duplicated(snp372[c(1,3,4,9,11)]),] # ensure no duplicates of category and snp
snp372_counts <- data.frame(counts = summary(as.factor(snp372$snp_cat)))
snp372_counts$pathways <- rownames(snp372_counts)

snp374 <- as.data.frame(snpsep[11])
snp374 <- unnest(snp374, cols = c(X374.new))
snp374$snp_cat <- gsub(" ", "", snp374$X374.new) # remove spaces
snp374 <- snp374[!duplicated(snp374[c(1,3,4,9,11)]),] # ensure no duplicates of category and snp
snp374_counts <- data.frame(counts = summary(as.factor(snp374$snp_cat)))
snp374_counts$pathways <- rownames(snp374_counts)

snp440 <- as.data.frame(snpsep[12])
snp440 <- unnest(snp440, cols = c(X440.new))
snp440$snp_cat <- gsub(" ", "", snp440$X440.new) # remove spaces
snp440 <- snp440[!duplicated(snp440[c(1,3,4,9,11)]),] # ensure no duplicates of category and snp
snp440_counts <- data.frame(counts = summary(as.factor(snp440$snp_cat)))
snp440_counts$pathways <- rownames(snp440_counts)

snp639 <- as.data.frame(snpsep[13])
snp639 <- unnest(snp639, cols = c(X639.new))
snp639$snp_cat <- gsub(" ", "", snp639$X639.new) # remove spaces
snp639 <- snp639[!duplicated(snp639[c(1,3,4,9,11)]),] # ensure no duplicates of category and snp
snp639_counts <- data.frame(counts = summary(as.factor(snp639$snp_cat)))
snp639_counts$pathways <- rownames(snp639_counts)

snp649 <- as.data.frame(snpsep[14])
snp649 <- unnest(snp649, cols = c(X649.new))
snp649$snp_cat <- gsub(" ", "", snp649$X649.new) # remove spaces
snp649 <- snp649[!duplicated(snp649[c(1,3,4,9,11)]),] # ensure no duplicates of category and snp
snp649_counts <- data.frame(counts = summary(as.factor(snp649$snp_cat)))
snp649_counts$pathways <- rownames(snp649_counts)

# get nsSNPs that affect a specific pathway
snp119[snp119$snp_cat == "fattyacid&beta;-oxidationI", ]
snp173[snp173$snp_cat == "fattyacid&beta;-oxidationI", ]
snp212[snp212$snp_cat == "fattyacid&beta;-oxidationI", ]
snp232[snp232$snp_cat == "fattyacid&beta;-oxidationI", ]
snp281[snp281$snp_cat == "fattyacid&beta;-oxidationI", ]
snp293[snp293$snp_cat == "fattyacid&beta;-oxidationI", ]
snp318[snp318$snp_cat == "fattyacid&beta;-oxidationI", ]
snp345[snp345$snp_cat == "fattyacid&beta;-oxidationI", ]

# combine similar categories
get_total <- function(strain_name){
  if(sum(strain_name[grep("adenineandadenosinesalvage*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(strain_name[grep("adenineandadenosinesalvage*", rownames(strain_name)), ]$counts), pathways = "adenineandadenosinesalvage") # add counts for everything matching this string
    x <- strain_name[-grep("adenineandadenosinesalvage*", rownames(strain_name)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("cholesteroldegradation*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("cholesteroldegradation*", rownames(x)), ]$counts), pathways = "cholesteroldegradation") # add counts for everything matching this string
    x <- x[-grep("cholesteroldegradation*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("aminoimidazoleribonucleotidebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("aminoimidazoleribonucleotidebiosynthesis*", rownames(x)), ]$counts), pathways = "aminoimidazoleribonucleotidebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("aminoimidazoleribonucleotidebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("adenosinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("adenosinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]$counts), pathways = "adenosinedeoxyribonucleotides_denovo_biosynthesis") # add counts for everything matching this string
    x <- x[-grep("adenosinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("CDP-diacylglycerolbiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("CDP-diacylglycerolbiosynthesis*", rownames(x)), ]$counts), pathways = "CDP-diacylglycerolbiosynthesis") # add counts for everything matching this string
    x <- x[-grep("CDP-diacylglycerolbiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("serineandglycinebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("serineandglycinebiosynthesis*", rownames(x)), ]$counts), pathways = "serineandglycinebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("serineandglycinebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("L-argininebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("L-argininebiosynthesis*", rownames(x)), ]$counts), pathways = "L-argininebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("L-argininebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("isoleucinebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("isoleucinebiosynthesis*", rownames(x)), ]$counts), pathways = "isoleucinebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("isoleucinebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("L-prolinebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("L-prolinebiosynthesis*", rownames(x)), ]$counts), pathways = "L-prolinebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("L-prolinebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("NADHtocytochrome<i>bd</i>oxidaseelectrontransfer*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("NADHtocytochrome<i>bd</i>oxidaseelectrontransfer*", rownames(x)), ]$counts), pathways = "NADHtocytochrome_bd_oxidaseelectrontransfer") # add counts for everything matching this string
    x <- x[-grep("NADHtocytochrome<i>bd</i>oxidaseelectrontransfer*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("pyrimidinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("pyrimidinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]$counts), pathways = "pyrimidinedeoxyribonucleotides_denovo_biosynthesis") # add counts for everything matching this string
    x <- x[-grep("pyrimidinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("pyruvatefermentationtoacetate", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("pyruvatefermentationtoacetate*", rownames(x)), ]$counts), pathways = "pyruvatefermentationtoacetate") # add counts for everything matching this string
    x <- x[-grep("pyruvatefermentationtoacetate*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("trehalosebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("trehalosebiosynthesis*", rownames(x)), ]$counts), pathways = "trehalosebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("trehalosebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("mycolatebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("mycolatebiosynthesis*", rownames(x)), ]$counts), pathways = "mycolatebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("mycolatebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("superpathwayofguanosinenucleotides<i>denovo</i>biosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("superpathwayofguanosinenucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]$counts), pathways = "superpathwayofguanosinenucleotides_denovo_biosynthesis") # add counts for everything matching this string
    x <- x[-grep("superpathwayofguanosinenucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("guanosinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("guanosinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]$counts), pathways = "guanosinedeoxyribonucleotides_denovo_biosynthesis") # add counts for everything matching this string
    x <- x[-grep("guanosinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("superpathwayofadenosinenucleotides<i>denovo</i>biosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("superpathwayofadenosinenucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]$counts), pathways = "superpathwayofadenosinenucleotides_denovo_biosynthesis") # add counts for everything matching this string
    x <- x[-grep("superpathwayofadenosinenucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("phosphatidylglycerolbiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("phosphatidylglycerolbiosynthesis*", rownames(x)), ]$counts), pathways = "phosphatidylglycerolbiosynthesis") # add counts for everything matching this string
    x <- x[-grep("phosphatidylglycerolbiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("pentosephosphatepathway*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("pentosephosphatepathway*", rownames(x)), ]$counts), pathways = "pentosephosphatepathway") # add counts for everything matching this string
    x <- x[-grep("pentosephosphatepathway*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("L-tyrosinebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("L-tyrosinebiosynthesis*", rownames(x)), ]$counts), pathways = "L-tyrosinebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("L-tyrosinebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("guanineandguanosinesalvage*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("guanineandguanosinesalvage*", rownames(x)), ]$counts), pathways = "guanineandguanosinesalvage") # add counts for everything matching this string
    x <- x[-grep("guanineandguanosinesalvage*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("L-methioninebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("L-methioninebiosynthesis*", rownames(x)), ]$counts), pathways = "L-methioninebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("L-methioninebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("L-tryptophanbiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("L-tryptophanbiosynthesis*", rownames(x)), ]$counts), pathways = "L-tryptophanbiosynthesis") # add counts for everything matching this string
    x <- x[-grep("L-tryptophanbiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("trehalosedegradation*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("trehalosedegradation*", rownames(x)), ]$counts), pathways = "trehalosedegradation") # add counts for everything matching this string
    x <- x[-grep("trehalosedegradation*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  print(x)
}

snp119_counts <- get_total(snp119_counts)
snp173_counts <- get_total(snp173_counts)
snp212_counts <- get_total(snp212_counts)
snp232_counts <- get_total(snp232_counts)
snp281_counts <- get_total(snp281_counts)
snp293_counts <- get_total(snp293_counts)
snp318_counts <- get_total(snp318_counts)
snp345_counts <- get_total(snp345_counts)
snp346_counts <- get_total(snp346_counts)
snp372_counts <- get_total(snp372_counts)
snp374_counts <- get_total(snp374_counts)
snp440_counts <- get_total(snp440_counts)
snp639_counts <- get_total(snp639_counts)
snp649_counts <- get_total(snp649_counts)

# add strain and lineage information to each object
snp119_counts$strain <- "119"
snp173_counts$strain <- "173"
snp212_counts$strain <- "212"
snp232_counts$strain <- "232"
snp281_counts$strain <- "281"
snp293_counts$strain <- "293"
snp318_counts$strain <- "318"
snp345_counts$strain <- "345"
snp346_counts$strain <- "346"
snp372_counts$strain <- "372"
snp374_counts$strain <- "374"
snp440_counts$strain <- "440"
snp639_counts$strain <- "639"
snp649_counts$strain <- "649"

snp119_counts$order_strain <- 5
snp173_counts$order_strain <- 13
snp212_counts$order_strain <- 6
snp232_counts$order_strain <- 3
snp281_counts$order_strain <- 1
snp293_counts$order_strain <- 12
snp318_counts$order_strain <- 14
snp345_counts$order_strain <- 8
snp346_counts$order_strain <- 4
snp372_counts$order_strain <- 2
snp374_counts$order_strain <- 7
snp440_counts$order_strain <- 10
snp639_counts$order_strain <- 11
snp649_counts$order_strain <- 9

snp_cats <- rbind(snp119_counts, snp173_counts, snp212_counts, snp232_counts, snp281_counts, snp293_counts, snp318_counts, snp345_counts, snp346_counts, snp372_counts, snp374_counts, snp440_counts, snp639_counts, snp649_counts)
snp_cats$lineage <- c(rep("lin2", length(snp119_counts$counts)), rep("lin4", length(snp173_counts$counts)), rep("lin2", length(snp212_counts$counts)), rep("lin1", length(snp232_counts$counts)), rep("lin1", length(snp281_counts$counts)), rep("lin4", length(snp293_counts$counts)), rep("lin4", length(snp318_counts$counts)), rep("lin2", length(snp345_counts$counts)), rep("lin1", length(snp346_counts$counts)), rep("lin1", length(snp372_counts$counts)), rep("lin2", length(snp374_counts$counts)), rep("lin4", length(snp440_counts$counts)), rep("lin4", length(snp639_counts$counts)), rep("lin2", length(snp649_counts$counts)))
# snp_cats2 <- cbind.fill(snp119_counts$counts, snp173_counts$counts, snp212_counts$counts, snp232_counts$counts, snp281_counts$counts, snp293_counts$counts, snp318_counts$counts, snp345_counts$counts, snp346_counts$counts, snp372_counts$counts, snp374_counts$counts, snp440_counts$counts, snp639_counts$counts, snp649_counts$counts)
# snp_cats$category <- rownames(snp_cats)
# write.csv(snp_cats2, "snp_functional_biocyc.csv")

condition <- c("119", "173", "212", "232", "281", "239", "318", "345", "346", "372", "374", "440", "639", "649")
order_lin <- c("346", "232", "372", "281", "119", "212", "345", "374", "649", "173", "239", "440", "318", "639")

snp_cats <- snp_cats[!snp_cats$counts == "NA's", ]
snp_cats <- snp_cats[!snp_cats$counts == "NA2", ]
snp_others <- snp_cats[snp_cats$pathways == "(Other)", ]
sum(snp_others$counts)
snp_cats <- snp_cats[!snp_cats$pathways == "(Other)", ]
# snp_cats <- snp_cats[snp_cats$counts > 3, ]
snp_cats <- snp_cats[snp_cats$counts > 5, ]

write.csv(summary(factor(snp_cats$pathways)), "functional_edit.csv")

snp_cats_edited <- read.csv("functional_edited_aminoacids.csv")

idx <- match(snp_cats$pathways, snp_cats_edited$X)

snp_cats$pathways_superfamily <- snp_cats_edited$Category[idx]

colourCount = length(unique(snp_cats$pathways_superfamily))
getPalette = colorRampPalette(brewer.pal(40, "Set1"))

# for lineages functional categories averages across strains
ggplot(snp_cats, aes(fill=pathways_superfamily, x = as.character(lineage) , y = as.numeric(counts)), legend=FALSE) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = getPalette(colourCount)) +
  ylab("Categories (% total)") +
  xlab("Lineage") +
  theme_minimal() +
  theme(axis.title.x = element_text(size=25, margin=margin(1, 1, 1, 1)),
        axis.title.y = element_text(size=25, margin=margin(1, 7, 1, 1))) +
  guides(fill=FALSE)

# for functional catergories for individual strains
png("functional_categories_strains_2_thesis.png", height = 2000, width = 2000)
ggplot(snp_cats, aes(fill=pathways_superfamily, x = reorder(as.character(strain), order_strain) , y = as.numeric(counts)), legend=FALSE) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = getPalette(colourCount)) +
  ylab("Categories (% total)") +
  xlab("Lineage") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title.x = element_text(size=60, margin=margin(25, 1, 1, 1)),
        axis.title.y = element_text(size=60, margin=margin(1, 25, 1, 1)),
        axis.text.x = element_text(size=58, angle=-45, hjust=0.1, vjust=0.75, margin=margin(10, 1, 1, 1)),
        axis.text.y = element_text(size=58, margin=margin(1, 10, 1, 1)),
        # legend.title = element_text(size=52),
        # legend.text = element_text(size=50),
        axis.line = element_line(colour="black", size = 1)) + guides(fill=FALSE)
dev.off()

png("functional_categories_strains_legend_thesis.png", height = 2850, width = 2850)
ggplot(snp_cats, aes(fill=pathways_superfamily, x = reorder(as.character(strain), order_strain) , y = as.numeric(counts)), legend=FALSE) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = getPalette(colourCount)) +
  ylab("Categories (% total)") +
  xlab("Lineage") +
  theme_minimal() +
  theme(legend.position = "left",
        axis.title.x = element_text(size=30, margin=margin(1, 1, 1, 1)),
        axis.title.y = element_text(size=30, margin=margin(1, 1, 1, 1)),
        axis.text.x = element_text(size=28, angle=-45, hjust=0.1, vjust=0.75, margin=margin(1, 1, 1, 1)),
        axis.text.y = element_text(size=28, margin=margin(1, 1, 1, 1)),
        legend.title = element_text(size=50),
        legend.text = element_text(size=48),
        axis.line = element_line(colour="black", size = 1),
        legend.key.size = unit(0.7, "cm"),
        legend.direction = "vertical") # + guides(fill=FALSE)
dev.off()

write.csv(snp_cats, "functional_pathways_thesis.csv")

##### dendrogram snp plots #####

# cluster by synonymous and nonsynonymous SNPs
snps_dend <- snps_nonsynonymous
snps_dend <- as.data.frame(t(snps_dend))
# rownames(snps_dend) <- gsub("X", "", rownames(snps_dend))
d <- dist(snps_dend, method = "euclidean") # distance matrix
fit_hc <- hclust(d, method="ward")

plot(fit_hc, hang = -1, 
     main = "Cluster dendrogram", sub = NULL,
     xlab = "Mtb Strain", ylab = "Height", cex=0.85,
     horiz = TRUE) # display dendogram

fit_hc <- as.dendrogram(fit_hc)

colors = c("gray52", "#003366", "#CC5555")
# colors = c("gray47", "#5566AA", "#F17B4B")
fit_hc <- hclust(d, method="ward")
clus3 = cutree(fit_hc, 3)

fit_hc <- hclust(d, method="ward")

png("dendrogram_nonsynonymous_snps_thesis.png", height = 750, width = 900)
plot(as.phylo(fit_hc), cex = 1.5, label.offset = 0.5, tip.color = colors[clus3])
dev.off()

png("dendrogram_nonsynonymous_thesis.png", height = 750, width = 350)
plot(as.phylo(fit_hc), cex = 1.5, label.offset = 0.5, tip.color = colors[clus3])
dev.off()

png("dendrogram_unrooted_nonsynonymous_thesis.png", height = 2750, width = 2750, res = 500)
plot(as.phylo(fit_hc), cex = 1.25, label.offset = 4.5, tip.color = colors[clus3],
     type = "unrooted", no.margin = TRUE)
dev.off()

##### tree with other lineages #####

all_lins <- read.csv("reference_dataset.csv")

# nonsynonymous snps
# df_lins <- data.frame(matrix(ncol = length(summary(factor(success_nonsyn$Strain))), nrow = dim(success_nonsyn)[1]))
# colnames(df_lins) <- colnames(snps_nonsynonymous)
# all_lins2 <- cbind(all_lins, df_lins)

# df1 <- cbind(success_nonsyn, dummy(success_nonsyn$Strain))
# df_ordered <- df1[order(df1$POS), ]
colnames(df_ordered) <- gsub("success_nonsyn_", "s", colnames(df_ordered))

mine <- df_ordered %>% distinct(POS, .keep_all = TRUE)
mine <- mine$POS

# select SNP from all variants
df_ordered <- subset(df_ordered, df_ordered$REF %in% c("A", "G", "T", "C"))
df_ordered <- subset(df_ordered, df_ordered$ALT %in% c("A", "G", "T", "C"))

empty_table2 <- list()
for (position in mine){ #for snp position
  zero_one <- list()
  my_table <- df_ordered[df_ordered$POS == position, ] #grab all rows with this position
  lineages_present <- my_table$Strain
  lineages_present <- paste("s", lineages_present, sep="") #grab all the lineages where this snp is present, note: grabs according to string which must be found in lineages_present and lineages_absent
  lineages_absent <- setdiff(colnames(df_ordered)[10:23], lineages_present) #grab all lineages where this snp is absent
  zeros <- rep(0, length(lineages_absent)) #for every lineage where its absent add a zero
  ones <- rep(1, length(lineages_present)) #for every lineage where its present add a one
  zeros_ones <- c(zeros, ones)
  names(zeros_ones) <- c(lineages_absent, lineages_present)
  zeros_ones <- zeros_ones[colnames(df_ordered)[10:23]] #reorders the columns to the order of linesgaes in the original table
  empty_table2 <- rbind(empty_table2, zeros_ones) #row = snp, add each row one at a time
}
rownames(empty_table2) <- mine


# add emtpy columns to dataset
df_lins <- data.frame(matrix(ncol = length(summary(factor(all_lins$Tuberculist.category))), nrow = dim(all_lins)[1]))
colnames(df_lins) <- unique(all_lins$Tuberculist.category)
all_lins2 <- cbind(all_lins, df_lins)

mine <- all_lins2 %>% distinct(Genome.position, .keep_all = TRUE)
mine <- mine$Genome.position

empty_table_ref <- list()
for (position in mine){ #for snp position
  zero_one <- list()
  my_table <- all_lins2[all_lins2$Genome.position == position, ] #grab all rows with this position
  print(my_table)
  lineages_present <- my_table$Tuberculist.category
  lineages_present <- paste(lineages_present, sep="") #grab all the lineages where this snp is present (rename to fit with names)
  print(lineages_present)
  lineages_absent <- setdiff(x = colnames(all_lins2)[14:112], y = lineages_present) #grab all lineages where this snp is absent
  print(lineages_absent)
  zeros <- rep(0, length(lineages_absent)) #for every lineage where its absent add a zero
  ones <- rep(1, length(lineages_present)) #for every lineage where its present add a one
  zeros_ones <- c(zeros, ones)
  names(zeros_ones) <- c(lineages_absent, lineages_present)
  zeros_ones <- zeros_ones[colnames(all_lins2)[14:112]] #reorders the columns to the order of linesgaesin the original table
  empty_table_ref <- rbind(empty_table_ref, zeros_ones) #row = snp, add each row one at a time
}
rownames(empty_table_ref) <- mine

reference <- as.data.frame(empty_table_ref)

# II.C.5 80625 check

# cluster by synonymous and nonsynonymous SNPs
# add reference snps
snps_add <- match(rownames(reference), rownames(snps_nonsynonymous)) # find the snps in both dataframes
snps_add <- snps_add[!is.na(snps_add)]
snps_add <- snps_nonsynonymous[snps_add,] # subset snps by common positions

reference2 <- reference
reference2 <- reference2[rownames(reference2) %in% rownames(snps_add), ]
reference2 <- cbind(reference2, snps_add)

snps_uncommon <- subset(reference, !rownames(reference) %in% rownames(snps_nonsynonymous))
snps_uncommon[, 100:113] <- 0
colnames(snps_uncommon)[100:113] <- colnames(snps_nonsynonymous)
snps_uncommon_again <- subset(snps_nonsynonymous, !rownames(snps_nonsynonymous) %in% rownames(reference))
snps_uncommon_again[, 15:113] <- 0
colnames(snps_uncommon_again)[15:113] <- colnames(reference)

# combine all
reference3 <- rbind(reference2, snps_uncommon, snps_uncommon_again)

##### dendrogram of global dataset #####

snps_dend <- reference3
snps_dend <- as.data.frame(t(snps_dend))
d <- dist(snps_dend, method = "euclidean") # distance matrix
fit_hc <- hclust(d, method="ward")

plot(fit_hc, hang = -1, 
     main = "Cluster dendrogram", sub = NULL,
     xlab = "Mtb Strain", ylab = "Height", cex=0.85,
     horiz = TRUE) # display dendogram

fit_hc <- as.dendrogram(fit_hc)

colors = c("gray52", "#003366", "#CC5555")
# colors = c("gray47", "#5566AA", "#F17B4B")
fit_hc <- hclust(d, method="ward")
clus3 = cutree(fit_hc, 3)

fit_hc <- hclust(d, method="ward")

png("dendrogram_alllineages_thesis.png", height = 1500, width = 900)
plot(as.phylo(fit_hc), cex = 1.5, label.offset = 0.5, tip.color = colors[clus3])
dev.off()

png("all_lineages_thesis1.png", height = 3250, width = 3750, res = 250)
plot(as.phylo(fit_hc), cex = 1.25, label.offset = 4.5, tip.color = colors[clus3],
     type = "unrooted", no.margin = TRUE)
dev.off()

# read in functional categories file from mycobase
func_ann <- read.delim("Biocyc_functionalcat.txt", header=TRUE)
# get gene information
snp_dataset3 <- snp_dataset
snp_dataset3$INFO <- sapply(snp_dataset3$INFO, `[`, 4)
snp_dataset3$new <- func_ann$Pathways.of.gene[match(snp_dataset3$INFO, func_ann$Name)]
# check func_ann[func_ann$INFO == "Rv3922c", ]

##### edit functional categories #####

# clean up category annotation
# split by backspaces
categories <- func_ann
modsplit=function(x)strsplit(x, "//")
rowmodsplit <- function(x){
  lapply(as.character(x), modsplit)
}
categories$cat <- lapply(categories$Pathways.of.gene, rowmodsplit)

categories$cat <- sapply(categories$cat, `[[`, 1)

# get gene information
snp_dataset3 <- snp_dataset
snp_dataset3$INFO <- sapply(snp_dataset3$INFO, `[`, 4)
snp_dataset3$new <- categories$cat[match(snp_dataset3$INFO, categories$Name)]

allsnps_functions <- snp_dataset3
# allsnps_function <- allsnps_functions[allsnps_functions$POS %in% rownames(snps_all2), ]
allsnps_functions <- allsnps_functions[allsnps_functions$POS %in% rownames(snps_all2), ]
snps_nonsynonymous_fun <- allsnps_functions[allsnps_functions$POS %in% rownames(snps_nonsynonymous),]

# split by strain
snpsep <- split(snps_nonsynonymous_fun, snps_nonsynonymous_fun$Strain)

combine_cats <- function(x){
  df <- as.data.frame(snpsep[x])
  df <- unnest(df, cols = 10)
  df$snp_cat <- sapply(df[, 10],function(y) gsub(" ","", y)) # remove spaces
  df <- df[!duplicated(df[c(1,3,4,9,11)]),] # ensure no duplicates of category and snp
  df_counts <- data.frame(counts = summary(as.factor(df$snp_cat)))
  df_counts$pathways <- rownames(df_counts)
  return(df_counts)
}

snp119_counts <- combine_cats(1)
snp173_counts <- combine_cats(2)
snp212_counts <- combine_cats(3)
snp232_counts <- combine_cats(4)
snp281_counts <- combine_cats(5)
snp293_counts <- combine_cats(6)
snp318_counts <- combine_cats(7)
snp345_counts <- combine_cats(8)
snp346_counts <- combine_cats(9)
snp372_counts <- combine_cats(10)
snp374_counts <- combine_cats(11)
snp440_counts <- combine_cats(12)
snp639_counts <- combine_cats(13)
snp649_counts <- combine_cats(14)

# combine similar categories
get_total <- function(strain_name){
  if(sum(strain_name[grep("adenineandadenosinesalvage*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(strain_name[grep("adenineandadenosinesalvage*", rownames(strain_name)), ]$counts), pathways = "adenineandadenosinesalvage") # add counts for everything matching this string
    x <- strain_name[-grep("adenineandadenosinesalvage*", rownames(strain_name)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("cholesteroldegradation*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("cholesteroldegradation*", rownames(x)), ]$counts), pathways = "cholesteroldegradation") # add counts for everything matching this string
    x <- x[-grep("cholesteroldegradation*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("aminoimidazoleribonucleotidebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("aminoimidazoleribonucleotidebiosynthesis*", rownames(x)), ]$counts), pathways = "aminoimidazoleribonucleotidebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("aminoimidazoleribonucleotidebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("adenosinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("adenosinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]$counts), pathways = "adenosinedeoxyribonucleotides_denovo_biosynthesis") # add counts for everything matching this string
    x <- x[-grep("adenosinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("CDP-diacylglycerolbiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("CDP-diacylglycerolbiosynthesis*", rownames(x)), ]$counts), pathways = "CDP-diacylglycerolbiosynthesis") # add counts for everything matching this string
    x <- x[-grep("CDP-diacylglycerolbiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("serineandglycinebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("serineandglycinebiosynthesis*", rownames(x)), ]$counts), pathways = "serineandglycinebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("serineandglycinebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("L-argininebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("L-argininebiosynthesis*", rownames(x)), ]$counts), pathways = "L-argininebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("L-argininebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("isoleucinebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("isoleucinebiosynthesis*", rownames(x)), ]$counts), pathways = "isoleucinebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("isoleucinebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("L-prolinebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("L-prolinebiosynthesis*", rownames(x)), ]$counts), pathways = "L-prolinebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("L-prolinebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("NADHtocytochrome<i>bd</i>oxidaseelectrontransfer*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("NADHtocytochrome<i>bd</i>oxidaseelectrontransfer*", rownames(x)), ]$counts), pathways = "NADHtocytochrome_bd_oxidaseelectrontransfer") # add counts for everything matching this string
    x <- x[-grep("NADHtocytochrome<i>bd</i>oxidaseelectrontransfer*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("pyrimidinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("pyrimidinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]$counts), pathways = "pyrimidinedeoxyribonucleotides_denovo_biosynthesis") # add counts for everything matching this string
    x <- x[-grep("pyrimidinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("pyruvatefermentationtoacetate", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("pyruvatefermentationtoacetate*", rownames(x)), ]$counts), pathways = "pyruvatefermentationtoacetate") # add counts for everything matching this string
    x <- x[-grep("pyruvatefermentationtoacetate*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("trehalosebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("trehalosebiosynthesis*", rownames(x)), ]$counts), pathways = "trehalosebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("trehalosebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("mycolatebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("mycolatebiosynthesis*", rownames(x)), ]$counts), pathways = "mycolatebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("mycolatebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("superpathwayofguanosinenucleotides<i>denovo</i>biosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("superpathwayofguanosinenucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]$counts), pathways = "superpathwayofguanosinenucleotides_denovo_biosynthesis") # add counts for everything matching this string
    x <- x[-grep("superpathwayofguanosinenucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("guanosinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("guanosinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]$counts), pathways = "guanosinedeoxyribonucleotides_denovo_biosynthesis") # add counts for everything matching this string
    x <- x[-grep("guanosinedeoxyribonucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("superpathwayofadenosinenucleotides<i>denovo</i>biosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("superpathwayofadenosinenucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]$counts), pathways = "superpathwayofadenosinenucleotides_denovo_biosynthesis") # add counts for everything matching this string
    x <- x[-grep("superpathwayofadenosinenucleotides<i>denovo</i>biosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("phosphatidylglycerolbiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("phosphatidylglycerolbiosynthesis*", rownames(x)), ]$counts), pathways = "phosphatidylglycerolbiosynthesis") # add counts for everything matching this string
    x <- x[-grep("phosphatidylglycerolbiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("pentosephosphatepathway*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("pentosephosphatepathway*", rownames(x)), ]$counts), pathways = "pentosephosphatepathway") # add counts for everything matching this string
    x <- x[-grep("pentosephosphatepathway*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("L-tyrosinebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("L-tyrosinebiosynthesis*", rownames(x)), ]$counts), pathways = "L-tyrosinebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("L-tyrosinebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("guanineandguanosinesalvage*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("guanineandguanosinesalvage*", rownames(x)), ]$counts), pathways = "guanineandguanosinesalvage") # add counts for everything matching this string
    x <- x[-grep("guanineandguanosinesalvage*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("L-methioninebiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("L-methioninebiosynthesis*", rownames(x)), ]$counts), pathways = "L-methioninebiosynthesis") # add counts for everything matching this string
    x <- x[-grep("L-methioninebiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("L-tryptophanbiosynthesis*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("L-tryptophanbiosynthesis*", rownames(x)), ]$counts), pathways = "L-tryptophanbiosynthesis") # add counts for everything matching this string
    x <- x[-grep("L-tryptophanbiosynthesis*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  if(sum(strain_name[grep("trehalosedegradation*", rownames(strain_name)), ]$counts) > 0){
    extra_counts <- data.frame(counts = sum(x[grep("trehalosedegradation*", rownames(x)), ]$counts), pathways = "trehalosedegradation") # add counts for everything matching this string
    x <- x[-grep("trehalosedegradation*", rownames(x)), ]
    x[length(rownames(x))+1,] <- extra_counts
  }
  print(x)
}

snp119_counts <- get_total(snp119_counts)
snp173_counts <- get_total(snp173_counts)
snp212_counts <- get_total(snp212_counts)
snp232_counts <- get_total(snp232_counts)
snp281_counts <- get_total(snp281_counts)
snp293_counts <- get_total(snp293_counts)
snp318_counts <- get_total(snp318_counts)
snp345_counts <- get_total(snp345_counts)
snp346_counts <- get_total(snp346_counts)
snp372_counts <- get_total(snp372_counts)
snp374_counts <- get_total(snp374_counts)
snp440_counts <- get_total(snp440_counts)
snp639_counts <- get_total(snp639_counts)
snp649_counts <- get_total(snp649_counts)

# add strain and lineage information to each object
snp119_counts$strain <- "119"
snp173_counts$strain <- "173"
snp212_counts$strain <- "212"
snp232_counts$strain <- "232"
snp281_counts$strain <- "281"
snp293_counts$strain <- "293"
snp318_counts$strain <- "318"
snp345_counts$strain <- "345"
snp346_counts$strain <- "346"
snp372_counts$strain <- "372"
snp374_counts$strain <- "374"
snp440_counts$strain <- "440"
snp639_counts$strain <- "639"
snp649_counts$strain <- "649"

snp119_counts$order_strain <- 5
snp173_counts$order_strain <- 13
snp212_counts$order_strain <- 6
snp232_counts$order_strain <- 3
snp281_counts$order_strain <- 1
snp293_counts$order_strain <- 12
snp318_counts$order_strain <- 14
snp345_counts$order_strain <- 8
snp346_counts$order_strain <- 4
snp372_counts$order_strain <- 2
snp374_counts$order_strain <- 7
snp440_counts$order_strain <- 10
snp639_counts$order_strain <- 11
snp649_counts$order_strain <- 9

cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
  rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

snp_cats <- rbind(snp119_counts, snp173_counts, snp212_counts, snp232_counts, snp281_counts, snp293_counts, snp318_counts, snp345_counts, snp346_counts, snp372_counts, snp374_counts, snp440_counts, snp639_counts, snp649_counts)
snp_cats$lineage <- c(rep("lin2", length(snp119_counts$counts)), rep("lin4", length(snp173_counts$counts)), rep("lin2", length(snp212_counts$counts)), rep("lin1", length(snp232_counts$counts)), rep("lin1", length(snp281_counts$counts)), rep("lin4", length(snp293_counts$counts)), rep("lin4", length(snp318_counts$counts)), rep("lin2", length(snp345_counts$counts)), rep("lin1", length(snp346_counts$counts)), rep("lin1", length(snp372_counts$counts)), rep("lin2", length(snp374_counts$counts)), rep("lin4", length(snp440_counts$counts)), rep("lin4", length(snp639_counts$counts)), rep("lin2", length(snp649_counts$counts)))
# snp_cats2 <- cbind.fill(snp119_counts$counts, snp173_counts$counts, snp212_counts$counts, snp232_counts$counts, snp281_counts$counts, snp293_counts$counts, snp318_counts$counts, snp345_counts$counts, snp346_counts$counts, snp372_counts$counts, snp374_counts$counts, snp440_counts$counts, snp639_counts$counts, snp649_counts$counts)
# snp_cats$category <- rownames(snp_cats)
# # write.csv(snp_cats2, "snp_functional_biocyc_thesis.csv")

condition <- c("119", "173", "212", "232", "281", "239", "318", "345", "346", "372", "374", "440", "639", "649")
order_lin <- c("346", "232", "372", "281", "119", "212", "345", "374", "649", "173", "239", "440", "318", "639")

snp_cats <- snp_cats[!snp_cats$counts == "NA's", ]
snp_cats <- snp_cats[!snp_cats$counts == "NA2", ]
sum(snp_cats$counts)
snp_others <- snp_cats[snp_cats$pathways == "(Other)", ]
sum(snp_others$counts)
snp_cats <- snp_cats[!snp_cats$pathways == "(Other)", ]
sum(snp_cats$counts)
snp_cats <- snp_cats[snp_cats$counts > 5, ]
sum(snp_cats$counts)

# write.csv(summary(factor(snp_cats$pathways)), "functional_edit_thesis.csv")

# get edited names for the function categories
snp_cats_edited <- read.csv("functional_edited_aminoacids.csv")

idx <- match(snp_cats$pathways, snp_cats_edited$X)

snp_cats$pathways_superfamily <- snp_cats_edited$Category[idx]

colourCount = length(unique(snp_cats$pathways_superfamily))

##### Load metabolite files #####
##### positive ion mode
metab_pos <- read.csv("metabolites_pos_R.csv", header=TRUE)

rownames(metab_pos) <- make.names(metab_pos$Mass, unique=TRUE)
colnames(metab_pos) <- gsub(".raw.", " ", colnames(metab_pos))
metab_pos <- metab_pos[, c(1:84, 86:92)] # remove 649_6
metab_rt_pos <- metab_pos
metab_pos <- metab_pos[, 2:91]
metab_pos2 <- log(metab_pos, 2)

##### negative ion mode
metab_neg <- read.csv("metabolites_neg_R.csv")
rownames(metab_neg) <- make.names(metab_neg$Mass, unique=TRUE)
colnames(metab_neg) <- gsub(".raw.", " ", colnames(metab_neg))
metab_rt_neg <- metab_neg
metab_neg <- metab_neg[, 2:90]
metab_neg2 <- log(metab_neg, 2)

# average metabolite levels for each strain
# avg <- rowMeans(metab_pos[grep("X119_", colnames(metab_pos))]) # test
pos_avg <- sapply(c("X119", "X173", "X212", "X232", "X281", "X293", "X318", "X345", "X346", "X372", "X374", "X440", "X639", "X649", "H37Rv"), function(x) rowMeans(metab_pos[grep(x, colnames(metab_pos))]))
neg_avg <- sapply(c("X119", "X173", "X212", "X232", "X281", "X293", "X318", "X345", "X346", "X372", "X374", "X440", "X639", "X649", "H37Rv"), function(x) rowMeans(metab_neg[grep(x, colnames(metab_neg))]))
plog_avg <- sapply(c("X119", "X173", "X212", "X232", "X281", "X293", "X318", "X345", "X346", "X372", "X374", "X440", "X639", "X649", "H37Rv"), function(x) rowMeans(metab_pos2[grep(x, colnames(metab_pos2))]))
nlog_avg <- sapply(c("X119", "X173", "X212", "X232", "X281", "X293", "X318", "X345", "X346", "X372", "X374", "X440", "X639", "X649", "H37Rv"), function(x) rowMeans(metab_neg2[grep(x, colnames(metab_neg2))]))

##### plot all targeted metabolites IGNORE #######

setwd("D:/Lineage/Data_run/results/")

metab_plot <- function(x, y){
  ggplot(targeted_avg_bar, aes(x = rownames(targeted_avg_bar), y = targeted_avg_bar[,x], fill = rownames(targeted_avg_bar))) +
    geom_errorbar(aes(ymin=targeted_avg_bar[,x]-targeted_avg_sd[,x], ymax=targeted_avg_bar[,x]+targeted_avg_sd[,x]), width=.2, size = 0.75, color = c("#CC6666", "gray67", "#003366")) +
    geom_bar(stat = "identity", width = 0.65) +
    xlab("") +
    ylab(y) +
    # geom_text(aes(label=non_syn), vjust=1.6, color="white", size=3.5) +
    scale_fill_manual("Lineages", values=c("#CC6666", "gray67", "#003366"), aesthetics = c("colour", "fill")) +
    theme_minimal() +
    theme(axis.title.x = element_text(size=30, margin=margin(7, 1, 1, 1)),
          axis.title.y = element_text(size=30, margin=margin(1, 10, 1, 1)),
          axis.text.x = element_text(size=28, hjust=0.5, vjust=0.5, colour = "gray25"),
          axis.text.y = element_text(size=28, colour = "gray25"),
          legend.title = element_text(size=30),
          legend.text = element_text(size=29),
          axis.line = element_line(colour="black"))
  # scale_y_continuous(expand = c(0.005, 10), limits = c(0, 29))
}

{
png("bar_transaconitate.png", width = 600, height = 575)
metab_plot(x = 1, y = "Abundance of trans aconitate")
dev.off()
png("bar_itaconicacid.png", width = 600, height = 575)
metab_plot(x = 2, y = "Abundance of itaconic acid")
dev.off()
png("bar_oxoglutaricacid.png", width = 600, height = 575)
metab_plot(x = 3, y = "Abundance of oxoglutaric acid")
dev.off()
png("bar_succinicacid.png", width = 600, height = 575)
metab_plot(x = 4, y = "Abundance of succinic acid")
dev.off()
png("bar_fumaricacid.png", width = 600, height = 575)
metab_plot(x = 5, y = "Abundance of fumaric acid")
dev.off()
png("bar_pyruvicacid.png", width = 600, height = 575)
metab_plot(x = 6, y = "Abundance of pyruvic acid")
dev.off()
png("bar_lactate.png", width = 600, height = 575)
metab_plot(x = 7, y = "Abundance of lactate")
dev.off()
png("bar_malicacid.png", width = 600, height = 575)
metab_plot(x = 8, y = "Abundance of malic acid")
dev.off()
png("bar_citricacid.png", width = 600, height = 575)
metab_plot(x = 9, y = "Abundance of citric acid")
dev.off()
png("bar_glycerylphosphotylethanolamine.png", width = 600, height = 575)
metab_plot(x = 10, y = "Abundance of glycerylphosphotylethanolamine")
dev.off()
png("bar_lasparticacid.png", width = 600, height = 575)
metab_plot(x =11, y = "Abundance of L-aspartic acid")
dev.off()
png("bar_lglutamate.png", width = 600, height = 575)
metab_plot(x = 12, y = "Abundance of L-glutmate")
dev.off()
png("bar_lproline.png", width = 600, height = 575)
metab_plot(x = 13, y = "Abundance of L-proline")
dev.off()
png("bar_lphenylalanine.png", width = 600, height = 575)
metab_plot(x = 14, y = "Abundance of L-phenylalanine")
dev.off()
png("bar_lserine.png", width = 600, height = 575)
metab_plot(x = 15, y = "Abundance of L-serine")
dev.off()
png("bar_ltyrosine.png", width = 600, height = 575)
metab_plot(x = 16, y = "Abundance of L-tyrosine")
dev.off()
png("bar_lthreonine.png", width = 600, height = 575)
metab_plot(x = 17, y = "Abundance of L-threonine")
dev.off()
png("bar_ltryptophan.png", width = 600, height = 575)
metab_plot(x = 18, y = "Abundance of L-tryptophan")
dev.off()
png("bar_glutamine.png", width = 600, height = 575)
metab_plot(x = 19, y = "Abundance of glutamine")
dev.off()
png("bar_glycine.png", width = 600, height = 575)
metab_plot(x = 20, y = "Abundance of glycine")
dev.off()
png("bar_lisoleucine.png", width = 600, height = 575)
metab_plot(x = 21, y = "Abundance of L-isoleucine")
dev.off()
png("bar_lvaline.png", width = 600, height = 575)
metab_plot(x = 22, y = "Abundance of L-valine")
dev.off()
png("bar_lleucine.png", width = 600, height = 575)
metab_plot(x = 23, y = "Abundance of L-leucine")
dev.off()
png("bar_lalanine.png", width = 600, height = 575)
metab_plot(x = 24, y = "Abundance of L-alanine")
dev.off()
png("bar_citrulline.png", width = 600, height = 575)
metab_plot(x = 25, y = "Abundance of citrulline")
dev.off()
png("bar_nacetyllornithine.png", width = 600, height = 575)
metab_plot(x = 26, y = "Abundance of N-acetyl L-ornithine")
dev.off()
png("bar_NAD.png", width = 600, height = 575)
metab_plot(x = 27, y = "Abundance of NAD")
dev.off()
png("bar_adenine.png", width = 600, height = 575)
metab_plot(x = 28, y = "Abundance of adenine")
dev.off()
png("bar_aminobutryicacid.png", width = 600, height = 575)
metab_plot(x = 29, y = "Abundance of aminobutryic acid")
dev.off()
png("bar_lhistidine.png", width = 600, height = 575)
metab_plot(x = 30, y = "Abundance of L-histidine")
dev.off()
png("bar_larginine.png", width = 600, height = 575)
metab_plot(x = 31, y = "Abundance of L-arginine")
dev.off()
png("bar_llysine.png", width = 600, height = 575)
metab_plot(x = 32, y = "Abundance of L-lysine")
dev.off()
png("bar_inosine.png", width = 600, height = 575)
metab_plot(x = 33, y = "Abundance of inosine")
dev.off()
png("bar_AMP.png", width = 600, height = 575)
metab_plot(x = 34, y = "Abundance of adenosine monophosphate")
dev.off()
png("bar_ergothionine.png", width = 600, height = 575)
metab_plot(x = 35, y = "Abundance of ergothionine")
dev.off()
png("bar_laminoadipicacid.png", width = 600, height = 575)
metab_plot(x = 36, y = "Abundance of L-aminoadipic acid")
dev.off()
png("bar_diaminopimelicacid.png", width = 600, height = 575)
metab_plot(x = 37, y = "Abundance of diaminopimelic acid")
dev.off()
png("bar_lalinine1.png", width = 600, height = 575)
metab_plot(x = 38, y = "Abundance of L-alinine 1")
dev.off()
png("bar_dalanyldalanine.png", width = 600, height = 575)
metab_plot(x = 39, y = "Abundance of D-alanyl D-alanine")
dev.off()
png("bar_llysine1.png", width = 600, height = 575)
metab_plot(x = 40, y = "Abundance of L-lysine 1")
dev.off()
}

png("metabolitelevels_diaminopimelicacid_avg_paper.png", width = 600, height = 575)

ggplot(targeted_avg_bar, aes(x = rownames(targeted_avg_bar), y = targeted_avg_bar[,37], fill = rownames(targeted_avg_bar))) +
  geom_bar(stat = "identity", width = 0.65) +
  geom_errorbar(aes(ymin=targeted_avg_bar[,37]-targeted_avg_sd[,37], ymax=targeted_avg_bar[,37]+targeted_avg_sd[,37]), width=.2, size = 0.75, color = c("#CC6666", "gray67", "#003366")) +
  xlab("") +
  ylab("Abundance of diaminopimelic acid") +
  # geom_text(aes(label=non_syn), vjust=1.6, color="white", size=3.5) +
  scale_fill_manual("Lineages", values=c("#CC6666", "gray67", "#003366"), aesthetics = c("colour", "fill")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size=30, margin=margin(7, 1, 1, 1)),
        axis.title.y = element_text(size=30, margin=margin(1, 10, 1, 1)),
        axis.text.x = element_text(size=28, hjust=0.5, vjust=0.5, colour = "gray25"),
        axis.text.y = element_text(size=28, colour = "gray25"),
        legend.title = element_text(size=30),
        legend.text = element_text(size=29),
        axis.line = element_line(colour="black"))
# scale_y_continuous(expand = c(0.005, 10), limits = c(0, 29))

dev.off()

targeted_avg_bar <- as.data.frame(t(targeted_avg_bar))
targeted_avg_bar$lin <- c("Lin1", "Lin2", "Lin4")
bar <- melt(targeted_avg_bar)

png("metabolitelevels_all_paper.png", width = 700, height = 575)

ggplot(bar, aes(x = lin, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  xlab("Strains") +
  ylab("Metabolite abundance (% of total)") +
  # geom_text(aes(label=non_syn), vjust=1.6, color="white", size=3.5) +
  # scale_fill_manual("Lineages", values=c("#CC6666", "gray67", "#003366")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size=20, margin=margin(1, 1, 1, 1)),
        axis.title.y = element_text(size=20, margin=margin(1, 7, 1, 1)),
        axis.text.x = element_text(size=20, hjust=0.1, vjust=0.5),
        axis.text.y = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=10),
        axis.line = element_line(colour="black"))
# scale_y_continuous(expand = c(0.005, 10), limits = c(0, 100))

dev.off()

##### Targeted metabolite DAA #####

###-- read the data #####

setwd("D:/Lineage/Data_run")
targeted <- read.csv("Targeted_all_strains_exp3_July20.csv")
targeted <- targeted[!targeted$X == 333, ]
targeted <- targeted[!targeted$X == 367, ]
targeted[sapply(targeted, is.factor)] <- lapply(targeted[sapply(targeted, is.factor)], function(x) as.numeric(as.character(x)))
rownames(targeted) <- make.names(targeted$X, unique = TRUE)
metab <- targeted
metab <- metab[!rownames(metab) == "X649.5", ]
metab[, 2:42] <- sapply(metab[, 2:42], as.numeric)

metab <- metab[!rownames(metab) == "X281.4", ]
metab <- metab[!rownames(metab) == "X372.4", ]
metab <- metab[!rownames(metab) == "X639.3", ]
metab <- metab[!rownames(metab) == "H37Rv.5", ]
metab <- metab[!rownames(metab) == "X119.5", ]

###-- log and normalise #####

# log2 transform
metab <- log(metab[, 2:42], 2)
metab[metab == "-Inf"] <- 0

# subset based on lineage number
x119 <- metab[tidyselect::vars_select(rownames(metab), starts_with('X119')), ]
x212 <- metab[tidyselect::vars_select(rownames(metab), starts_with('X212')), ]
x173 <- metab[tidyselect::vars_select(rownames(metab), starts_with('X173')), ]
x232 <- metab[tidyselect::vars_select(rownames(metab), starts_with('X232')), ]
x281 <- metab[tidyselect::vars_select(rownames(metab), starts_with('X281')), ]
x293 <- metab[tidyselect::vars_select(rownames(metab), starts_with('X293')), ]
x318 <- metab[tidyselect::vars_select(rownames(metab), starts_with('X318')), ]
x345 <- metab[tidyselect::vars_select(rownames(metab), starts_with('X345')), ]
x346 <- metab[tidyselect::vars_select(rownames(metab), starts_with('X346')), ]
x372 <- metab[tidyselect::vars_select(rownames(metab), starts_with('X372')), ]
x374 <- metab[tidyselect::vars_select(rownames(metab), starts_with('X374')), ]
x440 <- metab[tidyselect::vars_select(rownames(metab), starts_with('X440')), ]
x639 <- metab[tidyselect::vars_select(rownames(metab), starts_with('X639')), ]
x649 <- metab[tidyselect::vars_select(rownames(metab), starts_with('X649')), ]
h37rv <- metab[tidyselect::vars_select(rownames(metab), starts_with('H37Rv')), ]

# differential expression of metabolite levels
names <- c("x119", "x119", "x119", "x119", "x119", "x119",
           "x173", "x173", "x173", "x173", "x173", "x173",
           "x212", "x212", "x212", "x212", "x212", "x212",
           "x232", "x232", "x232", "x232", "x232", "x232",
           "x281", "x281", "x281", "x281", "x281", "x281",
           "x293", "x293", "x293", "x293", "x293", "x293",
           "x318", "x318", "x318", "x318", "x318", "x318",
           "x345", "x345", "x345", "x345", "x345", "x345",
           "x346", "x346", "x346", "x346", "x346", "x346",
           "x372", "x372", "x372", "x372", "x372", "x372",
           "x374", "x374", "x374", "x374", "x374", "x374",
           "x440", "x440", "x440", "x440", "x440", "x440",
           "x639", "x639", "x639", "x639", "x639", "x639",
           "x649", "x649", "x649", "x649", "x649",
           "h37rv", "h37rv", "h37rv", "h37rv", "h37rv", "h37rv")

names <- c("x119", "x119", "x119", "x119", "x119",
           "x173", "x173", "x173", "x173", "x173", "x173",
           "x212", "x212", "x212", "x212", "x212", "x212",
           "x232", "x232", "x232", "x232", "x232", "x232",
           "x281", "x281", "x281", "x281", "x281",
           "x293", "x293", "x293", "x293", "x293", "x293",
           "x318", "x318", "x318", "x318", "x318", "x318",
           "x345", "x345", "x345", "x345", "x345", "x345",
           "x346", "x346", "x346", "x346", "x346", "x346",
           "x372", "x372", "x372", "x372", "x372",
           "x374", "x374", "x374", "x374", "x374", "x374",
           "x440", "x440", "x440", "x440", "x440", "x440",
           "x639", "x639", "x639", "x639", "x639",
           "x649", "x649", "x649", "x649", "x649",
           "h37rv", "h37rv", "h37rv", "h37rv", "h37rv")

###-- comparison against h37rv heatmap for paper #####

# only look at comparing against h37rv
# targeted ion mode
design <- model.matrix(~0 + factor(names))
colnames(design) <- c("h37rv", "x119", "x173", "x212", "x232", "x281", "x293", "x318", "x345", "x346", "x372", "x374", "x440", "x639", "x649")
fit_targeted <- lmFit(t(metab), design) #fit a linear model
contrast.matrix <- makeContrasts(x119-h37rv, x173-h37rv, x212-h37rv, x232-h37rv, x281-h37rv, x293-h37rv, x318-h37rv, x345-h37rv, x346-h37rv, x372-h37rv, x374-h37rv, x440-h37rv, x639-h37rv, x649-h37rv,
                                 levels=design)
fit.cb.targeted <- eBayes(contrasts.fit(fit_targeted, contrast.matrix))
p_values <- fit.cb.targeted$p.value
res_targeted <- topTable(fit.cb.targeted, sort="none", n=Inf) # gets all results

# select differentially expressed genes for each comparison individually and subset original data frame by these
significant_filter <- function(fit_data, x, dge_res){
  y <- topTable(fit_data, coef = x, sort="none", n=Inf) # gets significant x119 - h37rv results
  y <- y[y$adj.P.Val < 0.05, ]
  y[, 7] <- "119"
  y1 <- topTable(fit_data, coef = x+1, sort="none", n=Inf)
  y1 <- y1[y1$adj.P.Val < 0.05, ]
  y1[, 7] <- "173"
  y2 <- topTable(fit_data, coef = x+2, sort="none", n=Inf)
  y2 <- y2[y2$adj.P.Val < 0.05, ]
  y2[, 7] <- "212"
  y3 <- topTable(fit_data, coef = x+3, sort="none", n=Inf)
  y3 <- y3[y3$adj.P.Val < 0.05, ]
  y3[, 7] <- "232"
  y4 <- topTable(fit_data, coef = x+4, sort="none", n=Inf)
  y4 <- y4[y4$adj.P.Val < 0.05, ]
  y4[, 7] <- "281"
  y5 <- topTable(fit_data, coef = x+5, sort="none", n=Inf)
  y5 <- y5[y5$adj.P.Val < 0.05, ]
  y5[, 7] <- "293"
  y6 <- topTable(fit_data, coef = x+6, sort="none", n=Inf)
  y6 <- y6[y6$adj.P.Val < 0.05, ]
  y6[, 7] <- "318"
  y7 <- topTable(fit_data, coef = x+7, sort="none", n=Inf)
  y7 <- y7[y7$adj.P.Val < 0.05, ]
  y7[, 7] <- "345"
  y8 <- topTable(fit_data, coef = x+8, sort="none", n=Inf)
  y8 <- y8[y8$adj.P.Val < 0.05, ]
  y8[, 7] <- "346"
  y9 <- topTable(fit_data, coef = x+9, sort="none", n=Inf)
  y9 <- y9[y9$adj.P.Val < 0.05, ]
  y9[, 7] <- "372"
  y10 <- topTable(fit_data, coef = x+10, sort="none", n=Inf)
  y10 <- y10[y10$adj.P.Val < 0.05, ]
  y10[, 7] <- "374"
  y11 <- topTable(fit_data, coef = x+11, sort="none", n=Inf)
  y11 <- y11[y11$adj.P.Val < 0.05, ]
  y11[, 7] <- "440"
  y12 <- topTable(fit_data, coef = x+12, sort="none", n=Inf)
  y12 <- y12[y12$adj.P.Val < 0.05, ]
  y12[, 7] <- "639"
  y13 <- topTable(fit_data, coef = x+13, sort="none", n=Inf)
  y13 <- y13[y13$adj.P.Val < 0.05, ]
  y13[, 7] <- "649"
  y_final <- list(y, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13) # groups everything
  fil1 <- rownames(y_final[[1]]) # get the metabolite number which is significant
  fil2 <- rownames(y_final[[2]])
  fil3 <- rownames(y_final[[3]])
  fil4 <- rownames(y_final[[4]])
  fil5 <- rownames(y_final[[5]])
  fil6 <- rownames(y_final[[6]])
  fil7 <- rownames(y_final[[7]])
  fil8 <- rownames(y_final[[8]])
  fil9 <- rownames(y_final[[9]])
  fil10 <- rownames(y_final[[10]])
  fil11 <- rownames(y_final[[11]])
  fil12 <- rownames(y_final[[12]])
  fil13 <- rownames(y_final[[13]])
  fil14 <- rownames(y_final[[14]])
  filter_by <- c(fil1, fil2, fil3, fil4, fil5, fil6, fil7, fil8, fil9, fil10, fil11, fil12, fil13, fil14) # groups metabolite numbers
  print(length(filter_by))
  filter_by <- unique(filter_by)
  print(length(filter_by))
  return(dge_res[rownames(dge_res) %in% filter_by, ])
}

# get the significant metabolites for the individual strains
significant_filter_metabolites <- function(fit_data, x, dge_res){
  y <- topTable(fit_data, coef = x, sort="none", n=Inf) # gets significant x119 - h37rv results
  y <- y[y$adj.P.Val < 0.05, ]
  y[, 7] <- "119"
  y1 <- topTable(fit_data, coef = x+1, sort="none", n=Inf)
  y1 <- y1[y1$adj.P.Val < 0.05, ]
  y1[, 7] <- "173"
  y2 <- topTable(fit_data, coef = x+2, sort="none", n=Inf)
  y2 <- y2[y2$adj.P.Val < 0.05, ]
  y2[, 7] <- "212"
  y3 <- topTable(fit_data, coef = x+3, sort="none", n=Inf)
  y3 <- y3[y3$adj.P.Val < 0.05, ]
  y3[, 7] <- "232"
  y4 <- topTable(fit_data, coef = x+4, sort="none", n=Inf)
  y4 <- y4[y4$adj.P.Val < 0.05, ]
  y4[, 7] <- "281"
  y5 <- topTable(fit_data, coef = x+5, sort="none", n=Inf)
  y5 <- y5[y5$adj.P.Val < 0.05, ]
  y5[, 7] <- "293"
  y6 <- topTable(fit_data, coef = x+6, sort="none", n=Inf)
  y6 <- y6[y6$adj.P.Val < 0.05, ]
  y6[, 7] <- "318"
  y7 <- topTable(fit_data, coef = x+7, sort="none", n=Inf)
  y7 <- y7[y7$adj.P.Val < 0.05, ]
  y7[, 7] <- "345"
  y8 <- topTable(fit_data, coef = x+8, sort="none", n=Inf)
  y8 <- y8[y8$adj.P.Val < 0.05, ]
  y8[, 7] <- "346"
  y9 <- topTable(fit_data, coef = x+9, sort="none", n=Inf)
  y9 <- y9[y9$adj.P.Val < 0.05, ]
  y9[, 7] <- "372"
  y10 <- topTable(fit_data, coef = x+10, sort="none", n=Inf)
  y10 <- y10[y10$adj.P.Val < 0.05, ]
  y10[, 7] <- "374"
  y11 <- topTable(fit_data, coef = x+11, sort="none", n=Inf)
  y11 <- y11[y11$adj.P.Val < 0.05, ]
  y11[, 7] <- "440"
  y12 <- topTable(fit_data, coef = x+12, sort="none", n=Inf)
  y12 <- y12[y12$adj.P.Val < 0.05, ]
  y12[, 7] <- "639"
  y13 <- topTable(fit_data, coef = x+13, sort="none", n=Inf)
  y13 <- y13[y13$adj.P.Val < 0.05, ]
  y13[, 7] <- "649"
  y[, 8] <- rownames(y) # get the metabolite number which is significant
  y1[, 8] <- rownames(y1)
  y2[, 8] <- rownames(y2)
  y3[, 8] <- rownames(y3)
  y4[, 8] <- rownames(y4)
  y5[, 8] <- rownames(y5)
  y6[, 8] <- rownames(y6)
  y7[, 8] <- rownames(y7)
  y8[, 8] <- rownames(y8)
  y9[, 8] <- rownames(y9)
  y10[, 8] <- rownames(y10)
  y11[, 8] <- rownames(y11)
  y12[, 8] <- rownames(y12)
  y13[, 8] <- rownames(y13)
  y_final <- rbind(y3, y4, y8, y9, y, y2, y10, y7, y13, y1, y5, y6, y11, y12) # groups everything
  print(y_final)
}

sig_res_targeted <- significant_filter(fit.cb.targeted, 1, res_targeted)
sig_res_metabs <- significant_filter_metabolites(fit.cb.targeted, 1, res_targeted)

index <- which(c(sig_res_metabs$logFC > 0) == TRUE)
sig_res_metabs2 <- sig_res_metabs
sig_res_metabs2$logFC[index] <- "POS"
index <- which(c(sig_res_metabs$logFC < 0) == TRUE)
sig_res_metabs2$logFC[index] <- "NEG"
colnames(sig_res_metabs2) <- c(colnames(sig_res_metabs2)[1:6], "strain", "metabolite")

sig_res_metabs3 <- pivot_wider(sig_res_metabs2[, c(1, 7:8)], names_from = strain, values_from = logFC)

write.csv(sig_res_metabs3, "sig_res_metabs3_thesis.csv")

ggplot(res_targeted) +
  geom_point(aes(x=x119...h37rv, y=-log10(adj.P.Val))) +
  #geom_text(aes(x = lin1...lin2, y = -log10(adj.P.Val), label = ifelse(res$sig2 == T, rownames(res),"")), size=2.5) +
  ggtitle("lin2vslin4") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_hline(yintercept=2, linetype="dashed", color="black", size=1) +
  geom_vline(xintercept=2, linetype="dashed", color="black", size=1) +
  geom_vline(xintercept= -2, linetype="dashed", color="black", size=1) +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

write.csv(sig_res_targeted, "sig_res_targeted_thesis.csv")

# plot as heatmap
mat <- sig_res_targeted[!rownames(sig_res_targeted) == "GABA", ]
colnames(mat) <- gsub("x", "s", colnames(mat))
colnames(mat) <- gsub("[..].*", "", colnames(mat))
mat <- mat[, 1:14]
rownames(mat)[rownames(mat) == "N2.Acetyl.L.ornithine"] <- "Acetyl.ornithine"
rownames(mat)[rownames(mat) == "Nicotinamide.adenine.dinucleotide..NAD."] <- "Nicotinamide.adeenine.dinucleotide"
rownames(mat)[rownames(mat) == "L.2.Aminoadipic.acid"] <- "Aminoadipic.acid"
rownames(mat)[rownames(mat) == "D.Alanyl.D.alanine"] <- "Alanyl.alanine"


colours_lin <- c("#999999", "#003366", "#999999", "#CC6666", "#CC6666", "#003366", "#003366", "#999999", "#CC6666", "#CC6666", "#999999", "#003366", "#003366", "#999999", "#003366")

colours2 <- as.data.frame(colours_lin[-15])
lineages_colours <- list(colours3=colours_lin[-15])

png("targeted_metabolites_DAA2.png", width=2250, height=2000, res = 300)

# pheatmap(mat)
# h_anno <- rowAnnotation(df=colours2, col = list(colours3=colours_lin[-15][[1]]))
h_anno <- HeatmapAnnotation(df=colours2, col = lineages_colours)
# Heatmap(as.matrix(mat), cluster_rows = FALSE, top_annotation = h_anno)
Heatmap(as.matrix(mat), top_annotation = h_anno)

# lgd = Legend(title = "Lineages", at = c("Lin1", "Lin2", "Lin4"), legend_gp = gpar(fill = 1:3))

dev.off()

write.csv(mat, "heatmap_data_2023.csv")

# ergothioneine

# plot as heatmap
mat <- sig_res_targeted[rownames(sig_res_targeted) == "Ergothioneine", ]
colnames(mat) <- gsub("x", "s", colnames(mat))
colnames(mat) <- gsub("[..].*", "", colnames(mat))
mat <- mat[, 1:14]

colours_lin <- c("#999999", "#003366", "#999999", "#CC6666", "#CC6666", "#003366", "#003366", "#999999", "#CC6666", "#CC6666", "#999999", "#003366", "#003366", "#999999", "#003366")

colours2 <- as.data.frame(colours_lin[-15])
lineages_colours <- list(colours3=colours_lin[-15])

png("targeted_ergothioneine_heatmap.png", width=1000, height=225, res = 150)

# pheatmap(mat)
# h_anno <- rowAnnotation(df=colours2, col = list(colours3=colours_lin[-15][[1]]))
h_anno <- HeatmapAnnotation(df=colours2, col = lineages_colours)
# Heatmap(as.matrix(mat), cluster_rows = FALSE, top_annotation = h_anno)
Heatmap(as.matrix(mat), top_annotation = h_anno)

# lgd = Legend(title = "Lineages", at = c("Lin1", "Lin2", "Lin4"), legend_gp = gpar(fill = 1:3))

dev.off()

###-- pca and log_avg #####

# pca plot
names_targeted <- rownames(sig_res_targeted) # get names of sig metabolites

log_avg <- sapply(c("X119", "X173", "X212", "X232", "X281", "X293", "X318", "X345", "X346", "X372", "X374", "X440", "X639", "X649", "H37Rv"), function(x) colMeans(metab[grep(x, rownames(metab)), ]))
log_avg <- as.data.frame(log_avg)
log_avg2 <- log_avg[, -15]
sub <- log_avg2[! rownames(log_avg2) == "γ.Aminobutryic.acid",] # remove outlier metabolite
sub <- sub[! rownames(sub) == "GABA",] # remove outlier metabolite

# log_avg_ordered <- log_avg2[c(1:9, 14:36, 38:39, 10:13, 37), ]
log_avg_ordered <- log_avg2
colnames(log_avg_ordered) <- gsub("X", "s", colnames(log_avg_ordered))

sub <- as.data.frame(sub[rownames(sub) %in% names_targeted, ]) # subset original dataframe logfc values by significant metabolites
data <- sub
data <- t(data)
pca <- prcomp(data)

###-- correlation figures #####

# correlation plots

# get the metabolites that are significant
data.df2 <- log_avg_ordered[rownames(log_avg_ordered) %in% rownames(sig_res_targeted), ]

data.df2 <- rbind(snp_blocks, data.df2)
data.df2 <- data.df2[-4, ] # remove 4013 as strange values in this

results_cor <- rcorr(as.matrix(t(data.df2)))

# filter for p-value below 0.05

# results_cor_p <- as.data.frame(results_cor$P[116:135, 1:115]) # create data frame with metabolites as rows and blocks as columns
is_numeric <- is.na(as.numeric(rownames(results_cor$P))) # select all values and convert to T/F responses
is_trues <- which(is_numeric == TRUE) # select all numeric values
is_falses <- which(is_numeric == FALSE) # select all character values
results_cor_p <- as.data.frame(results_cor$P[(is_falses[length(is_falses)]+1):is_trues[length(is_trues)], 1:is_falses[length(is_falses)]]) # create data frame with metabolites as rows and blocks as columns
# results_cor_p <- as.data.frame(results_cor$P[116:151, 1:115]) # create data frame with metabolites as rows and blocks as columns
results_cor_p$metabolites <- rownames(results_cor_p) # add metabolites
results_cor_p <- as.data.frame(results_cor_p)
results_p <- melt(results_cor_p, "metabolites") # create list for filtering based on metabolites
length(unique(results_p$variable)) # number of snps
results_p <- results_p[results_p$value < 0.05, ]
length(unique(results_p$variable)) # number of snps with significant correlation
pos_cor_blocks <- unique(results_p$variable)

data.df2 <- log_avg_ordered[rownames(log_avg_ordered) %in% rownames(sig_res_targeted), ]

data.df2 <- rbind(snp_blocks[rownames(snp_blocks) %in% pos_cor_blocks, ], data.df2) # [, -15] to remove h37rv as not in snp data

results_cor2 <- t(data.df2)
results_circoplot <- rcorr(as.matrix(results_cor2))

# get correlation values
cor_df <- results_circoplot$r
cor_df <- as.data.frame(cor_df)

# get p values
cor_df_p <- results_circoplot$P
cor_df_p <- as.data.frame(cor_df_p)

filter_correlation <- function(x){
  indexing <- which(x < 0.5 & x > -0.5)
  y <- as.data.frame(x[indexing])
  y$blocks <- rownames(x[indexing])
  print(y)
}

lapply(cor_df, filter_correlation)

write.csv(results_circoplot$r, "results_circoplot_targeted_r_thesis.csv")
write.csv(results_circoplot$P, "results_circoplot_targeted_p_thesis.csv")

# find the SNP block and correlated metabolite(s)

results_cor_p2 <- as.data.frame(results_circoplot$P)
results_cor_r2 <- as.data.frame(results_circoplot$r)
results_cor_p2$`340132`
results_cor_r2$`340132`
snpblock <- results_cor_r2

# select for only high correlation values

# results_cor_sig_r <- as.data.frame(results_circoplot$r[102:137,1:101])
# results_cor_sig_r <- as.data.frame(results_circoplot$r[92:111, 1:91])
is_numeric <- is.na(as.numeric(rownames(results_circoplot$r))) # select all values and convert to T/F responses
is_trues <- which(is_numeric == TRUE) # select all numeric values i.e. metabolites
is_falses <- which(is_numeric == FALSE) # select all character values i.e. nsSNP blocks
results_cor_sig_r <- as.data.frame(results_circoplot$r[(is_falses[length(is_falses)]+1):is_trues[length(is_trues)], 1:is_falses[length(is_falses)]]) # create data frame with metabolites as rows and blocks as columns
# results_cor_sig_r <- as.data.frame(results_circoplot$r[103:138, 1:102])
results_cor_sig_r$metabolites <- rownames(results_cor_sig_r)
results_cor_sig_r <- as.data.frame(results_cor_sig_r)
results_r <- melt(results_cor_sig_r, "metabolites")
length(unique(results_r$variable)) # number of snps
results_r_pos50 <- results_r[results_r$value >= 0.75, ]
results_r_neg50 <- results_r[results_r$value <= -0.75, ]
results_r_50 <- rbind(results_r_pos50, results_r_neg50)
length(unique(results_r_50$variable)) # number of snps with significant correlation
pos_cor_sig_blocks <- unique(results_r_50$variable)

data.df2 <- log_avg_ordered[rownames(log_avg_ordered) %in% rownames(sig_res_targeted), ]

data.df2 <- rbind(snp_blocks[rownames(snp_blocks) %in% pos_cor_sig_blocks, ], data.df2)
# data.df2 <- data.df2[-4, ] # remove 4013 as all the same value

blocks_filtered <- t(data.df2)
blocks_filtered <- rcorr(as.matrix(blocks_filtered))
jpeg("correlation_plot_targeted_pvalue_cor_filtered_thesis.jpeg", height = 750, width = 2000)
corrplot(blocks_filtered$r[41:76,1:40], p.mat = blocks_filtered$P[41:76,1:40], sig.level = 0.05, insig = "blank", outline = TRUE, tl.col = "black", addgrid.col = NA, col = brewer.pal(n = 8, name = "RdBu"),
         mar = c(0.5, 1, 0.5, 1),
         tl.cex = 1.5,
         cl.cex = 1.5,
         cl.pos = "r",
         cl.ratio = 1)
dev.off()

##### correlation plot for lineages paper #####

jpeg("correlation_plot_targeted_pvalue_cor_filtered_large_thesis_corrected.jpeg", height = 1200, width = 2500)
corrplot(blocks_filtered$r, p.mat = blocks_filtered$P, sig.level = 0.05, insig = "blank", outline = TRUE, type = "upper",
         order = "hclust",  tl.col = "black", addgrid.col = NA, col = brewer.pal(n = 8, name = "RdBu"),
         mar = c(0.5, 1, 0.5, 1),
         tl.cex = 1.5,
         cl.cex = 1.5,
         cl.pos = "r",
         cl.ratio = 1)
dev.off()

# this one

jpeg("correlation_plot_targeted_pvalue_cor_filtered_large_thesis_corrected2.jpeg", height = 1600, width = 2800)
corrplot(blocks_filtered$r, p.mat = blocks_filtered$P, sig.level = 0.05, insig = "blank", outline = TRUE, type = "upper",
         tl.col = "black", addgrid.col = NA, col = brewer.pal(n = 8, name = "RdBu"),
         mar = c(0.5, 1, 0.5, 1),
         tl.cex = 1.5,
         cl.cex = 1.5,
         cl.pos = "r",
         cl.ratio = 1)
dev.off()

###### ergothioneine ######

data.df2 <- log_avg_ordered["Ergothioneine", ]

data.df2 <- rbind(snp_blocks, data.df2)
data.df2 <- data.df2[-4, ] # remove 4013 as strange values in this

results_cor <- rcorr(as.matrix(t(data.df2)))

# filter for p-value below 0.05

# results_cor_p <- as.data.frame(results_cor$P[116:135, 1:115]) # create data frame with metabolites as rows and blocks as columns
is_numeric <- is.na(as.numeric(rownames(results_cor$P))) # select all values and convert to T/F responses
is_trues <- which(is_numeric == TRUE) # select all numeric values
is_falses <- which(is_numeric == FALSE) # select all character values
results_cor_p <- as.data.frame(results_cor$P[(is_falses[length(is_falses)]+1):is_trues[length(is_trues)], 1:is_falses[length(is_falses)]]) # create data frame with metabolites as rows and blocks as columns
# results_cor_p <- as.data.frame(results_cor$P[116:151, 1:115]) # create data frame with metabolites as rows and blocks as columns
results_cor_p$metabolites <- "Ergothioneine" # add metabolites
results_cor_p <- as.data.frame(results_cor_p)
results_cor_p$snp_block <- rownames(results_cor_p)
# results_p <- melt(results_cor_p, "metabolites", "snp_block") # create list for filtering based on metabolites
results_p <- results_cor_p
results_p$value <- results_p$`results_cor$P[(is_falses[length(is_falses)] + 1):is_trues[length(is_trues)], 1:is_falses[length(is_falses)]]`
length(unique(results_p$snp_block)) # number of snps
results_p <- results_p[results_p$value < 0.05, ]
length(unique(results_p$snp_block)) # number of snps with significant correlation
pos_cor_blocks <- unique(results_p$snp_block)

data.df2 <- log_avg_ordered["Ergothioneine", ]

data.df2 <- rbind(snp_blocks[rownames(snp_blocks) %in% pos_cor_blocks, ], data.df2) # [, -15] to remove h37rv as not in snp data

results_cor2 <- t(data.df2)
results_circoplot <- rcorr(as.matrix(results_cor2))

# write.csv(results_circoplot$r, "results_circoplot_targeted_r_thesis.csv")
# write.csv(results_circoplot$P, "results_circoplot_targeted_p_thesis.csv")

# find the SNP block and correlated metabolite(s)

results_cor_p2 <- as.data.frame(results_circoplot$P)
results_cor_r2 <- as.data.frame(results_circoplot$r)
results_cor_p2$`340132`
results_cor_r2$`340132`
snpblock <- results_cor_r2

# select for only high correlation values

# results_cor_sig_r <- as.data.frame(results_circoplot$r[102:137,1:101])
# results_cor_sig_r <- as.data.frame(results_circoplot$r[92:111, 1:91])
is_numeric <- is.na(as.numeric(rownames(results_circoplot$r))) # select all values and convert to T/F responses
is_trues <- which(is_numeric == TRUE) # select all numeric values i.e. metabolites
is_falses <- which(is_numeric == FALSE) # select all character values i.e. nsSNP blocks
results_cor_sig_r <- as.data.frame(results_circoplot$r[(is_falses[length(is_falses)]+1):is_trues[length(is_trues)], 1:is_falses[length(is_falses)]]) # create data frame with metabolites as rows and blocks as columns
# results_cor_sig_r <- as.data.frame(results_circoplot$r[103:138, 1:102])
results_cor_sig_r$metabolites <- "Ergothioneine"
results_cor_sig_r <- as.data.frame(results_cor_sig_r)
results_cor_sig_r$snp_block <- rownames(results_cor_sig_r)
# results_r <- melt(results_cor_sig_r, "metabolites")
results_r <- results_cor_sig_r
length(unique(results_r$snp_block)) # number of snps
results_r$value <- results_r$`results_circoplot$r[(is_falses[length(is_falses)] + 1):is_trues[length(is_trues)], 1:is_falses[length(is_falses)]]`
results_r_pos50 <- results_r[results_r$value >= 0.5, ]
results_r_neg50 <- results_r[results_r$value <= -0.5, ]
results_r_50 <- rbind(results_r_pos50, results_r_neg50)
length(unique(results_r_50$snp_block)) # number of snps with significant correlation
pos_cor_sig_blocks <- unique(results_r_50$snp_block)

data.df2 <- log_avg_ordered["Ergothioneine", ]

data.df2 <- rbind(snp_blocks[rownames(snp_blocks) %in% pos_cor_sig_blocks, ], data.df2)
# data.df2 <- data.df2[-4, ] # remove 4013 as all the same value

blocks_filtered <- t(data.df2)
blocks_filtered <- rcorr(as.matrix(blocks_filtered))

jpeg("correlation_plot_targeted_pvalue_cor_filtered_large_thesis_erothioneine.jpeg", height = 600, width = 800, res = 100)
corrplot(blocks_filtered$r, p.mat = blocks_filtered$P, sig.level = 0.05, insig = "blank", outline = TRUE, type = "upper",
         tl.col = "black", addgrid.col = NA, col = brewer.pal(n = 8, name = "RdBu"),
         mar = c(0.5, 1, 0.5, 1),
         tl.cex = 1.5,
         cl.cex = 1.5,
         cl.pos = "r",
         cl.ratio = 1)
dev.off()

jpeg("correlation_plot_targeted_pvalue_cor_filtered_large_thesis_erothioneine.jpeg", height = 3750, width = 4250, res = 300)
corrplot(blocks_filtered$r, p.mat = blocks_filtered$P, sig.level = 0.05, insig = "blank", outline = TRUE, type = "upper",
         tl.col = "black", addgrid.col = NA, col = brewer.pal(n = 8, name = "RdBu"),
         mar = c(0.5, 1, 0.5, 1),
         tl.cex = 1.5,
         cl.cex = 1.5,
         cl.pos = "r",
         cl.ratio = 1)
dev.off()

##### correlation matrices of nsSNP blocks and untargeted metabolites #####

untargeted_mdf <- rbind(snp_blocks, nlog_avg2)

results_untargeted <- rcorr(as.matrix(t(untargeted_mdf)))

is_numeric <- is.na(as.numeric(rownames(results_untargeted$P))) # select all values and convert to T/F responses
is_trues <- which(is_numeric == TRUE) # select all numeric values
is_falses <- which(is_numeric == FALSE) # select all character values
results_untargeted_p <- as.data.frame(results_untargeted$P[(is_falses[length(is_falses)]+1):is_trues[length(is_trues)], 1:is_falses[length(is_falses)]]) # create data frame with metabolites as rows and blocks as columns
# results_cor_p <- as.data.frame(results_cor$P[116:151, 1:115]) # create data frame with metabolites as rows and blocks as columns
results_untargeted_p$metabolites <- rownames(results_untargeted_p) # add metabolites
results_untargeted_p <- as.data.frame(results_untargeted_p)
results_u <- melt(results_untargeted_p, "metabolites") # create list for filtering based on metabolites
length(unique(results_u$variable)) # number of snps
results_u <- results_u[results_u$value < 0.05, ]
length(unique(results_u$variable)) # number of snps with significant correlation
pos_untargeted_blocks <- unique(results_u$variable)

untargeted_mdf <- rbind(snp_blocks[rownames(snp_blocks) %in% pos_untargeted_blocks, ], untargeted_mdf) # [, -15] to remove h37rv as not in snp data

results_untargeted2 <- t(untargeted_mdf)
results_circoplot <- rcorr(as.matrix(results_untargeted2))

cor_df <- results_circoplot$r
cor_df <- as.data.frame(cor_df)
cor_df_p <- results_circoplot$P
cor_df_p <- as.data.frame(cor_df_p)

filter_correlation <- function(x){
  indexing <- which(x < 0.5 & x > -0.5)
  y <- as.data.frame(x[indexing])
  y$blocks <- rownames(x[indexing])
  print(y)
}

lapply(cor_df, filter_correlation)

# find the SNP block and correlated metabolite(s)

results_cor_p2 <- as.data.frame(results_circoplot$P)
results_cor_r2 <- as.data.frame(results_circoplot$r)
results_cor_p2$`340132`
results_cor_r2$`340132`
snpblock <- results_cor_r2

# select for only high correlation values

# results_cor_sig_r <- as.data.frame(results_circoplot$r[102:137,1:101])
# results_cor_sig_r <- as.data.frame(results_circoplot$r[92:111, 1:91])
is_numeric <- is.na(as.numeric(rownames(results_circoplot$r))) # select all values and convert to T/F responses
is_trues <- which(is_numeric == TRUE) # select all numeric values i.e. metabolites
is_falses <- which(is_numeric == FALSE) # select all character values i.e. nsSNP blocks
results_cor_sig_r <- as.data.frame(results_circoplot$r[(is_falses[length(is_falses)]+1):is_trues[length(is_trues)], 1:is_falses[length(is_falses)]]) # create data frame with metabolites as rows and blocks as columns
# results_cor_sig_r <- as.data.frame(results_circoplot$r[103:138, 1:102])
results_cor_sig_r$metabolites <- rownames(results_cor_sig_r)
results_cor_sig_r <- as.data.frame(results_cor_sig_r)
results_r <- melt(results_cor_sig_r, "metabolites")
length(unique(results_r$variable)) # number of snps
results_r_pos50 <- results_r[results_r$value >= 0.75, ]
results_r_neg50 <- results_r[results_r$value <= -0.75, ]
results_r_50 <- rbind(results_r_pos50, results_r_neg50)
length(unique(results_r_50$variable)) # number of snps with significant correlation
pos_cor_sig_blocks <- unique(results_r_50$variable)

# data.df2 <- log_avg_ordered[rownames(log_avg_ordered) %in% rownames(sig_res_targeted), ]

untargeted_mdf <- rbind(snp_blocks[rownames(snp_blocks) %in% pos_cor_sig_blocks, ], untargeted_mdf)

blocks_filtered <- t(untargeted_mdf)
blocks_filtered <- rcorr(as.matrix(blocks_filtered))

blocks_filtered$r <- blocks_filtered$r[-4, -4]
blocks_filtered$P <- blocks_filtered$P[-4, -4]

blocks_filtered$r[which(blocks_filtered$r[, 9] > 0.74), ]
blocks_filtered$r[which(blocks_filtered$r[, 9] < -0.74), ]

##### Untargeted and targeted metabolites DAA #####

##### POS vs NEG FC differences

# targeted <- read.csv("Targeted_all_strains_exp3_July20.csv")
# targeted <- targeted[!targeted$X == 333, ]
# targeted <- targeted[!targeted$X == 367, ]
# targeted[sapply(targeted, is.factor)] <- lapply(targeted[sapply(targeted, is.factor)], function(x) as.numeric(as.character(x)))
# rownames(targeted) <- make.names(targeted$X, unique = TRUE)
# metab <- targeted
# metab <- metab[!rownames(metab) == "X649.5", ]
# metab[, 2:42] <- sapply(metab[, 2:42], as.numeric)

# metab <- metab[!rownames(metab) == "X281.4", ]
# metab <- metab[!rownames(metab) == "X372.4", ]
# metab <- metab[!rownames(metab) == "X639.3", ]
# metab <- metab[!rownames(metab) == "H37Rv.5", ]
# metab <- metab[!rownames(metab) == "X119.5", ]

# log2 transform
# metab <- log(metab[, 2:42], 2)
# metab[metab == "-Inf"] <- 0
# metab <- as.data.frame(t(metab))

daa_getfit <- function(metabolite_abundance){
# differential expression of metabolite levels
  names <- c("x119", "x119", "x119", "x119", "x119", "x119",
             "x173", "x173", "x173", "x173", "x173", "x173",
             "x212", "x212", "x212", "x212", "x212", "x212",
             "x232", "x232", "x232", "x232", "x232", "x232",
             "x281", "x281", "x281", "x281", "x281", "x281",
             "x293", "x293", "x293", "x293", "x293", "x293",
             "x318", "x318", "x318", "x318", "x318", "x318",
             "x345", "x345", "x345", "x345", "x345", "x345",
             "x346", "x346", "x346", "x346", "x346", "x346",
             "x372", "x372", "x372", "x372", "x372", "x372",
             "x374", "x374", "x374", "x374", "x374", "x374",
             "x440", "x440", "x440", "x440", "x440", "x440",
             "x639", "x639", "x639", "x639", "x639", "x639",
             "x649", "x649", "x649", "x649", "x649",
             "h37rv", "h37rv", "h37rv", "h37rv", "h37rv", "h37rv")

metabolite <- log(metabolite_abundance, 2)
metabolite[metabolite == "-Inf"] <- 0
# only look at comparing against h37rv
design <- model.matrix(~0 + factor(names))
colnames(design) <- c("h37rv", "x119", "x173", "x212", "x232", "x281", "x293", "x318", "x345", "x346", "x372", "x374", "x440", "x639", "x649")
fit <- lmFit(metabolite, design) #fit a linear model
contrast.matrix <- makeContrasts(x119-h37rv, x173-h37rv, x212-h37rv, x232-h37rv, x281-h37rv, x293-h37rv, x318-h37rv, x345-h37rv, x346-h37rv, x372-h37rv, x374-h37rv, x440-h37rv, x639-h37rv, x649-h37rv,
                                 levels=design)
fit.cb <- eBayes(contrasts.fit(fit, contrast.matrix))
return(fit.cb)
}
# remove logging for targeted metabolites
daa_getfit2 <- function(metabolite_abundance){
  # differential expression of metabolite levels
  names <- c("x119", "x119", "x119", "x119", "x119",
             "x173", "x173", "x173", "x173", "x173", "x173",
             "x212", "x212", "x212", "x212", "x212", "x212",
             "x232", "x232", "x232", "x232", "x232", "x232",
             "x281", "x281", "x281", "x281", "x281",
             "x293", "x293", "x293", "x293", "x293", "x293",
             "x318", "x318", "x318", "x318", "x318", "x318",
             "x345", "x345", "x345", "x345", "x345", "x345",
             "x346", "x346", "x346", "x346", "x346", "x346",
             "x372", "x372", "x372", "x372", "x372",
             "x374", "x374", "x374", "x374", "x374", "x374",
             "x440", "x440", "x440", "x440", "x440", "x440",
             "x639", "x639", "x639", "x639", "x639",
             "x649", "x649", "x649", "x649", "x649",
             "h37rv", "h37rv", "h37rv", "h37rv", "h37rv")
  
  # only look at comparing against h37rv
  design <- model.matrix(~0 + factor(names))
  colnames(design) <- c("h37rv", "x119", "x173", "x212", "x232", "x281", "x293", "x318", "x345", "x346", "x372", "x374", "x440", "x639", "x649")
  fit <- lmFit(metabolite_abundance, design) #fit a linear model
  contrast.matrix <- makeContrasts(x119-h37rv, x173-h37rv, x212-h37rv, x232-h37rv, x281-h37rv, x293-h37rv, x318-h37rv, x345-h37rv, x346-h37rv, x372-h37rv, x374-h37rv, x440-h37rv, x639-h37rv, x649-h37rv,
                                   levels=design)
  fit.cb <- eBayes(contrasts.fit(fit, contrast.matrix))
  return(fit.cb)
}

metab_pos2 <- metab_pos[, -90]

# get fit
fit_p <- daa_getfit(metab_pos2)
fit_n <- daa_getfit(metab_neg)
# fit_targeted <- daa_getfit2(metab)

# get results
res_p <- topTable(fit_p, sort="none", n=Inf) # gets all results
res_n <- topTable(fit_n, sort="none", n=Inf)
# res_targeted <- topTable(fit_targeted, sort="none", n=Inf)

# filter for significance
filter_pvalue <- function(x){
  res_pvalues <- x[x$adj.P.Val < 0.01, ]
  res_pvalues <- res_pvalues[order(res_pvalues$adj.P.Val, decreasing = FALSE), ]
  return(res_pvalues)
}

sig_res_pos <- filter_pvalue(res_p)
sig_res_neg <- filter_pvalue(res_n)
# sig_res_targeted <- filter_pvalue(res_targeted)

# select differentially expressed genes for each comparison individually and subset original data frame by these
significant_filter <- function(fit_data, x, dge_res){
  y <- topTable(fit_data, coef = x, sort="none", n=Inf) # gets significant x119 - h37rv results
  y <- y[y$adj.P.Val < 0.05, ]
  y[, 7] <- "119"
  y1 <- topTable(fit_data, coef = x+1, sort="none", n=Inf)
  y1 <- y1[y1$adj.P.Val < 0.05, ]
  y1[, 7] <- "173"
  y2 <- topTable(fit_data, coef = x+2, sort="none", n=Inf)
  y2 <- y2[y2$adj.P.Val < 0.05, ]
  y2[, 7] <- "212"
  y3 <- topTable(fit_data, coef = x+3, sort="none", n=Inf)
  y3 <- y3[y3$adj.P.Val < 0.05, ]
  y3[, 7] <- "232"
  y4 <- topTable(fit_data, coef = x+4, sort="none", n=Inf)
  y4 <- y4[y4$adj.P.Val < 0.05, ]
  y4[, 7] <- "281"
  y5 <- topTable(fit_data, coef = x+5, sort="none", n=Inf)
  y5 <- y5[y5$adj.P.Val < 0.05, ]
  y5[, 7] <- "293"
  y6 <- topTable(fit_data, coef = x+6, sort="none", n=Inf)
  y6 <- y6[y6$adj.P.Val < 0.05, ]
  y6[, 7] <- "318"
  y7 <- topTable(fit_data, coef = x+7, sort="none", n=Inf)
  y7 <- y7[y7$adj.P.Val < 0.05, ]
  y7[, 7] <- "345"
  y8 <- topTable(fit_data, coef = x+8, sort="none", n=Inf)
  y8 <- y8[y8$adj.P.Val < 0.05, ]
  y8[, 7] <- "346"
  y9 <- topTable(fit_data, coef = x+9, sort="none", n=Inf)
  y9 <- y9[y9$adj.P.Val < 0.05, ]
  y9[, 7] <- "372"
  y10 <- topTable(fit_data, coef = x+10, sort="none", n=Inf)
  y10 <- y10[y10$adj.P.Val < 0.05, ]
  y10[, 7] <- "374"
  y11 <- topTable(fit_data, coef = x+11, sort="none", n=Inf)
  y11 <- y11[y11$adj.P.Val < 0.05, ]
  y11[, 7] <- "440"
  y12 <- topTable(fit_data, coef = x+12, sort="none", n=Inf)
  y12 <- y12[y12$adj.P.Val < 0.05, ]
  y12[, 7] <- "639"
  y13 <- topTable(fit_data, coef = x+13, sort="none", n=Inf)
  y13 <- y13[y13$adj.P.Val < 0.05, ]
  y13[, 7] <- "649"
  y_final <- list(y, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13) # groups everything
  fil1 <- rownames(y_final[[1]]) # get the metabolite number which is significant
  fil2 <- rownames(y_final[[2]])
  fil3 <- rownames(y_final[[3]])
  fil4 <- rownames(y_final[[4]])
  fil5 <- rownames(y_final[[5]])
  fil6 <- rownames(y_final[[6]])
  fil7 <- rownames(y_final[[7]])
  fil8 <- rownames(y_final[[8]])
  fil9 <- rownames(y_final[[9]])
  fil10 <- rownames(y_final[[10]])
  fil11 <- rownames(y_final[[11]])
  fil12 <- rownames(y_final[[12]])
  fil13 <- rownames(y_final[[13]])
  fil14 <- rownames(y_final[[14]])
  filter_by <- c(fil1, fil2, fil3, fil4, fil5, fil6, fil7, fil8, fil9, fil10, fil11, fil12, fil13, fil14) # groups metabolite numbers
  print(length(filter_by))
  filter_by <- unique(filter_by)
  print(length(filter_by))
  return(dge_res[rownames(dge_res) %in% filter_by, ])
}

# get the significant metabolites for the individual strains
significant_filter_metabolites <- function(fit_data, x, dge_res){
  y <- topTable(fit_data, coef = x, sort="none", n=Inf) # gets significant x119 - h37rv results
  y <- y[y$adj.P.Val < 0.05, ]
  y[, 7] <- "119"
  y1 <- topTable(fit_data, coef = x+1, sort="none", n=Inf)
  y1 <- y1[y1$adj.P.Val < 0.05, ]
  y1[, 7] <- "173"
  y2 <- topTable(fit_data, coef = x+2, sort="none", n=Inf)
  y2 <- y2[y2$adj.P.Val < 0.05, ]
  y2[, 7] <- "212"
  y3 <- topTable(fit_data, coef = x+3, sort="none", n=Inf)
  y3 <- y3[y3$adj.P.Val < 0.05, ]
  y3[, 7] <- "232"
  y4 <- topTable(fit_data, coef = x+4, sort="none", n=Inf)
  y4 <- y4[y4$adj.P.Val < 0.05, ]
  y4[, 7] <- "281"
  y5 <- topTable(fit_data, coef = x+5, sort="none", n=Inf)
  y5 <- y5[y5$adj.P.Val < 0.05, ]
  y5[, 7] <- "293"
  y6 <- topTable(fit_data, coef = x+6, sort="none", n=Inf)
  y6 <- y6[y6$adj.P.Val < 0.05, ]
  y6[, 7] <- "318"
  y7 <- topTable(fit_data, coef = x+7, sort="none", n=Inf)
  y7 <- y7[y7$adj.P.Val < 0.05, ]
  y7[, 7] <- "345"
  y8 <- topTable(fit_data, coef = x+8, sort="none", n=Inf)
  y8 <- y8[y8$adj.P.Val < 0.05, ]
  y8[, 7] <- "346"
  y9 <- topTable(fit_data, coef = x+9, sort="none", n=Inf)
  y9 <- y9[y9$adj.P.Val < 0.05, ]
  y9[, 7] <- "372"
  y10 <- topTable(fit_data, coef = x+10, sort="none", n=Inf)
  y10 <- y10[y10$adj.P.Val < 0.05, ]
  y10[, 7] <- "374"
  y11 <- topTable(fit_data, coef = x+11, sort="none", n=Inf)
  y11 <- y11[y11$adj.P.Val < 0.05, ]
  y11[, 7] <- "440"
  y12 <- topTable(fit_data, coef = x+12, sort="none", n=Inf)
  y12 <- y12[y12$adj.P.Val < 0.05, ]
  y12[, 7] <- "639"
  y13 <- topTable(fit_data, coef = x+13, sort="none", n=Inf)
  y13 <- y13[y13$adj.P.Val < 0.05, ]
  y13[, 7] <- "649"
  y[, 8] <- rownames(y) # get the metabolite number which is significant
  y1[, 8] <- rownames(y1)
  y2[, 8] <- rownames(y2)
  y3[, 8] <- rownames(y3)
  y4[, 8] <- rownames(y4)
  y5[, 8] <- rownames(y5)
  y6[, 8] <- rownames(y6)
  y7[, 8] <- rownames(y7)
  y8[, 8] <- rownames(y8)
  y9[, 8] <- rownames(y9)
  y10[, 8] <- rownames(y10)
  y11[, 8] <- rownames(y11)
  y12[, 8] <- rownames(y12)
  y13[, 8] <- rownames(y13)
  y_final <- rbind(y3, y4, y8, y9, y, y2, y10, y7, y13, y1, y5, y6, y11, y12) # groups everything
  print(y_final)
}

# sig_res_targeted <- significant_filter(fit_data=fit_targeted, x=1, dge_res=res_targeted)
sig_res_pos <- significant_filter(fit_data=fit_p, x=1, dge_res=res_p)
sig_res_neg <- significant_filter(fit_data=fit_n, x=1, dge_res=res_n)
# sig_res_metabs <- significant_filter_metabolites(fit_data=fit_p, x=1, dge_res=res_p)

index <- which(c(sig_res_metabs$logFC > 0) == TRUE)
sig_res_metabs2 <- sig_res_metabs
sig_res_metabs2$logFC[index] <- "POS"
index <- which(c(sig_res_metabs$logFC < 0) == TRUE)
sig_res_metabs2$logFC[index] <- "NEG"
colnames(sig_res_metabs2) <- c(colnames(sig_res_metabs2)[1:6], "strain", "metabolite")

sig_res_metabs3 <- pivot_wider(sig_res_metabs2[, c(1, 7:8)], names_from = strain, values_from = logFC)

##### Correlation plots #####

snps2 <- snps_nonsynonymous
rownames(snps2) <- gsub("s", "", rownames(snps2))

# select blocks in all the combinations
snps_noncorrelated <- unique(snps2)

nlog_avg2 <- nlog_avg[rownames(nlog_avg) %in% rownames(sig_res_neg), ]
plog_avg2 <- plog_avg[rownames(plog_avg) %in% rownames(sig_res_pos), ]
nlog_avg2 <- nlog_avg2[,1:14]
plog_avg2 <- plog_avg2[,1:14]
colnames(nlog_avg2) <- gsub("X", "s", colnames(nlog_avg2))
colnames(plog_avg2) <- gsub("X", "s", colnames(plog_avg2))

data.df <- rbind(snp_blocks, nlog_avg2)
data.df <- as.data.frame(t(data.df))
cor_matrix <- cor(data.df[,unlist(lapply(data.df, is.numeric))])
corrplot(cor_matrix[1:100,1:100], type="upper", method="color")

# log_avg_ordered <- log_avg2[c(1:9, 14:36, 38:39, 10:13, 37), ]
log_avg_ordered <- log_avg2
colnames(log_avg_ordered) <- c(colnames(snp_blocks))

data.df <- rbind(snp_blocks, log_avg_ordered) # remove h37rv as not in snp data
cor_matrix <- cor(t(data.df))

# create plot function
correlation_plot <- function(df, tl_cex, cl_cex){
  corrplot(df, outline = TRUE, tl.col = "black", addgrid.col = NA, col = brewer.pal(n = 8, name = "RdBu"),
           mar = c(0.5, 1, 0.5, 1),
           tl.cex = tl_cex,
           cl.cex = cl_cex,
           cl.pos = "r",
           cl.ratio = 1)
}

correlation_plot_pvalues <- function(df, pvalues_df, sig, method, tl_ex, cl_cex){
  corrplot(df, p.mat = pvalues_df, sig.level = sig, type = method, utline = TRUE, tl.col = "black", addgrid.col = NA, col = brewer.pal(n = 8, name = "RdBu"),
           mar = c(0.5, 1, 0.5, 1),
           tl.cex = tl_cex,
           cl.cex = cl_cex,
           cl.pos = "r",
           cl.ratio = 1)
}

# plot
correlation_plot(cor_matrix[115:155,1:100], 1, 1)

# in significantly abundant metabolites only (against H37Rv ref) 
data.df2 <- log_avg_ordered[rownames(log_avg_ordered) %in% rownames(sig_res_targeted), ]
data.df2 <- rbind(snp_blocks, data.df2) # [, -15] to remove h37rv as not in snp data

# remove 4013 as all NA values
data.df2 <- data.df2[-4, ]

sig_metab_circoplot <- t(data.df2)
cor_matrix <- cor(sig_metab_circoplot)

# show p-values
results_cor <- rcorr(as.matrix(sig_metab_circoplot))

# results_cor <- rcorr(as.matrix(sig_metab_circoplot))

# filter for p-value above 0.05
# results_cor_p <- as.data.frame(results_cor$P[117:136, 1:116])
# results_cor_p <- as.data.frame(results_cor$P[116:153, 1:115])
is_numeric <- is.na(as.numeric(rownames(results_cor$P))) # select all values and convert to T/F responses
is_trues <- which(is_numeric == TRUE) # select all numeric values i.e. metabolites
is_falses <- which(is_numeric == FALSE) # select all character values i.e. nsSNP blocks
results_cor_p <- as.data.frame(results_cor$P[(is_falses[length(is_falses)]+1):is_trues[length(is_trues)], 1:is_falses[length(is_falses)]]) # create data frame with metabolites as rows and blocks as columns
results_cor_p$metabolites <- rownames(results_cor_p)
results_cor_p <- as.data.frame(results_cor_p)
results_p <- melt(results_cor_p, "metabolites")
length(unique(results_p$variable)) # number of snps
results_p <- results_p[results_p$value < 0.05, ]
length(unique(results_p$variable)) # number of snps with significant correlation
pos_cor_blocks <- unique(results_p$variable)

data.df2 <- log_avg_ordered[rownames(log_avg_ordered) %in% rownames(sig_res_targeted), ]

# snp data frame plus metabolite abundance
data.df2 <- rbind(snp_blocks[rownames(snp_blocks) %in% pos_cor_blocks, ], data.df2) # [, -15] to remove h37rv as not in snp data

results_cor2 <- t(data.df2)
results_circoplot <- rcorr(as.matrix(results_cor2))
# jpeg("correlation_plot_targeted_pvalue_filtered_thesis.jpeg", height = 1000, width = 3500)
# corrplot(results_circoplot$r[96:115,1:94], p.mat = results_circoplot$P[96:115,1:94], sig.level = 0.05, outline = TRUE, tl.col = "black", addgrid.col = NA, col = brewer.pal(n = 8, name = "RdBu"),
#          mar = c(0.5, 1, 0.5, 1),
#          tl.cex = 1.5,
#          cl.cex = 1.5,
#          cl.pos = "r",
#          cl.ratio = 1)
# dev.off()

# select for only high correlation values

# results_cor_sig_r <- as.data.frame(results_circoplot$r[96:115,1:94])
# results_cor_sig_r <- as.data.frame(results_circoplot$r[105:142,1:104])
is_numeric <- is.na(as.numeric(rownames(results_circoplot$r))) # select all values and convert to T/F responses
is_trues <- which(is_numeric == TRUE) # select all numeric values i.e. metabolites
is_falses <- which(is_numeric == FALSE) # select all character values i.e. nsSNP blocks
results_cor_sig_r <- as.data.frame(results_circoplot$r[(is_falses[length(is_falses)]+1):is_trues[length(is_trues)], 1:is_falses[length(is_falses)]]) # create data frame with metabolites as rows and blocks as columns
results_cor_sig_r$metabolites <- rownames(results_cor_sig_r)
results_cor_sig_r <- as.data.frame(results_cor_sig_r)
results_r <- melt(results_cor_sig_r, "metabolites")
length(unique(results_r$variable)) # number of snps
results_r_pos50 <- results_r[results_r$value >= 0.75, ]
results_r_neg50 <- results_r[results_r$value <= -0.75, ]
results_r_50 <- rbind(results_r_pos50, results_r_neg50)
length(unique(results_r_50$variable)) # number of snps with significant correlation
pos_cor_sig_blocks <- unique(results_r_50$variable)

positive_correlations <- log_avg_ordered[rownames(log_avg_ordered) %in% rownames(sig_res_targeted), ]

positive_correlations <- rbind(snps_noncorrelated[rownames(snps_noncorrelated) %in% pos_cor_sig_blocks, ], positive_correlations) # [, -15] to remove h37rv as not in snp data

blocks_filtered <- t(positive_correlations)
blocks_filtered <- rcorr(as.matrix(blocks_filtered))

# change names

colnames(blocks_filtered$r)[colnames(blocks_filtered$r) == "N2.Acetyl.L.ornithine"] <- "Acetyl.ornithine"
colnames(blocks_filtered$r)[colnames(blocks_filtered$r) == "Nicotinamide.adenine.dinucleotide..NAD."] <- "Nicotinamide.adenine.dinucleotide"
colnames(blocks_filtered$r)[colnames(blocks_filtered$r) == "L.2.Aminoadipic.acid"] <- "Aminoadipic.acid"
colnames(blocks_filtered$r)[colnames(blocks_filtered$r) == "D.Alanyl.D.alanine"] <- "Alanyl.alanine"

colnames(blocks_filtered$P)[colnames(blocks_filtered$P) == "N2.Acetyl.L.ornithine"] <- "Acetyl.ornithine"
colnames(blocks_filtered$P)[colnames(blocks_filtered$P) == "Nicotinamide.adenine.dinucleotide..NAD."] <- "Nicotinamide.adenine.dinucleotide"
colnames(blocks_filtered$P)[colnames(blocks_filtered$P) == "L.2.Aminoadipic.acid"] <- "Aminoadipic.acid"
colnames(blocks_filtered$P)[colnames(blocks_filtered$P) == "D.Alanyl.D.alanine"] <- "Alanyl.alanine"

colnames(blocks_filtered$n)[colnames(blocks_filtered$n) == "N2.Acetyl.L.ornithine"] <- "Acetyl.ornithine"
colnames(blocks_filtered$n)[colnames(blocks_filtered$n) == "Nicotinamide.adenine.dinucleotide..NAD."] <- "Nicotinamide.adenine.dinucleotide"
colnames(blocks_filtered$n)[colnames(blocks_filtered$n) == "L.2.Aminoadipic.acid"] <- "Aminoadipic.acid"
colnames(blocks_filtered$n)[colnames(blocks_filtered$n) == "D.Alanyl.D.alanine"] <- "Alanyl.alanine"

# jpeg("correlation_plot_targeted_pvalue_filtered_repeat_thesis.jpeg", height = 1000, width = 3500)
corrplot(blocks_filtered$r[43:80,1:42], p.mat = blocks_filtered$P[43:80,1:42], sig.level = 0.05, outline = TRUE, tl.col = "black", addgrid.col = NA, col = brewer.pal(n = 8, name = "RdBu"),
         mar = c(0.5, 1, 0.5, 1),
         tl.cex = 1.5,
         cl.cex = 1.5,
         cl.pos = "r",
         cl.ratio = 1)
# dev.off()

# for lineages paper

jpeg("correlationplot_targeted_pvalue_filtered_repeat_thesis.jpeg", height = 6500, width = 9000, res = 300)
corrplot(blocks_filtered$r, p.mat = blocks_filtered$P, sig.level = 0.05, insig = "blank", outline = TRUE, type = "upper",
         tl.col = "black", addgrid.col = NA, col = brewer.pal(n = 8, name = "RdBu"),
         mar = c(0.5, 1, 0.5, 1),
         tl.cex = 1.5,
         cl.cex = 1.5,
         cl.pos = "r",
         cl.ratio = 0.5)
dev.off()

# without legend

jpeg("correlationplot_targeted_pvalue_filtered_repeat_thesis.jpeg", height = 3500, width = 5500, res = 300)
corrplot(blocks_filtered$r, p.mat = blocks_filtered$P, sig.level = 0.05, insig = "blank", outline = TRUE, type = "upper",
         tl.col = "black", addgrid.col = NA, col = brewer.pal(n = 8, name = "RdBu"),
         mar = c(0.5, 1, 0.5, 1),
         tl.cex = 1.5,
         cl.cex = 1.5,
         cl.pos = "n",
         cl.ratio = 1)
dev.off()

jpeg("correlationplot_targeted_pvalue_filtered_repeat_thesis2.jpeg", height = 3750, width = 7000, res = 150)
corrplot(blocks_filtered$r, p.mat = blocks_filtered$P, sig.level = 0.05, insig = "blank", outline = TRUE, type = "upper",
         order = "hclust",  tl.col = "black", addgrid.col = NA, col = brewer.pal(n = 8, name = "RdBu"),
         mar = c(0.5, 1, 0.5, 1),
         tl.cex = 1.5,
         cl.cex = 1.5,
         cl.pos = "r",
         cl.ratio = 1)
dev.off()

write.csv(blocks_filtered$r, "blocks_filtered_targeted.csv")

##### mixomics targeted metabolites plus others #####

X <- list(SNP = t(snp_blocks),
          metab_neg = t(nlog_avg2),
          metab_pos = t(plog_avg2))
Y <- factor(c("lin2", "lin4", "lin2", "lin1", "lin1", "lin4", "lin4", "lin2", "lin1", "lin1", "lin2", "lin4", "lin4", "lin2"))
list.keepX <- list(SNP = c(3,3), metab_pos = c(150,150), metab_neg = c(150,150))
res.diablo <- block.splsda(X, Y, keepX=list.keepX)
vars_negpos <- selectVar(res.diablo) # check variables selected
plotIndiv(res.diablo)
plotVar(res.diablo, var.names=c(FALSE, FALSE), legend=TRUE, style="lattice") # plot=FALSE gives matrix
plotDiablo(res.diablo, ncomp = 1)

log_avg2 <- log_avg[, 1:14]
colnames(log_avg2) <- colnames(snps_nonsynonymous)
X <- list(SNP = t(snps_noncorrelated[, 1:14]),
          metab_neg = t(log_avg2))
Y <- factor(c("lin2", "lin4", "lin2", "lin1", "lin1", "lin4", "lin4", "lin2", "lin1", "lin1", "lin2", "lin4", "lin4", "lin2"))
list.keepX <- list(SNP = c(116,116), metab_neg = c(40,40))
res.diablo <- block.splsda(X, Y, keepX=list.keepX)
blah <- selectVar(res.diablo) # check variables selected
plotIndiv(res.diablo)
plotVar(res.diablo, var.names=c(FALSE, FALSE), legend=TRUE, style="lattice") # plot=FALSE gives matrix
plotDiablo(res.diablo, ncomp = 1)
plotLoadings(res.diablo, comp = 1, contrib = "max")
network(res.diablo, blocks = c(1,2),
        color.node = c('darkorchid', 'brown1'), 
        cutoff = 0.75, save = 'jpeg', name.save = 'DIABLOnetwork')

png("diablo_circplot_targeted_neg_blocks_thesis.png", height = 1000, width = 1000, res = 150)
circ <- circosPlot(res.diablo, cutoff=0.5, showIntraLinks=FALSE, comp=1, size.variables=0.7, size.legend = 1.29, line = FALSE,
                   color.cor = c("#3388BB", "black"),
                   color.blocks = c("#CC6666", "gray67"))
dev.off()

png("diablo_circplot_targeted_neg_component2_thesis.png", height = 750, width = 750)
circ <- circosPlot(res.diablo, cutoff=0.75, showIntraLinks=FALSE, comp=2, size.variables=0.7, size.legend = 1.29, line = FALSE,
                   color.cor = c("#3388BB", "black"),
                   color.blocks = c("#CC6666", "gray67"))
dev.off()

# get metabolites that correlate with SNP blocks
# metab_rt_pos and metab_rt_neg are the objects with all metabolites and retention times included
circ <- as.data.frame(circ)
correlated <- circ$`6112` < -0.85
cor_085 <- which(correlated == TRUE)
correlated <- circ$`6112` > 0.85
cor_085 <- c(cor_085, which(correlated == TRUE))
correlated_085 <- circ[cor_085, ]
metab_correlated_neg <- metab_rt_neg[rownames(metab_rt_neg) %in% rownames(correlated_085), ]
metab_correlated_neg <- metab_correlated_neg[, c(1, 91)]

circ <- circosPlot(res.diablo, cutoff=0.85, showIntraLinks=FALSE, comp=2, size.variables=0.7, size.legend = 1.29)

circ <- as.data.frame(circ)
correlated <- circ$`1849` < -0.85
cor_085 <- which(correlated == TRUE)
correlated <- circ$`1849` > 0.85
cor_085 <- c(cor_085, which(correlated == TRUE))
correlated_085_comp2 <- circ[cor_085, ]
metab_correlated_neg_comp2 <- metab_rt_neg[rownames(metab_rt_neg) %in% rownames(correlated_085_comp2), ]
metab_correlated_neg_comp2 <- metab_correlated_neg_comp2[, c(1, 91)]

X <- list(SNP = t(selected_noncorrelated),
          metab_pos = t(plog_avg2))
Y <- factor(c("lin2", "lin4", "lin2", "lin1", "lin1", "lin4", "lin4", "lin2", "lin1", "lin1", "lin2", "lin4", "lin4", "lin2"))

list.keepX <- list(SNP = c(3,3), metab_pos = c(150,7))
res.diablo <- block.splsda(X, Y, keepX=list.keepX)

circ <- circosPlot(res.diablo, cutoff=0.85, showIntraLinks=FALSE, comp=1, size.variables=0.7, size.legend = 1.29)

circ <- as.data.frame(circ)
correlated <- circ$`6112` < -0.85
cor_085 <- which(correlated == TRUE)
correlated <- circ$`6112` > 0.85
cor_085 <- c(cor_085, which(correlated == TRUE))
correlated_085 <- circ[cor_085, ]
metab_correlated_pos <- metab_rt_pos[rownames(metab_rt_pos) %in% rownames(correlated_085), ]
metab_correlated_pos <- metab_correlated_pos[, c(1, 91)]

circ <- circosPlot(res.diablo, cutoff=0.85, showIntraLinks=FALSE, comp=2, size.variables=0.7, size.legend = 1.29)

circ <- as.data.frame(circ)
correlated <- circ$`1849` < -0.85
cor_085 <- which(correlated == TRUE)
correlated <- circ$`1849` > 0.85
cor_085 <- c(cor_085, which(correlated == TRUE))
correlated_085_comp2 <- circ[cor_085, ]
metab_correlated_pos_comp2 <- metab_rt_pos[rownames(metab_rt_pos) %in% rownames(correlated_085_comp2), ]
metab_correlated_pos_comp2 <- metab_correlated_pos_comp2[, c(1, 91)]

metab_correlated <- cbind.fill(metab_correlated_neg, metab_correlated_neg_comp2, metab_correlated_pos, metab_correlated_pos_comp2, fill = NA)
colnames(metab_correlated) <- c("rt_neg", "mass_neg", "rt_neg2", "mass_neg2", "rt_pos", "mass_pos", "rt_pos2", "mass_pos2")
metab_correlated$rt_neg <- gsub(".*@", "", metab_correlated$rt_neg)
metab_correlated$rt_neg2 <- gsub(".*@", "", metab_correlated$rt_neg2)
metab_correlated$rt_pos <- gsub(".*@", "", metab_correlated$rt_pos)
metab_correlated$rt_pos2 <- gsub(".*@", "", metab_correlated$rt_pos2)
# write.csv(metab_correlated, "metabolites_correlated_snps_diablo_26062020.csv")

# differentially expressed metabolites
dem <- c("X464.1376", "X155.0399", "X364.0603", "X313.005", "X471.1584", "X430.1321", "X412.0978", "X406.9771", "X328.9998", "X170.0147")
metab_correlated_dem <- metab_rt_neg[rownames(metab_rt_neg) %in% dem, c(1, 91)]
metab_correlated_dem$rt <- gsub(".*@", "", metab_correlated_dem$Compound)
# write.csv(metab_correlated_dem, "metabolites_correlated_dem_26062020.csv")

##### Figures for Paper #####

library(ggpubr)

setwd("C:/Users/ashle/OneDrive/Documents/PhD/Lineages Paper/Figures")

sum3$lin <- c("lin4", "lin2", "lin1", "lin4", "lin4", "lin2", "lin1", "lin1", "lin2", "lin4", "lin2", "lin1", "lin4", "lin2")
sum3$Order <- c(13, 6, 3, 12, 14, 8, 4, 2, 7, 10, 5, 1, 11, 9)

sum01$lin <- c("lin2", "lin4", "lin2", "lin1", "lin1", "lin4", "lin4", "lin2", "lin1", "lin1", "lin2", "lin4", "lin4", "lin2")
sum01$Order <- c(5, 13, 6, 3, 1, 12, 14, 8, 4, 2, 7, 10, 11, 9)

# added
sum_all <- cbind(sum01, sum3)
sum_all <- sum_all[, c(1:2, 5:7)]
colnames(sum_all) <- c("sum_all", "Strains", "sum_non", "sum_syn", "lin")
#

sum_all$Order <- c("5", "13", "6", "3", "1", "12", "14", "8", "4", "2", "7", "10", "11", "9")
sum_all <- sum_all[order(sum_all$Order), ]
sum_all2 <- melt(sum_all[, -1], id.vars = "Strains")
sum_all2$Order <- rep(sum_all$Order, 2)
sum_all2 <- sum_all2[1:28, ]
# changed
sum_all2 <- sum_all2[c(1:14, 29:42), ]
sum_all2$value <- as.numeric(sum_all2$value)
ggplot(sum_all2, aes(fill = variable, x = reorder(Strains, Order), y = value)) + geom_bar(stat = "identity")

png("number_allsnps_andnonsynonymous_thesis.png", width = 11500, height = 10000, res = 300)

ggplot(sum_all2, aes(fill = variable, x = reorder(Strains, as.numeric(Order)), y = value)) +
  geom_bar(stat="identity") + xlab("Strains") +
  ylab("Number of SNPs") +
  # geom_text(aes(label=sum), vjust=1.6, color="white", size=3.5) +
  scale_fill_manual(labels = c("All SNPs", "Non-synonymous SNPs"), "Lineages", values=c("gray67", "#CC6666", "#003366")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size=90, margin=margin(1, 1, 1, 1)),
        axis.title.y = element_text(size=90, margin=margin(1, 7, 1, 1)),
        axis.text.x = element_text(size=70, angle=-45, hjust=0.1, vjust=0.5),
        axis.text.y = element_text(size=70),
        legend.title = element_text(size=90),
        legend.text = element_text(size=65),
        axis.line = element_line(colour="black")) +
  scale_y_continuous(expand = c(0.005, 10), limits = c(0, 2250))

dev.off()

ggplot(sum_all2, aes(fill = variable, x = reorder(Strains, as.numeric(Order)), y = value)) +
  geom_bar(stat="identity") + xlab("Strains") +
  ylab("Number of SNPs") +
  # geom_text(aes(label=sum), vjust=1.6, color="white", size=3.5) +
  scale_fill_manual(labels = c("All SNPs", "Non-synonymous SNPs"), "Lineages", values=c("gray67", "#CC6666", "#003366")) +
  theme_minimal() +
  scale_y_continuous(expand = c(0.005, 10), limits = c(0, 2250))

colors = c("blue", "red", "#e391c0")
# colors = c("gray47", "#5566AA", "#F17B4B")
fit_hc <- hclust(d, method="ward")
clus3 = cutree(fit_hc, 3)

png("dendrogram_nonsynonymous_snps_thesis.png", height = 10000, width = 5000, res = 300)
plot(as.phylo(fit_hc), cex = 7, edge.width = 3, label.offset = 0.5, tip.color = colors[clus3])
dev.off()

# for functional catergories for individual strains
png("functional_categories_strains_3_thesis.png", height = 20000, width = 20000, res = 300)
ggplot(snp_cats, aes(fill=pathways_superfamily, x = reorder(as.character(strain), order_strain) , y = as.numeric(counts)), legend=FALSE) + 
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_fill_manual(values = getPalette(colourCount)) +
  ylab("Functional categories (% total)") +
  xlab("Lineage") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title.x = element_text(size=100, margin=margin(25, 1, 1, 1)),
        axis.title.y = element_text(size=100, margin=margin(1, 25, 1, 1)),
        axis.text.x = element_text(size=100, angle=-45, hjust=0.1, vjust=0.75, margin=margin(10, 1, 1, 1)),
        axis.text.y = element_text(size=100, margin=margin(1, 10, 1, 1)),
        axis.line = element_line(colour="black", size = 3),
        panel.grid.major = element_line(size = 3),
        panel.grid.minor = element_line(size = 2.5)) +
        coord_flip() +
        guides(fill=FALSE)
dev.off()

png("functional_categories_strains_legend_2_thesis.png", height = 10000, width = 20000, res = 300)
ggplot(snp_cats, aes(fill=pathways_superfamily, x = reorder(as.character(strain), order_strain) , y = as.numeric(counts)), legend=FALSE) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = getPalette(colourCount)) +
  ylab("Categories (% total)") +
  xlab("Lineage") +
  theme_minimal() +
  theme(legend.position = "left",
        axis.title.x = element_text(size=40, margin=margin(1, 1, 1, 1)),
        axis.title.y = element_text(size=40, margin=margin(1, 1, 1, 1)),
        axis.text.x = element_text(size=38, angle=-45, hjust=0.1, vjust=0.75, margin=margin(1, 1, 1, 1)),
        axis.text.y = element_text(size=38, margin=margin(1, 1, 1, 1)),
        legend.title = element_text(size=80),
        legend.text = element_text(size=80),
        axis.line = element_line(colour="black", size = 1),
        legend.key.size = unit(0.7, "cm"),
        legend.direction = "vertical") +
        coord_flip()
dev.off()

png("Doubling_times_paper.png", height = 4750, width = 4700, res = 300)
ggplot(double, aes(x = reorder(Strain, Order), y = Hour, color = as.character(Lineage))) + geom_point(pch = "-", size = 9) +
  geom_errorbar(aes(ymin=Hour-SD, ymax=Hour+SD), width=.35, size = 3) +
  xlab("Strains") +
  ylab("Doubling time (hours)") +
  scale_color_manual("Lineages", values=c("#CC6666", "gray67", "#003366")) +
  theme_minimal() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.title.x = element_text(size=37, margin=margin(-17, 1, 1, 1)),
        axis.title.y = element_text(size=37, margin=margin(1, 7, 1, 1)),
        axis.text.x = element_text(size=37, angle=-45, hjust=0.1, vjust=0.5),
        axis.text.y = element_text(size=37),
        legend.title = element_text(size=37),
        legend.text = element_text(size=37),
        axis.line = element_line(colour="black", size = 1),
        panel.grid.major = element_line(size = 1.75),
        panel.grid.minor = element_line(size = 1.5)) +
  scale_y_continuous(limits = c(0, 58))
dev.off()

# pca plot

names_targeted <- rownames(sig_res_targeted) # get names of sig metabolites

log_avg <- sapply(c("X119", "X173", "X212", "X232", "X281", "X293", "X318", "X345", "X346", "X372", "X374", "X440", "X639", "X649", "H37Rv"), function(x) colMeans(metab[grep(x, rownames(metab)), ]))
log_avg <- as.data.frame(log_avg)
# log_avg2 <- log_avg[, -15]
sub <- log_avg[! rownames(log_avg) == "γ.Aminobutryic.acid",] # remove outlier metabolite
sub <- sub[! rownames(sub) == "GABA",] # remove outlier metabolite

# sub <- sub[c(1:9, 14:36, 38:39, 10:13, 37), ]
colnames(sub) <- gsub("X", "s", colnames(sub))

sub <- as.data.frame(sub[rownames(sub) %in% names_targeted, ]) # subset original dataframe logfc values by significant metabolites
data <- sub
data <- as.data.frame(t(data))
# add groups
groups <- c("#999999", "#003366", "#999999", "#CC6666", "#CC6666", "#003366", "#003366", "#999999", "#CC6666", "#CC6666", "#999999", "#003366", "#003366", "#999999", "#003366")
data$group <- groups
# create pca object
pca <- prcomp(data[, 1:37])

png("pca_sig_targeted_thesis_repeat_2023.png", width = 5000, height = 5000, res = 300)
autoplot(pca, data=data, label=TRUE, shape=FALSE, label.size=12,
  colour=c("#999999", "#003366", "#999999", "#CC6666", "#CC6666", "#003366", "#003366", "#999999", "#CC6666", "#CC6666", "#999999", "#003366", "#003366", "#999999", "#003366")) +
  theme_bw() +
  # stat_ellipse(data = pca$group, aes(group = group)) +
  theme(plot.margin = margin(2, 2, 2, 2, "cm"),
        axis.title.x = element_text(size=25, margin=margin(7, 1, 1, 1)),
        axis.title.y = element_text(size=25, margin=margin(1, 7, 1, 1)),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=25),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(size = 1.5),
        panel.grid.minor = element_line(size = 1.25))
dev.off()

colours_lin <- c("#999999", "#003366", "#999999", "#CC6666", "#CC6666", "#003366", "#003366", "#999999", "#CC6666", "#CC6666", "#999999", "#003366", "#003366", "#999999", "#003366")

# all metabolites

log_avg <- sapply(c("X119", "X173", "X212", "X232", "X281", "X293", "X318", "X345", "X346", "X372", "X374", "X440", "X639", "X649", "H37Rv"), function(x) colMeans(metab[grep(x, rownames(metab)), ]))
log_avg <- as.data.frame(log_avg)
sub <- log_avg[! rownames(log_avg) == "γ.Aminobutryic.acid",] # remove outlier metabolite
sub <- sub[! rownames(sub) == "GABA",] # remove outlier metabolite

sub_all <- sub[c(1:9, 14:36, 38:39, 10:13, 37), ]
colnames(sub_all) <- gsub("X", "s", colnames(sub_all))

data <- sub_all
data <- as.data.frame(t(data))
# add groups
groups <- c("#999999", "#003366", "#999999", "#CC6666", "#CC6666", "#003366", "#003366", "#999999", "#CC6666", "#CC6666", "#999999", "#003366", "#003366", "#999999", "#003366")
data$group <- groups
pca <- prcomp(data[, 1:39])

png("pca_targeted_thesis_repeat_2023.png", width = 5000, height = 5000, res = 300)
autoplot(pca, data=data, label=TRUE, shape=FALSE, label.size=12,
  colour=c("#999999", "#003366", "#999999", "#CC6666", "#CC6666", "#003366", "#003366", "#999999", "#CC6666", "#CC6666", "#999999", "#003366", "#003366", "#999999", "#003366")) +
  theme_bw() +
  # stat_ellipse(data = pca$group, aes(group = group)) +
  theme(plot.margin = margin(2, 2, 2, 2, "cm"),
        axis.title.x = element_text(size=35, margin=margin(7, 1, 1, 1)),
        axis.title.y = element_text(size=35, margin=margin(1, 7, 1, 1)),
        axis.text.x = element_text(size=35),
        axis.text.y = element_text(size=35),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(size = 1.5),
        panel.grid.minor = element_line(size = 1.25))
dev.off()

##### different colour palette #####

# for functional catergories for individual strains
png("functional_categories_strains_3_thesis_colourscheme2.png", height = 20000, width = 20000, res = 300)
ggplot(snp_cats, aes(fill=pathways_superfamily, x = reorder(as.character(strain), order_strain) , y = as.numeric(counts)), legend=FALSE) + 
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_fill_viridis(discrete=TRUE) +
  # scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(colourCount)) +
  ylab("Functional categories (% total)") +
  xlab("Lineage") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title.x = element_text(size=100, margin=margin(25, 1, 1, 1)),
        axis.title.y = element_text(size=100, margin=margin(1, 25, 1, 1)),
        axis.text.x = element_text(size=100, angle=-45, hjust=0.1, vjust=0.75, margin=margin(10, 1, 1, 1)),
        axis.text.y = element_text(size=100, margin=margin(1, 10, 1, 1)),
        axis.line = element_line(colour="black", size = 3),
        panel.grid.major = element_line(size = 3),
        panel.grid.minor = element_line(size = 2.5)) +
  coord_flip() +
  guides(fill=FALSE)
dev.off()

png("functional_categories_strains_3_thesis_colourscheme3.png", height = 20000, width = 20000, res = 300)
ggplot(snp_cats, aes(fill=pathways_superfamily, x = reorder(as.character(strain), order_strain) , y = as.numeric(counts)), legend=FALSE) + 
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_fill_manual(values=c("#30123BFF", "#3A2E7CFF", "#4249B1FF", "#4662D7FF", "#477AF2FF", "#4392FEFF", "#36AAF9FF", "#25C0E7FF", "#1AD5CEFF", "#1AE4B6FF", "#2DF09DFF",
                    "#4EF97DFF", "#72FE5EFF", "#95FE45FF", "#AFFA37FF", "#C7EF34FF", "#DDE037FF", "#EFCE3AFF", "#FABA39FF", "#FEA331FF", "#FC8725FF", "#F66B19FF",
                     "#EB510EFF", "#DD3D08FF", "#CB2A04FF", "#B41B01FF", "#990E01FF", "#7A0403FF")) +
  # scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(colourCount)) +
  ylab("Functional categories (% total)") +
  xlab("Lineage") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title.x = element_text(size=100, margin=margin(25, 1, 1, 1)),
        axis.title.y = element_text(size=100, margin=margin(1, 25, 1, 1)),
        axis.text.x = element_text(size=100, angle=-45, hjust=0.1, vjust=0.75, margin=margin(10, 1, 1, 1)),
        axis.text.y = element_text(size=100, margin=margin(1, 10, 1, 1)),
        axis.line = element_line(colour="black", size = 3),
        panel.grid.major = element_line(size = 3),
        panel.grid.minor = element_line(size = 2.5)) +
  coord_flip() +
  guides(fill=FALSE)
dev.off()

png("functional_categories_strains_legend_2_thesis_colour3.png", height = 10000, width = 20000, res = 300)
ggplot(snp_cats, aes(fill=pathways_superfamily, x = reorder(as.character(strain), order_strain) , y = as.numeric(counts)), legend=FALSE) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values=c("#30123BFF", "#3A2E7CFF", "#4249B1FF", "#4662D7FF", "#477AF2FF", "#4392FEFF", "#36AAF9FF", "#25C0E7FF", "#1AD5CEFF", "#1AE4B6FF", "#2DF09DFF",
                             "#4EF97DFF", "#72FE5EFF", "#95FE45FF", "#AFFA37FF", "#C7EF34FF", "#DDE037FF", "#EFCE3AFF", "#FABA39FF", "#FEA331FF", "#FC8725FF", "#F66B19FF",
                             "#EB510EFF", "#DD3D08FF", "#CB2A04FF", "#B41B01FF", "#990E01FF", "#7A0403FF")) +
  ylab("Categories (% total)") +
  xlab("Lineage") +
  theme_minimal() +
  theme(legend.position = "left",
        axis.title.x = element_text(size=40, margin=margin(1, 1, 1, 1)),
        axis.title.y = element_text(size=40, margin=margin(1, 1, 1, 1)),
        axis.text.x = element_text(size=38, angle=-45, hjust=0.1, vjust=0.75, margin=margin(1, 1, 1, 1)),
        axis.text.y = element_text(size=38, margin=margin(1, 1, 1, 1)),
        legend.title = element_text(size=80),
        legend.text = element_text(size=80),
        axis.line = element_line(colour="black", size = 1),
        legend.key.size = unit(0.7, "cm"),
        legend.direction = "vertical") +
  coord_flip()
dev.off()

##### global plot #####

# loads all files in current directory ending with .csv
temp = list.files(pattern="*.vcf")
for (i in 1:length(temp)) assign(temp[i], read.table(temp[i]))

files <- Sys.glob("*.vcf")
# files <- files[1:15]

lineage_list <- vector(mode="list", length=length(files))
names(lineage_list) <- files
for (filename in files) {
  # read data
  sample <- read.csv(filename, sep = "\t")
  lineage_list[[filename]] <- sample %>% dplyr::select(POS, ALT, REF, QUAL) %>%
    mutate(lin=filename)
  
}

new_lin <- bind_rows(lineage_list)

global_lins <- new_lin
global_lins$lin <- gsub(".fastq.sam.bam.sorted.bam.bcf.variants.vcf.final_variants.vcf", "", global_lins$lin)
global_lins$lin <- gsub(".fastq.sam.bam.sorted.bam.sorted.bam.bcf.variants.vcf.final_variants.vcf", "", global_lins$lin)

###

names <- grep('.vcf$', ls(), value=T)

newDF <- mapply(get, grep('.vcf', ls(), value=T))

# get annotation information

allglobal <- global_lins

gl_metadat <- read.csv("Global_dataset_meta.csv")
global <- read.csv("Global_dataset_meta.csv")

global$Filename.Accession.number

# filter for high quality score
allglobal <- subset(allglobal, allglobal$QUAL > 50)

# clean up SNP annotation based on first annotation
ref_dataset <- allglobal

# keep only single nucleotide variants
ref_dataset <- subset(ref_dataset, ref_dataset$REF %in% c("A", "G", "T", "C"))
ref_dataset <- subset(ref_dataset, ref_dataset$ALT %in% c("A", "G", "T", "C"))

# synonymous snps

df_ordered <- ref_dataset
df_ordered$Strain <- df_ordered$lin

df_lins <- data.frame(matrix(ncol = length(summary(factor(df_ordered$Strain))), nrow = dim(df_ordered)[1]))
colnames(df_lins) <- unique(df_ordered$Strain)
global_ordered <- cbind(df_ordered, df_lins)

colnames(global_ordered)[7:101] <- unique(df_ordered$Strain)

mine <- global_ordered %>% distinct(POS, .keep_all = TRUE)
mine <- mine$POS
#'%ni%' <- Negate('%in%')

# select SNP from all variants
global_ordered <- subset(global_ordered, global_ordered$REF %in% c("A", "G", "T", "C"))
global_ordered <- subset(global_ordered, global_ordered$ALT %in% c("A", "G", "T", "C"))

empty_table <- list()
for (position in mine){ #for snp position
  zero_one <- list()
  my_table <- global_ordered[global_ordered$POS == position, ] #grab all rows with this position
  lineages_present <- my_table$Strain
  print(lineages_present)
  # lineages_present <- paste("s", lineages_present, sep="") #grab all the lineages where this snp is present (rename to fit with names)
  lineages_absent <- setdiff(x = colnames(global_ordered)[7:101], y = lineages_present) #grab all lineages where this snp is absent
  zeros <- rep(0, length(lineages_absent)) #for every lineage where its absent add a zero
  print(zeros)
  ones <- rep(1, length(lineages_present)) #for every lineage where its present add a one
  print(ones)
  zeros_ones <- c(zeros, ones)
  print(zeros_ones)
  names(zeros_ones) <- c(lineages_absent, lineages_present)
  zeros_ones <- zeros_ones[colnames(global_ordered)[7:101]] #reorders the columns to the order of linesgaesin the original table
  empty_table <- rbind(empty_table, zeros_ones) #row = snp, add each row one at a time
}
rownames(empty_table) <- mine

global <- empty_table
global <- as.data.frame(global)

# add global snps
snps_add <- match(rownames(global), rownames(snps_all2)) # find the snps in both dataframes
snps_add <- snps_add[!is.na(snps_add)]
snps_add <- snps_all2[snps_add,] # subset snps by common positions

global2 <- global
global2 <- global2[rownames(global2) %in% rownames(snps_add), ]
global2 <- cbind(global2, snps_add)

# get snps that are in global dataset only
snps_uncommon <- subset(global, !rownames(global) %in% rownames(snps_all2))
snps_uncommon[, 96:109] <- 0
colnames(snps_uncommon)[96:109] <- colnames(snps_all2)

# add snps that are in snps_all2 dataset only
snps_uncommon_again <- subset(snps_all2, !rownames(snps_all2) %in% rownames(global))
snps_uncommon_again[, 15:109] <- 0
colnames(snps_uncommon_again)[15:109] <- colnames(global)

# combine all
global3 <- rbind(global2, snps_uncommon, snps_uncommon_again)

colnames(global3) <- gsub("_sequence", "", colnames(global3))
colnames(global3) <- gsub("_1$", "", colnames(global3))

# find metadata for the strains
gldf <- gl_metadat
gldf$Filename.Accession.number <- gsub("_sequence.fastq", "", gldf$Filename.Accession.number)
gldf$Filename.Accession.number <- gsub("_1/2.fastq", "", gldf$Filename.Accession.number)
gldf$Filename.Accession.number <- gsub(".fastq", "", gldf$Filename.Accession.number)
gldf$Filename.Accession.number <- gsub("MTB_V", "s", gldf$Filename.Accession.number)
gldf$Filename.Accession.number <- gsub("MTB_1_V", "s", gldf$Filename.Accession.number)
gldf$Filename.Accession.number <- gsub("MTB_2_V", "s", gldf$Filename.Accession.number)
gldf$Filename.Accession.number <- gsub("MTB_3_V", "s", gldf$Filename.Accession.number)
gldf$Filename.Accession.number <- gsub("MTB_4_V", "s", gldf$Filename.Accession.number)
gldf$Filename.Accession.number <- gsub("MTB_5_V", "s", gldf$Filename.Accession.number)
gldf$Filename.Accession.number <- gsub("MTB_6_V", "s", gldf$Filename.Accession.number)
gldf$Filename.Accession.number <- gsub("MTB_7_V", "s", gldf$Filename.Accession.number)
gldf$Filename.Accession.number <- gsub("MTB_GT_", "s", gldf$Filename.Accession.number)
gldf$Filename.Accession.number <- gsub("BJ", "", gldf$Filename.Accession.number)
gldf$Filename.Accession.number <- gsub("IO", "", gldf$Filename.Accession.number)
gldf$Filename.Accession.number <- gsub("EA", "", gldf$Filename.Accession.number)

gldf2 <- gldf[gldf$Filename.Accession.number %in% colnames(global3), ]
gldf3 <- global3[!colnames(global3) %in% gldf2$Filename.Accession.number, ]

# order global dataset by my dataset
index <- match(colnames(global3), gldf2$Filename.Accession.number)
gldf2 <- gldf2[index,]

# cluster by synonymous and all2 SNPs
snps_dend <- global3
snps_dend <- as.data.frame(t(snps_dend))
# rownames(snps_dend) <- gsub("X", "", rownames(snps_dend))
d <- dist(snps_dend, method = "euclidean") # distance matrix
fit_hc <- hclust(d, method="ward")

plot(fit_hc, hang = -1, 
     main = "Cluster dendrogram", sub = NULL,
     xlab = "Mtb Strain", ylab = "Height", cex=0.85,
     horiz = TRUE) # display dendogram

fit_hc <- as.dendrogram(fit_hc)

# colors = c("gray52", "#003366", "#CC5555", "red", "blue", "green")
colors = gldf2$Lineage
colors[colors == "Lineage 1"] <- "#e37f7f"
colors[colors == "Lineage 2"] <- "#0066CC"
colors[colors == "Lineage 3"] <- "#663399"
colors[colors == "Lineage 4"] <- "red"
colors[colors == "Lineage 5"] <- "#660000"
colors[colors == "Lineage 6"] <- "#339900"
colors[colors == "Lineage 7"] <- "#FFCC33"

# colors = c("gray47", "#5566AA", "#F17B4B")
fit_hc <- hclust(d, method="ward")
clus3 = cutree(fit_hc, 6)

fit_hc <- hclust(d, method="ward")

png("dendrogram_all_strains.png", height = 750, width = 900)
plot(as.phylo(fit_hc), cex = 1.5, label.offset = 0.5, tip.color = colors[clus3])
dev.off()

png("dendrogram_all_strains2.png", height = 750, width = 350)
plot(as.phylo(fit_hc), cex = 1.5, label.offset = 0.5, tip.color = colors[clus3])
dev.off()

png("dendrogram_all_strains3.png", height = 5750, width = 5750, res = 500)
plot(as.phylo(fit_hc), cex = 0.75, label.offset = 4.5, tip.color = colors, lab4ut="axial",
     type = "unrooted", no.margin = TRUE)
dev.off()

png("dendrogram_all_strains_circle.png", height = 3750, width = 3850, res = 300)
plot(as.phylo(fit_hc), cex = 1.25, label.offset = 4.5, tip.color = colors,
     type = "fan", no.margin = TRUE)
dev.off()


