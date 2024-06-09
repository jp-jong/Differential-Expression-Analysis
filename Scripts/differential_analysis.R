install.packages("BiocManager")
install.packages("tidyverse")
BiocManager::install(c("ballgown", "genefilter"))

library(ballgown)
library(genefilter)
library(tidyverse)

pheno_data <- (read.csv(file.choose()))
bg <- ballgown(dataDir="ballgown", samplePattern="test", pData=pheno_data)
bg_filtered <- subset(bg, "rowVars(texpr(bg))>=1")
results_transcripts <- stattest(bg_filtered, feature="transcript", covariate="condition", getFC=TRUE, meas="FPKM")
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filtered), geneIDs=ballgown::geneIDs(bg_filtered), transcriptIDs=ballgown::transcriptNames(bg_filtered), results_transcripts)
results_transcripts = arrange(results_transcripts, qval)
head(results_transcripts)
subset(results_transcripts, results_transcripts$qval<0.01)

# https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
de <- data.frame(results_transcripts["fc"], results_transcripts["pval"], results_transcripts["geneNames"])
de$diffexpressed <- "NO"
de$diffexpressed[log2(de$fc) > 1.0 & de$pval < 0.01] <- "UP"
de$diffexpressed[log2(de$fc) < -1.0 & de$pval < 0.01] <- "DOWN"
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
ggplot(data=de, aes(x=log2(fc), y=-log10(pval), col=diffexpressed)) +
  geom_point() + theme_minimal() +
  scale_colour_manual(values=mycolors) +
  geom_vline(xintercept=c(-1.0, 1.0), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")

dev.new()
de$diffexpressed <- "NO"
de$diffexpressed[log2(de$fc) > 1.0 & de$pval < 0.0001] <- "UP"
de$diffexpressed[log2(de$fc) < -1.0 & de$pval < 0.0001] <- "DOWN"
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$geneNames[de$diffexpressed != "NO"]
mycolors <- c("turquoise", "firebrick", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")
ggplot(data=de, aes(x=log2(fc), y=-log10(pval), col=diffexpressed, label=delabel)) +
  geom_point() + theme_minimal() + geom_text() +
  scale_colour_manual(values=mycolors) +
  geom_vline(xintercept=c(-1.0, 1.0), col="red") +
  geom_hline(yintercept=-log10(0.0001), col="red")

# down_regulated
library(dplyr)
down_diff <- de %>%
  filter(diffexpressed == "DOWN") %>%
  arrange(pval)

up_diff <- de %>%
  filter(diffexpressed == "UP") %>%
  arrange(pval)

# take the first 20 rows
head(down_diff,20)
head(up_diff,20)

# transpose

down_pval <- as.numeric(down_diff$pval)
names(down_pval) <- down_diff$geneNames
down_pval <- as.numeric(down_pval)

down_diff_names <- as.data.frame(t(down_diff$geneNames))
down_diff_pvalue <- as.data.frame(t(down_diff$pval))
down_diff_names <- as.numeric(t(down_diff$geneNames))
numeric_matrix <- as.matrix(as.numeric(down_diff_pvalue))
# down_trans <- rbind(down_diff_names, down_diff_pvalue)
names(down_diff_pvalue) <- down_diff_names
# convert to list
#down_trans_list <- as.list(down_trans)

BiocManager::install("topGO", force = TRUE)
BiocManager::install("org.Sc.sgd.db") # yeast database
library(BiocManager)
library(org.Sc.sgd.db)
library(topGO)
library(AnnotationDbi)
require(org.Sc.sgd.db)
require(topGO)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Sc.sgd.db")
BiocManager::install("org.Sc.sgd.db", force = TRUE)

# 
selection <- function(allScore){return(allScore<0.0001)} # function to return genes with p values less than 0.0001
allGO2genes <- annFUN.org(whichOnto = "BP", 
                          feasibleGenes = NULL, 
                          mapping="org.Sc.sgd.db", 
                          ID="entrez")
GOdata <- new ("topGOdata", 
               ontology = "BP",
               allGenes = down_pval,
               annot=annFUN.GO2genes,
               GO2genes = allGO2genes,
               geneSel=selection,
               nodeSize=10)
