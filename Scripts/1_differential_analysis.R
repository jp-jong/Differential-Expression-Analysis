install.packages("BiocManager")
install.packages("tidyverse")
BiocManager::install(c("ballgown", "genefilter"))

library(ballgown)
library(genefilter)
library(tidyverse)

# Use a csv file "phenodata.csv" available on NTU COOL.
pheno_data <- (read.csv(file.choose()))
bg <- ballgown(dataDir="ballgown", samplePattern="test", pData=pheno_data)

# https://github.com/CBC-UCONN/RNA-Seq-Model-Organism-Arabidopsis-thaliana#Seventh_Point_Header

# Select a subset where the variance is ≥1 to discard those with small changes.
bg_filtered <- subset(bg, "rowVars(texpr(bg))>=1")
# Conduct a statistical test for differential expressions between samples
results_transcripts <- stattest(bg_filtered, feature="transcript", covariate="condition", getFC=TRUE, meas="FPKM")
# Add identifiers to the results.
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filtered), geneIDs=ballgown::geneIDs(bg_filtered), transcriptIDs=ballgown::transcriptNames(bg_filtered), results_transcripts)
# Sort the result by q-values.
results_transcripts = arrange(results_transcripts, qval)
results_transcripts
# View a subset where q < 0.1.
subset(results_transcripts, results_transcripts$qval<0.1)


# VISUALIZATION
# Select transcripts with ≥1 mean expression level.
bg.expressed <- subset(bg, "rowMeans(texpr(bg))>=1")
# Take log2 of FPKM.
log2fpkm <- log2(texpr(bg.expressed, meas="FPKM") + 0.01)
# Change the format for ggplot.
log2fpkm.long <- as.data.frame(log2fpkm) %>% tidyr::gather(Dataset, FPKM)
# Draw a boxplot.
ggplot(log2fpkm.long, aes(y=FPKM, x=Dataset)) + geom_boxplot() + theme(axis.text.x = element_text(angle=90, hjust=0, vjust=.5)) + xlab("") + ylab("log2(FPKM+0.01)")
# Draw a violin plot.
dev.new()
ggplot(log2fpkm.long, aes(y=FPKM, x=Dataset)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + theme(axis.text.x = element_text(angle=90, hjust=0, vjust=.5)) + xlab("") + ylab("log2(FPKM+0.01)")

## DTD1 

log2fpkm["892", ]

# Adjust the format.
log2fpkm.892.long <- as.data.frame(log2fpkm["892", ]) %>% tidyr::gather(Dataset, FPKM)
log2fpkm.892.long <- data.frame(log2fpkm.892.long, Condition=c(rep("Ctrl", 3), rep("Exp", 3)))
dev.new()
# Draw a boxplot.
ggplot(log2fpkm.892.long, aes(x=Condition, y=FPKM, fill=Condition)) + geom_boxplot() + geom_point() + xlab("") + ylab("log2(FPKM+0.01)") + theme(legend.position = "bottom")

dev.new()
# Display one gene and its expression levels.
plotTranscripts(ballgown::geneIDs(bg)[892], bg, main=c('DTD1'), sample=c('test_ERR4626657','test_ERR4626658','test_ERR4626659','test_ERR4626660','test_ERR4626661','test_ERR4626662'))

## ADE8 

log2fpkm["1649", ]

# Adjust the format.
log2fpkm.1649.long <- as.data.frame(log2fpkm["1649", ]) %>% tidyr::gather(Dataset, FPKM)
log2fpkm.1649.long <- data.frame(log2fpkm.1649.long, Condition=c(rep("Ctrl", 3), rep("Exp", 3)))
dev.new()
# Draw a boxplot.
ggplot(log2fpkm.1649.long, aes(x=Condition, y=FPKM, fill=Condition)) + geom_boxplot() + geom_point() + xlab("") + ylab("log2(FPKM+0.01)") + theme(legend.position = "bottom")

dev.new()
# Display one gene and its expression levels.
plotTranscripts(ballgown::geneIDs(bg)[1649], bg, main=c('ADE8'), sample=c('test_ERR4626657','test_ERR4626658','test_ERR4626659','test_ERR4626660','test_ERR4626661','test_ERR4626662'))

# hsp26
log2fpkm["365", ]

# Adjust the format.
log2fpkm.365.long <- as.data.frame(log2fpkm["365", ]) %>% tidyr::gather(Dataset, FPKM)
log2fpkm.365.long <- data.frame(log2fpkm.365.long, Condition=c(rep("Ctrl", 3), rep("Exp", 3)))
dev.new()
# Draw a boxplot.
ggplot(log2fpkm.365.long, aes(x=Condition, y=FPKM, fill=Condition)) + geom_boxplot() + geom_point() + xlab("") + ylab("log2(FPKM+0.01)") + theme(legend.position = "bottom")

dev.new()
# Display one gene and its expression levels.
plotTranscripts(ballgown::geneIDs(bg)[365], bg, main=c('HSP26'), sample=c('test_ERR4626657','test_ERR4626658','test_ERR4626659','test_ERR4626660','test_ERR4626661','test_ERR4626662'))
