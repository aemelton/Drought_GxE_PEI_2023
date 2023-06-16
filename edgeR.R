### AE Melton, 2021
# Based on http://www.nathalievialaneix.eu/doc/html/solution_edgeR-tomato-withcode.html

### Set working environment ###
#
getwd()
project.folder <- getwd()
#

#
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("edgeR")
# BiocManager::install("HTSFilter")
# BiocManager::install('mixOmics')
#

# Load the required libraries
#
library(edgeR)
library(limma)
library(RColorBrewer)
library(mixOmics)
library(HTSFilter)
library(statmod)
#

### Load Data ###
# The "idxstats" files come from the outputs of the 'samtools' idxstats function
# setwd("idxstats/")
# file.list <- list.files(pattern = ".tsv")
# tissues <- data.frame(nrow = 938307) # nrow = number of transcripts in reference transcriptome
# for (i in 1:length(file.list)) {
#   tmp <- read.csv(file.list[i], sep = "\t", header = F)[,c(3)]
#   tissues <- cbind(tissues, tmp)
# }
# SampleNames <- gsub(pattern = "_idxstats.tsv",replacement = "", x = file.list)
# #
# boop <- read.csv(file.list[1], sep = "\t", header = F)
# TranscriptIDs <- boop$V1
# #
# tissues <- tissues[,-1]
# colnames(tissues) <- c(SampleNames)
# rownames(tissues) <- TranscriptIDs
# head(tissues)
# 
# setwd(project.folder)
# write.csv(x = tissues, file = "idxstats_data.csv", row.names = T, col.names = T)

tissues <- read.csv(file = "idxstats_data_RemoveOutliers.csv", row.names = 1)
colnames(tissues) <- gsub(pattern = "X", replacement = "", x = colnames(tissues)) # Had an X in front of all the names
head(tissues)
nrow(tissues)
ncol(tissues)
#

#
# setwd(project.folder)
#

# Set up a table with treatment and population data for all samples
# For T1 vs T2
# treats <- read.csv("Treatments.csv", header = T)
# treats <- unique(treats)
# Samples <- colnames(tissues)
# seedlingID <- gsub(pattern = "_.*", replacement = "", x = Samples)
# #Treatment <- treats[which(treats$ID %in% seedlingID),3]
# 
# Treatment <- NULL
# for (i in 1:length(seedlingID)) {
#   boop <- treats[which(treats$ID == seedlingID[i]),]
#   Treatment <- rbind(Treatment, as.character(boop$Treatment))
# }
# 
# Pop <- NULL
# for (i in 1:length(seedlingID)) {
#   tmp <- treats[which(treats$ID == seedlingID[i]),]
#   Pop <- rbind(Pop, as.character(tmp$Pop))
# }
# 
# # Combine 
# Big.Data <- cbind(Pop, Treatment)
# head(Big.Data)
# 
# # For leaf vs root
# # Samples <- colnames(tissues)
# # Tissue <- gsub(pattern = "[[:digit:]]*_", replacement = "", x = Samples) # Just for leaf vs root
# # Groups <- data.frame(Samples, Tissue)
# #
# 
# # Assmbled big data frame
# # Group <- paste(gsub(pattern = "[[:digit:]]*_", replacement = "", x = Samples), Treatment, sep = "")
# Tissue <- gsub(pattern = "[[:digit:]]*_", replacement = "", x = Samples)
# Groups <- data.frame(Samples, Tissue, Big.Data[,1], Big.Data[,2]) # Currently need to remake "Groups" for each comparison
# colnames(Groups) <- c("Samples", "Tissue", "Pop", "Treatment")
# head(Groups)

# write.csv(x = Groups, file = "Sample_and_Group_Data.csv", row.names = F)
Groups <- read.csv("Sample_and_Group_Data_RemoveOutliers.csv")
colnames(Groups) <- c("Sample", "Tissue", "Population", "Treatment")
#

### Set up for analyses, explore data ###
#
setwd(project.folder)
Group <- factor(paste(Groups$Tissue, Groups$Population, Groups$Treatment, sep = "."))
beep <- cbind(Groups, Group)
head(beep)

dgeFull <- DGEList(counts = tissues, genes = rownames(tissues))
dgeFull

design.mat <- model.matrix(~0+Group) #model.matrix(~Groups$Tissue+Groups$Population+Groups$Treatment) # model.matrix(~ 0 + dgeFull$samples$group)
colnames(design.mat) <- levels(Group) 
rownames(design.mat) <- Groups$Sample
head(design.mat)

keep <- filterByExpr(y = dgeFull, design = design.mat)
table(keep)
dgeFilt <- dgeFull[keep, , keep.lib.sizes = FALSE]

dgeFilt <- calcNormFactors(object = dgeFilt, method = "TMM")
dgeFilt$samples

# pdf("Plots/mds.pdf")
# plotMDS(dgeFilt)
# dev.off()
#

# Make a more colorful plot to evaluate patterns and check for outliers
# Run once to find outliers, then re-run after outliers removed
# Groupy <- NULL
# for (i in 1:nrow(design.mat)) {
#   for (j in 1:ncol(design.mat)) {
#     if(design.mat[i,j] == 1){
#       Groupy[i] <- names(which(design.mat[i,] == 1))
#     }
#   }
# }
# 
# pdf("Plots/colorful_mds.pdf")
# plotMDS(dgeFilt, pch = 20, col = as.numeric(as.factor(Groupy))) #method = "bcv", 
# legend("bottomleft", 
#        as.character(unique(Groupy)), 
#        col = unique(as.numeric(as.factor(Groupy))), 
#        pch = 20,
#        cex = 0.65)
# dev.off()
# # Ok, we have an outlier.... 133_Leaves
# 
# km.res <- kmeans(dgeFilt$counts, centers = 2, nstart = 1234) # Using kvalue from previous function
# res.pca <- prcomp(t(dgeFilt$counts))
# 
# design.mat[order(res.pca$x[,1]),]
# rownames(res.pca$x)
#

# Re-do the pretty plot once the dissedent has been removed
Groupy <- NULL
for (i in 1:nrow(design.mat)) {
  for (j in 1:ncol(design.mat)) {
    if(design.mat[i,j] == 1){
      Groupy[i] <- names(which(design.mat[i,] == 1))
    }
  }
}

pdf("Plots/colorful_mds_RemoveOutliers.pdf")
plotMDS(dgeFilt, pch = 20, col = as.numeric(as.factor(Groupy))) #method = "bcv", 
legend("bottomleft", 
       as.character(unique(Groupy)), 
       col = unique(as.numeric(as.factor(Groupy))), 
       pch = 20,
       cex = 0.65)
dev.off()
#

#
d2 <- estimateDisp(dgeFilt, design.mat, robust = TRUE) 
pdf("Plots/bcv_RemoveOutliers.pdf")
plotBCV(d2)
dev.off()
#

### Define contrasts, perform GLM
#
my.contrasts <- makeContrasts(
  Leaves_IDT3_T = Leaves.IDT3.T2 - Leaves.IDT3.T1,
  Leaves_UTT2_T = Leaves.UTT2.T2 - Leaves.UTT2.T1,
  Roots_IDT3_T = Roots.IDT3.T2 - Roots.IDT3.T1,
  Roots_UTT2_T = Roots.UTT2.T2 - Roots.UTT2.T1,
  Leaves_POP_T1 = Leaves.UTT2.T1 - Leaves.IDT3.T1,
  Leaves_POP_T2 = Leaves.UTT2.T2 - Leaves.IDT3.T2,
  Roots_POP_T1 = Roots.UTT2.T1 - Roots.IDT3.T1,
  Roots_POP_T2 = Roots.UTT2.T2 - Roots.IDT3.T2,
  levels = design.mat)
#

#
# fit <- glmQLFit(d2, design.mat)
# save(fit, file = "qlfit_RemoveOutliers.RData")
qlf <- glmQLFTest(glmfit = fit, contrast = my.contrasts[, "Roots_POP_T2"])
save(qlf, file = "Roots_POP_T2_qlf_RemoveOutliers.RData")
topTags(qlf)
#

#
de2 <- decideTestsDGE(qlf, adjust.method = "BH", p.value = 0.05)
summary(decideTests(qlf))
de2tags <- rownames(d2)[as.logical(de2)]

pdf("Plots/Roots_POP_T2_Smear_RemoveOutliers.pdf")
plotSmear(qlf, de.tags = de2tags)
abline(h = c(-2, 2), col = "blue")
dev.off()
#

pdf("Plots/Roots_POP_T2_md_RemoveOutliers.pdf")
plotMD(qlf)
abline(h=c(-1, 1), col = "blue")
dev.off()
#

### Save plots ###
# plot an histogram of unadjusted p-values
pdf("Plots/Roots_POP_T2_value_hist_RemoveOutliers.pdf") ### WRITE FILE
hist(qlf$table[,"PValue"], breaks = 50)
dev.off()

# extract a summary of the differential expression statistics
resFilt <- topTags(qlf, n = nrow(qlf$table))
head(resFilt)

# resFilt <- topTags(dgeTestFilt, n=nrow(dgeTest$table))
# head(resFilt)

sum(resFilt$table$FDR < 0.05) # Changed from 0.01 in the example

# extract and sort differentially expressed genes
Sig <- resFilt$table[resFilt$table$FDR < 0.05,]

sigDownReg <- Sig[Sig$logFC < 0,]
#sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
head(sigDownReg)

sigUpReg <- Sig[Sig$logFC > 0,]
#sigUpReg <- sigDownReg[order(sigDownReg$logFC, decreasing=TRUE),]
head(sigUpReg)

write.csv(sigDownReg, file = "CSVs/Roots_POP_T2_sigDownReg_RemoveOutliers.csv") ### WRITE FILE
write.csv(sigUpReg, file = "CSVs/Roots_POP_T2_sigUpReg_RemoveOutliers.csv") ### WRITE FILE

# # IF USING TRANSCRIPTOME

# up <- read.csv("sigUpReg.csv")
# down <- read.csv("sigDownReg.csv")

trinotation <- read.csv("Trinotate.tsv", sep = "\t")
head(trinotation)

# UpRegTrinotate <- trinotation[trinotation$transcript_id %in% up$X,]
# DownRegTrinotate <- trinotation[trinotation$transcript_id %in% down$X,]
UpRegTrinotate <- trinotation[trinotation$transcript_id %in% rownames(sigUpReg),]
DownRegTrinotate <- trinotation[trinotation$transcript_id %in% rownames(sigDownReg),]

write.csv(x = UpRegTrinotate, file = "CSVs/Roots_POP_T2_UpRegTrinotate_RemoveOutliers.csv") ### WRITE FILE
write.csv(x = DownRegTrinotate, file = "CSVs/Roots_POP_T2_DownRegTrinotate_RemoveOutliers.csv") ### WRITE FILE
# 

# Volcano plot
volcanoData <- cbind(resFilt$table$logFC, -log10(resFilt$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
head(volcanoData)

pdf("Plots/Roots_POP_T2_VolcanoPlot_RemoveOutliers.pdf") ### WRITE FILE
plot(volcanoData, pch = 19)
dev.off()

# transform the normalized counts in log-counts-per-million
y <- cpm(y = dgeFilt, log = TRUE, prior.count = 1)
head(y)

# select 1% differentially expressed genes and produce a heatmap
zzz <-y[,which(Group %in% c("Roots.IDT3.T2", "Roots.UTT2.T2"))] # hOW CAN i MAKE THIS WORK JUST BY CHANGING CONTRAST AND MATCHING EVERYWHERE?
selY <- zzz[rownames(resFilt$table)[resFilt$table$PValue < 0.05 & 
                                    abs(resFilt$table$logFC) > 1.5],]
#selY <- selY[,which(Group %in% c("Roots.UTT2.T1", "Roots.IDT3.T1"))] # hOW CAN i MAKE THIS WORK JUST BY CHANGING CONTRAST AND MATCHING EVERYWHERE?
head(selY)
nrow(selY)

pdf("Plots/Roots_POP_T2_Heatmap_RemoveOutliers.pdf") ### WRITE FILE
cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)[255:1]
finalHM <- cim(t(selY), color = cimColor, symkey = FALSE) # error for L_POP_T1 & L_POP_T2
dev.off()

pdf("Plots/Roots_POP_T2_Dendrogram_RemoveOutliers.pdf") ### WRITE FILE
plot(finalHM$ddc, leaflab = "none")
abline(h = 10, lwd = 2, col = "pink")
dev.off()


pheatmap(t(selY),
         scale = "none", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         fontsize_row = 8,
         annotation_names_col = FALSE,
         gaps_col = c(3,6),
         display_numbers = TRUE,
         number_format = "%.2f",         
         height=12,
         width=6)
