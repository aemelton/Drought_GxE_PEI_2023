### AE Melton, 2022
# Use LEA package to run a Stucture type analysis
####### CHECK ONLINE TUTORIAL TO UPDATE

# 
setwd("FilePath")

# Load libraries
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("LEA")

library(LEA)
library(vcfR)

library(adegenet)
library(SNPRelate)

library(poppr)

# read in VCF file
vcf_file <- "merged.IndelGap.BiAllel.calls.vcf.gz" # An unprocessed VCF
vcf <- read.vcfR( vcf_file, verbose = FALSE )

# queryMETA(vcf, element = 'DP')
# queryMETA(vcf, element = 'GT')
vcf

# ***** Object of Class vcfR *****
# 50 samples
# 226022 CHROMs
# 1,982,691 variants
# Object size: 1547.6 Mb
# 0 percent missing data
# *****        *****         *****

# mask by read depth and write a cleaned vcf file
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
rownames(dp) <- 1:nrow(dp)
tail(dp)
heatmap.bp(dp[1:1000,], rlabels = FALSE)

# 
quants <- apply(dp, MARGIN=2, quantile, probs=c(0.05, 0.95), na.rm = TRUE)
dp2 <- sweep(dp, MARGIN=2, FUN = "-", quants[1,])
dp[dp2 < 0] <- NA

dp2 <- sweep(dp, MARGIN=2, FUN = "-", quants[2,])
dp[dp2 > 0] <- NA

dp[dp < 2] <- NA

vcf@gt[,-1][ is.na(dp) == TRUE ] <- NA
vcf

# ***** Object of Class vcfR *****
# 50 samples
# 226022 CHROMs
# 1,982,691 variants
# Object size: 1493.1 Mb
# 79.57 percent missing data
# *****        *****         **********
heatmap.bp(dp[1:1000,], rlabels = FALSE)

# Omit samples
myMiss <- apply(dp, MARGIN = 2, function(x){ sum( is.na(x) ) } )
myMiss <- myMiss / nrow(dp)
# length(vcf@gt[, c(TRUE, myMiss < 0.9)])

vcf@gt <- as.matrix(vcf@gt[, c(TRUE, myMiss < 0.9)])
vcf

# ***** Object of Class vcfR *****
# 41 samples
# 226022 CHROMs
# 1,982,691 variants
# Object size: 1356.3 Mb
# 76.02 percent missing data
# *****        *****         *****

dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
heatmap.bp(dp[1:1000,], rlabels = FALSE)

# Omit variants
myMiss <- apply(dp, MARGIN = 1, function(x){ sum( is.na(x) ) } )
myMiss <- myMiss / ncol(dp)
vcf <- vcf[myMiss < 0.5, ]
vcf

# ***** Object of Class vcfR *****
# 41 samples
# 24441 CHROMs
# 346,919 variants
# Object size: 396.8 Mb
# 31.62 percent missing data
# *****        *****         *****

dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
heatmap.bp(dp[1:1000,], rlabels = FALSE)

write.vcf(x = vcf, file = "cleaned.vcf.gz")

# assess SNPs for linkage disequilibrium
snpgdsVCF2GDS(vcf.fn = "cleaned.vcf.gz", out.fn = "ccm.gds") # vcf.fn = vcf_file
genofile <- openfn.gds("ccm.gds")
snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.2, autosome.only = F)
snpset.id <- unlist(unname(snpset))
head(snpset.id)
write.csv(x = snpset.id, file = "snpset_LD_checked.csv", row.names = F) # Save the snps you want to put in PCA. Can use these to filter later if needed so you dont have to do the whoel LD thing again.

vcf.trim <- vcf[snpset.id,]
vcf.trim

# ***** Object of Class vcfR *****
# 41 samples
# 24289 CHROMs
# 107,777 variants
# Object size: 201.8 Mb
# 31.64 percent missing data
# *****        *****         *****

write.vcf(x = vcf.trim, file = "cleaned_subset.vcf.gz", APPEND = F) # Save the final cleaned and trimmed VCF file

# Convert VCF to geno
system("gzip -dk cleaned_subset.vcf.gz")
vcf2geno(input.file = "cleaned_subset.vcf", output.file = "SNPs.geno", force = T)

# Calculate Fst
sample_names <- gsub("_merged.bam_sorted.bam","", colnames(vcf.trim@gt)[-1]) # Already set correctly. Maybe check first?
sample_names <- gsub("_sorted.bam","", sample_names) # Already set correctly. Maybe check first?
pop_code <- as.factor(c(rep("IDT3", 16), rep("UTT2", 20), rep("IDT3", 5)))
sample.id <- sample_names
cbind(sample.id, pop_code) # Visually inspect the pop codes

# Two populations: IDT3 and UTT2
flag <- pop_code %in% c("IDT3", "UTT2")
samp.sel <- sample.id[flag]
pop.sel <- pop_code[flag]
snpgdsVCF2GDS(vcf.fn = "cleaned_subset.vcf.gz", out.fn = "ccmSUB.gds",  method = "biallelic.only") # vcf.fn = vcf_file
genofile <- snpgdsOpen(filename = "ccmSUB.gds")
#snpset.id <- read.csv("snpset_LD_checked.csv")
v <- snpgdsFst(gdsobj = genofile, sample.id = colnames(vcf.trim@gt)[-1], population = pop.sel,
               method="W&C84", autosome.only = F)
v
v$Fst
v$MeanFst
summary(v$FstSNP)

# > v$Fst
# [1] 0.1079868
# > v$MeanFst
# [1] 0.0375043
# > summary(v$FstSNP)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.309216 -0.004266  0.005461  0.037504  0.046369  1.000000 

v$FstSNP[which(v$FstSNP < 0)] <- 0
summary(v$FstSNP)
sd(x = v$FstSNP, na.rm = T)

# > summary(v$FstSNP)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000000 0.005461 0.042192 0.046369 1.000000 
# > sd(x = v$FstSNP, na.rm = T)
# [1] 0.08795518
#

# # pca
# pc <- pca("SNPs.geno", scale = TRUE)
# summary(pc)
# 
# 
# # Tracy-Widom test
# tw <- tracy.widom(pc)
# tw$pvalues[1:5]
# 
# pdf("TracyWidom_Test.pdf")
# plot(tw$percentage)
# dev.off()

pca <- snpgdsPCA(gdsobj = genofile, autosome.only = F, remove.monosnp = T)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

pop_code <- as.factor(c(rep("IDT3", 16), rep("UTT2", 20), rep("IDT3", 5)))
sample.id <- sample_names
head(cbind(sample.id, pop_code))

tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

# seedlingID <- gsub(pattern = "_[A-z]*.[A-z]*.[A-z]*.", replacement = "", sample.id$V1)
x.label <- paste0("eigenvector 1 (", head(round(pc.percent, 2))[1], "%)")
y.label <- paste0("eigenvector 2 (", head(round(pc.percent, 2))[2], "%)")

pdf("Figures/2023_gds_PCA.pdf")
plot(y = tab$EV2, x = tab$EV1, col = c(rep("red", 16), rep("blue", 20), rep("red", 5)), ylab = y.label, xlab = x.label, pch = 19)
#text(y = tab$EV2, x = tab$EV1, labels = seedlingID, cex = 0.5)
abline(h = 0, v = 0, col = "light grey", lty = 2)
legend("bottomright", legend = levels(tab$pop), pch = 19, col = c("red", "blue"))
dev.off()
##### ##### #####

##### ##### #####
# snmf()
# main options
# K = number of ancestral populations
# entropy = TRUE: computes the cross-entropy criterion,
# CPU = 4 the number of CPUs.
project <- NULL
project <- snmf("SNPs.geno",
               K = 1:10,
               entropy = TRUE,
               repetitions = 10,
               iterations = 1000,
               project = "new")

# summary of the project
summary(project)

# plot cross-entropy criterion for all runs in the snmf project
pdf("Figures/CrossEntropy.pdf")
plot(project, col = "blue", pch = 19, cex = 1.2)
dev.off()

# select the best run for K = x
best <- which.min(cross.entropy(project, K = 2))
my.colors <- c("tomato", "lightblue")

# Need to re-organize this by likelihood scores. Need to loop over "best" run to get score for each sample and then sort.
barchart(project, K = 2, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix", 
         sort.by.Q = T) -> bp

# Sort the sample names
sample_names.sorted <- sample_names[match(bp$order, rownames(data.frame(sample_names)))]

pdf("Figures/StructureLikePlot.pdf")
barchart(project, K = 2, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix", 
         sort.by.Q = T) -> bp
axis(1, 
     at = 1:length(bp$order),
     labels = sample_names.sorted, 
     las = 1,
     cex.axis = .35,
     srt = 35, 
     las = 2)
dev.off()

# get the ancestral genotype frequency matrix, G, for the best run for best K. 
G.matrix <- G(project, K = 2, run = best)
Q.matrix <- Q(project, K = 2, run = best)
rownames(Q.matrix) <- sample_names
write.csv(Q.matrix, "Qmatrix.csv", row.names = T)

# Population differentiation tests
p <- snmf.pvalues(project,
                 entropy = TRUE,
                 ploidy = 2,
                 K = 2)
pvalues <- p$pvalues

pdf("Figures/p_values.pdf")
par(mfrow = c(2,1))
hist(pvalues, col = "orange")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .5)
dev.off()


### Read in the final cleaned 
vcf.file <- "cleaned_subset.vcf.gz"
vcf.for.snp <- read.vcfR( vcf.file, verbose = FALSE )
snp <- vcfR2genind(vcf.for.snp)

# You can access the data using the "@" sign:
snp
snp.genclone <- as.genclone(snp)
snp.genclone@pop <- pop.sel

### Calculate some popgen metrics
pairwise.fst.out <- hierfstat::pairwise.fst.dosage(dos = snp.genclone, pop = snp.genclone@pop, matching = F)
basicStats.out <- hierfstat::basic.stats(data = snp.genclone)
genet.dist.out <- hierfstat::genet.dist(dat = snp.genclone, diploid = T)

# Print out the results
pairwise.fst.out
basicStats.out
genet.dist.out

