### S Buerki, 2021
# Edited by AEM for file paths

#######
#DE analyses
#######
# #Up regulated leaves
SigRegR <- read.csv("Leaves_POP_T2_DownRegTrinotate_RemoveOutliers.csv")

#Only select rows=transcripts with GO data
SigRegR <- SigRegR[-which(SigRegR$gene_ontology_BLASTX == "."),]
GO <- SigRegR$gene_ontology_BLASTX

GOs <- NULL
for(i in 1:length(GO)){
  int <- strsplit(as.vector(GO[i]), split = "GO:")[[1]][-1]
  GOtmp <- paste("GO:", sapply(strsplit(int, split="\\^"), "[[", 1), sep = "") 
  GOs <- c(GOs, GOtmp)
}
head(GOs)
write.csv(GOs, "Leaves_POP_T2_DownReg_GOs.txt", quote = F, row.names = F, col.names = F)

#Genes
genes <- SigRegR$sprot_Top_BLASTX_hit
# genes <- UpRegR$sprot_Top_BLASTX_hit

x <- sapply(strsplit(as.vector(genes), split = "_"), "[[", 1)
x <- x[which(x != ".")]
write.csv(x, "Leaves_POP_T2_DownReg_genes.txt", col.names = F, row.names = F, quote = F)

#Open output of Go Enrichment Analysis
# http://geneontology.org
# Using Arabidopsis thaliana as reference

GOEnr <- read.csv("Leaves_POP_T2_DownReg_GO_analysis_BP.txt", skip = 11, sep='\t')
head(GOEnr)
colnames(GOEnr) <- c("GO", "ref", "upload", "upload_exp", "over.under", "fold.change", "P", "FDR")
# GOEnr <- read.csv("UpRegR_analysis.txt", skip = 11, sep='\t')

#Extract GO numbers
GOEnr <- GOEnr[grep("\\(GO", as.vector(GOEnr$GO)),]
#Filter by fold.enrichment (>=10)
GOEnr <- GOEnr[which(as.numeric(GOEnr$fold.change) >= 2),]
GONum <- sapply(strsplit(as.vector(GOEnr$GO), split = "\\(GO:"), "[[", 2)
GONum <- paste("GO:", sapply(strsplit(GONum, split = ")"), "[[", 1), sep='')
#write out
head(GOEnr)
write.table(GONum, file = "Leaves_POP_T2_DownReg_GO_analysis_BP_GONum.txt", row.names = F, col.names = F, quote=F)

#Create GO graph
# http://amigo.geneontology.org/visualize?mode=client_amigo

# Create matrix of association between GO terms
# http://amigo.geneontology.org/matrix#order
