#QBiC Tuebingen - Damien Dos Santos Luis
#08/09/11

#Script for ChIPseq analysis and data prepartion for indexing in Qdataminer

library("ChIPseeker", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.2")
library("BSgenome.Mmusculus.UCSC.mm10", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.2")
require("org.Mm.eg.db")
require('GenomicFeatures')
library("TxDb.Mmusculus.UCSC.mm10.knownGene", lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.2")

setwd('/home/damiendsl/data_analysis_toolbox/datasets/GSE71250/macs2_srr24/')
files=list.files()
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene  #Reference for annotation

peaks= readPeakFile(files[[4]], as = "GRanges", header=FALSE)
colnames(mcols(peaks))=c("id", "score", "pileup_heigth", "fold_change", "-log10(pvalue)", "-log10(qvalue)", "summit_position_to_peak_start") 


#Profile of ChIP binding to TSS regions

promoter = getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrix = getTagMatrix(peaks, windows = promoter)
peakHeatmap(files[[4]], TxDb = txdb, upstream = 3000, downstream = 3000, color = "red")

#Average profile of ChIP peaks binding to TSS region

plotAvgProf(tagMatrix, xlim = c(-3000, 3000), xlab = "Genomic Region (5'->3')", ylab = "Read Count Frequency")
#plotAvgProf(tagMatrix, xlim = c(-3000, 3000), conf = 0.95, resample = 500)

#Peaks annotation

peakAnno = annotatePeak(peaks, tssRegion = c(-3000, 3000), TxDb = txdb,  
                        genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"), 
                        annoDb="org.Mm.eg.db", verbose = TRUE)

write.table(peakAnno, file = "peakAnno.txt")

#Choose the right genome
require('BSgenome')
require("BSgenome.Mmusculus.UCSC.mm10")
df_peakAnno=read.table('peakAnno.txt')

#load table containing Macs statistics
#macs_stats=read.table(files[[5]], skip=23, header=TRUE)

#creation of metadata fields
df_peakAnno$sample <- rep("sample_name",nrow(df_peakAnno)) 
df_peakAnno$type_experiment <- rep("type_exp",nrow(df_peakAnno))
df_peakAnno$TF_study <- rep("TF",nrow(df_peakAnno)) 
df_peakAnno$replicate <- rep(num_rep,nrow(df_peakAnno)) 

#Finding DNA sequences corresponding to window
genome = BSgenome.Mmusculus.UCSC.mm10
seqPeaks=data.frame(chr=df_peakAnno$geneChr, start=df_peakAnno$start, end=df_peakAnno$end, strand=df_peakAnno$strand)
Peak_DNA_sequence=as.data.frame(getSeq(genome, seqPeaks$chr, start=seqPeaks$start, end=seqPeaks$end, strand=seqPeaks$strand))
write.csv(Peak_DNA_sequence, 'Peak_DNA_sequence.csv')

#and eventually the entire gene sequence
#seqGenes=data.frame(chr=df_peakAnno$geneChr, start=df_peakAnno$geneStart, end=df_peakAnno$geneEnd)
#Gene_sequence=getSeq(genome, seqGenes$chr, start=seqGenes$start, end=seqGenes$end)

#Creation of the report
pdf(file = 'peekSeeker_report.pdf', width=8.50, height=11)
par(mfrow=(c(1,4)))

#ChIP peaks coverage plot
#Covplot requires the 10log10_pvalue column > The name of the column can change depending on the species
#V5 has been controlled for mice and humans - V4 for Drosophila
covplot(peaks, weightCol = "score")
#Average profile of ChIP peaks binding to TSS region
plotAvgProf(tagMatrix, xlim = c(-3000, 3000), xlab = "Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAnnoPie(peakAnno)
plotDistToTSS(peakAnno, title = "Distribution of transcrition factor-binding loci \n relative to TSS")
dev.off()
