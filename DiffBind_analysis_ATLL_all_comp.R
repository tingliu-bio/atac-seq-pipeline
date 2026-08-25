#load the package
library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# read-in sampleSheet.csv
#samplesheet_path="/cluster/projects/kridelgroup/ATAC_seq_SUDHL2_HBL1_YX968_24h/pipeline-chromatin-accessibility/data/R_scripts/sampleSheet"
#samplesheet_path="/cluster/projects/kridelgroup/Noorhan/ATAC-seq/pipeline-chromatin-accessibility/R_scripts/sampleSheet"
#samplesheet_path="/cluster/projects/kridelgroup/Noorhan/ATAC-seq/pipeline-chromatin-accessibility/data/R_scripts_ATLL1K/sampleSheet"
samplesheet_path="/cluster/projects/kridelgroup/Noorhan/ATAC-seq/pipeline-chromatin-accessibility/data/R_scripts_analysis_HTLV1_run/ATLL1K/sampleSheet_ATL1K"
##Val2contro comp
#load sample sheet
dbObj <- dba(sampleSheet= paste(samplesheet_path, "sampleSheet_Val2control_ATLL.csv", sep="/"))

#setwd('/cluster/projects/kridelgroup/Noorhan/ATAC-seq/pipeline-chromatin-accessibility/R_scripts/ATLL_output')
#setwd('/cluster/projects/kridelgroup/Noorhan/ATAC-seq/pipeline-chromatin-accessibility/data/R_scripts_ATLL1K/ATLL1K_output')
setwd('/cluster/projects/kridelgroup/Noorhan/ATAC-seq/pipeline-chromatin-accessibility/data/R_scripts_analysis_HTLV1_run/ATLL1K')
##DMSO2PROTAC_HBL1
# gen the Heatmap
pdf("Val2control_ATLL.pdf")
par(oma = c(4, 4, 2, 2) +1.5)
dba.plotHeatmap(dbObj)
dev.off()

# cal the diff
## Counts reads for each sample
dbObj <-dba.blacklist(dbObj)
dbObj <- dba.count(dbObj, bParallel=FALSE)  # read in the counts, turn off the Parallel computing
dbObj <-dba.normalize(dbObj)
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION ,minMembers=2)

##terminal mode, run these two command to only use one core 
#library(BiocParallel)
#register(MulticoreParam(1)) 

#check the samples
cat ("OBJ samples","\n")
dbObj$samples
#DBA_CONDITION =c("DMSO", "PROTAC")

#db/OBJ samplesObj <- dba.contrast(dbObj, contrast=c("Condition","Control","Treatment"), categories=DBA_CONDITION, minMembers=2)  # set up the contrast
#dbObj <- dba.contrast(dbObj, design="~Condition",minMembers=2) 
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION ,minMembers=2) 
cat ("print dbObj")
cat ("\n")
dbObj

cat ("check if contrasts is true","\n")
dba.show(dbObj, bContrasts=TRUE)

cat ("conduct dba.analyze","\n")
#dbObj <- dba.analyze(dbObj, design="~Condition")  # analyze
#force to run the dba.analyze() to get away with the error message below
#Analyze error: Error in pv.DBA(DBA, method, bTagwise = bTagwise, minMembers = 3, bParallel = bParallel)
dbObj <- dba.analyze(dbObj,method=DBA_DESEQ2)

cat ("check contrasts again","\n")
dba.show(dbObj, bContrasts=TRUE)

# diff peaks
diff_peaks <- dba.report(dbObj, th=0.01)  # choose FDR < 0.05
head(diff_peaks)
write.table(diff_peaks, file="diff_peaks_FDR_sampleSheet_Val2control_ATLL_th0.01.txt", quote=F, row.names = F, sep="\t")

#annotated_peaks <- annotatePeak(diff_peaks, tssRegion = c(-1000, 1000),TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
annotated_peaks <- annotatePeak(diff_peaks, tssRegion = c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")

write.table(annotated_peaks, file="annotated_diff_peaks_FDR_Val2control_ATLL_th0.01.txt", quote=F, row.names = F, sep="\t")

pdf("DiffBind_diff_peaks_annotated_peaks_Val2control_ATLL.pdf")
plotAnnoPie(annotated_peaks)
dev.off()


##
cat ("show diff peaks","\n")
dba.show(diff_peaks, bContrasts=TRUE)


diff_peaks <- dba.contrast(dbObj, categories=DBA_CONDITION, minMembers=2)
diff_peaks <- dba.analyze(diff_peaks, method=DBA_DESEQ2)

# check the contrast again
cat ("show diff peaks again","\n")
dba.show(diff_peaks, bContrasts=TRUE)

cat ("print(DBA_CONDITION)","\n")
print(DBA_CONDITION)
cat ("print(DBA_ID)","\n")
print(DBA_ID)

# gen PCA plot
cat ("gen PCA plot","\n")
pdf("DiffBind_diff_peaks_PCA_Val2control_ATLL.pdf")
dba.plotPCA(diff_peaks, label=DBA_ID,correlations=FALSE)
dev.off()

# gen Volcano plot
cat ("generate the Volcano plot","\n")
pdf("DiffBind_Volcano_diff_peaks_Val2control_ATLL.pdf")
dba.plotVolcano(diff_peaks,contrast=1)
dev.off()


##combi2control comp
#load sample sheet
dbObj <- dba(sampleSheet= paste(samplesheet_path, "sampleSheet_combi2control_ATLL.csv", sep="/"))

#########################################################3
#combi2control_ATLL
# gen the Heatmap
pdf("DiffBind_Heatmap_combi2control_ATLL.pdf")
par(oma = c(4, 4, 2, 2) +1.5)
dba.plotHeatmap(dbObj)
dev.off()

# cal the diff
## Counts reads for each sample
dbObj <-dba.blacklist(dbObj)
dbObj <- dba.count(dbObj, bParallel=FALSE)  # read in the counts, turn off the Parallel computing
dbObj <-dba.normalize(dbObj)
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION ,minMembers=2)

##terminal mode, run these two command to only use one core 
#library(BiocParallel)
#register(MulticoreParam(1)) 

#check the samples
cat ("OBJ samples","\n")
dbObj$samples


#db/OBJ samplesObj <- dba.contrast(dbObj, contrast=c("Condition","Control","Treatment"), categories=DBA_CONDITION, minMembers=2)  # set up the contrast
#dbObj <- dba.contrast(dbObj, design="~Condition",minMembers=2) 
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION ,minMembers=2) 
cat ("print dbObj")
cat ("\n")
dbObj

cat ("check if contrasts is true","\n")
dba.show(dbObj, bContrasts=TRUE)

cat ("conduct dba.analyze","\n")
#dbObj <- dba.analyze(dbObj, design="~Condition")  # analyze
#force to run the dba.analyze() to get away with the error message below
#Analyze error: Error in pv.DBA(DBA, method, bTagwise = bTagwise, minMembers = 3, bParallel = bParallel)
dbObj <- dba.analyze(dbObj,method=DBA_DESEQ2)

cat ("check contrasts again","\n")
dba.show(dbObj, bContrasts=TRUE)

# diff peaks
diff_peaks <- dba.report(dbObj, th=0.01)  # choose FDR < 0.05
head(diff_peaks)
write.table(diff_peaks, file="diff_peaks_FDR_combi2control_ATLL_th0.01.txt", quote=F, row.names = F, sep="\t")
#annotated_peaks <- annotatePeak(diff_peaks, tssRegion = c(-1000, 1000),TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene) 

annotated_peaks <- annotatePeak(diff_peaks, tssRegion = c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
write.table(annotated_peaks, file="annotated_diff_peaks_FDR_combi2control_ATLL_th0.01.txt", quote=F, row.names = F, sep="\t")

pdf("DiffBind_diff_peaks_annotated_peaks_combi2control_ATLL.pdf")
plotAnnoPie(annotated_peaks)
dev.off()
##
cat ("show diff peaks","\n")
dba.show(diff_peaks, bContrasts=TRUE)


diff_peaks <- dba.contrast(dbObj, categories=DBA_CONDITION, minMembers=2)
diff_peaks <- dba.analyze(diff_peaks, method=DBA_DESEQ2)

# check the contrast again
cat ("show diff peaks again","\n")
dba.show(diff_peaks, bContrasts=TRUE)

cat ("print(DBA_CONDITION)","\n")
print(DBA_CONDITION)
cat ("print(DBA_ID)","\n")
print(DBA_ID)

# gen PCA plot
cat ("gen PCA plot","\n")
pdf("DiffBind_diff_peaks_PCA_combi2control_ATLL.pdf")
dba.plotPCA(diff_peaks, label=DBA_ID,correlations=FALSE)
dev.off()

# gen Volcano plot
cat ("generate the Volcano plot","\n")
pdf("DiffBind_Volcano_diff_peaks_combi2control_ATLL.pdf")
dba.plotVolcano(diff_peaks,contrast=1)
dev.off()




##Val2combi comp
#load sample sheet
dbObj <- dba(sampleSheet= paste(samplesheet_path, "sampleSheet_Val2combi_ATLL.csv", sep="/"))

#########################################################3
#Val2combi_ATLL
# gen the Heatmap
pdf("DiffBind_Heatmap_Val2combi_ATLL.pdf")
par(oma = c(4, 4, 2, 2) +1.5)
dba.plotHeatmap(dbObj)
dev.off()

# cal the diff
## Counts reads for each sample
dbObj <-dba.blacklist(dbObj)
dbObj <- dba.count(dbObj, bParallel=FALSE)  # read in the counts, turn off the Parallel computing
dbObj <-dba.normalize(dbObj)
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION ,minMembers=2)

##terminal mode, run these two command to only use one core 
#library(BiocParallel)
#register(MulticoreParam(1)) 

#check the samples
cat ("OBJ samples","\n")
dbObj$samples


#db/OBJ samplesObj <- dba.contrast(dbObj, contrast=c("Condition","Control","Treatment"), categories=DBA_CONDITION, minMembers=2)  # set up the contrast
#dbObj <- dba.contrast(dbObj, design="~Condition",minMembers=2) 
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION ,minMembers=2) 
cat ("print dbObj")
cat ("\n")
dbObj

cat ("check if contrasts is true","\n")
dba.show(dbObj, bContrasts=TRUE)

cat ("conduct dba.analyze","\n")
#dbObj <- dba.analyze(dbObj, design="~Condition")  # analyze
#force to run the dba.analyze() to get away with the error message below
#Analyze error: Error in pv.DBA(DBA, method, bTagwise = bTagwise, minMembers = 3, bParallel = bParallel)
dbObj <- dba.analyze(dbObj,method=DBA_DESEQ2)

cat ("check contrasts again","\n")
dba.show(dbObj, bContrasts=TRUE)

# diff peaks
diff_peaks <- dba.report(dbObj, th=0.01)  # choose FDR < 0.05
head(diff_peaks)
write.table(diff_peaks, file="diff_peaks_FDR_Val2combi_ATLL_th0.01.txt", quote=F, row.names = F, sep="\t")
#annotated_peaks <- annotatePeak(diff_peaks, tssRegion = c(-1000, 1000),TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene) 

annotated_peaks <- annotatePeak(diff_peaks, tssRegion = c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
write.table(annotated_peaks, file="annotated_diff_peaks_FDR_Val2combi_ATLL_th0.01.txt", quote=F, row.names = F, sep="\t")

pdf("DiffBind_diff_peaks_annotated_peaks_Val2combi_ATLL.pdf")
plotAnnoPie(annotated_peaks)
dev.off()
##
cat ("show diff peaks","\n")
dba.show(diff_peaks, bContrasts=TRUE)


diff_peaks <- dba.contrast(dbObj, categories=DBA_CONDITION, minMembers=2)
diff_peaks <- dba.analyze(diff_peaks, method=DBA_DESEQ2)

# check the contrast again
cat ("show diff peaks again","\n")
dba.show(diff_peaks, bContrasts=TRUE)

cat ("print(DBA_CONDITION)","\n")
print(DBA_CONDITION)
cat ("print(DBA_ID)","\n")
print(DBA_ID)

# gen PCA plot
cat ("gen PCA plot","\n")
pdf("DiffBind_diff_peaks_PCA_Val2combi_ATLL.pdf")
dba.plotPCA(diff_peaks, label=DBA_ID,correlations=FALSE)
dev.off()

# gen Volcano plot
cat ("generate the Volcano plot","\n")
pdf("DiffBind_Volcano_diff_peaks_Val2combi_ATLL.pdf")
dba.plotVolcano(diff_peaks,contrast=1)
dev.off()


##Ent2control comp
#load sample sheet
dbObj <- dba(sampleSheet= paste(samplesheet_path, "sampleSheet_Ent2control_ATLL.csv", sep="/"))

#########################################################3
#Ent2control_ATLL
# gen the Heatmap
pdf("DiffBind_Heatmap_Ent2control_ATLL.pdf")
par(oma = c(4, 4, 2, 2) +1.5)
dba.plotHeatmap(dbObj)
dev.off()

# cal the diff
## Counts reads for each sample
dbObj <-dba.blacklist(dbObj)
dbObj <- dba.count(dbObj, bParallel=FALSE)  # read in the counts, turn off the Parallel computing
dbObj <-dba.normalize(dbObj)
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION ,minMembers=2)

##terminal mode, run these two command to only use one core 
#library(BiocParallel)
#register(MulticoreParam(1)) 

#check the samples
cat ("OBJ samples","\n")
dbObj$samples


#db/OBJ samplesObj <- dba.contrast(dbObj, contrast=c("Condition","Control","Treatment"), categories=DBA_CONDITION, minMembers=2)  # set up the contrast
#dbObj <- dba.contrast(dbObj, design="~Condition",minMembers=2) 
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION ,minMembers=2) 
cat ("print dbObj")
cat ("\n")
dbObj

cat ("check if contrasts is true","\n")
dba.show(dbObj, bContrasts=TRUE)

cat ("conduct dba.analyze","\n")
#dbObj <- dba.analyze(dbObj, design="~Condition")  # analyze
#force to run the dba.analyze() to get away with the error message below
#Analyze error: Error in pv.DBA(DBA, method, bTagwise = bTagwise, minMembers = 3, bParallel = bParallel)
dbObj <- dba.analyze(dbObj,method=DBA_DESEQ2)

cat ("check contrasts again","\n")
dba.show(dbObj, bContrasts=TRUE)

# diff peaks
diff_peaks <- dba.report(dbObj, th=0.01)  # choose FDR < 0.05
head(diff_peaks)
write.table(diff_peaks, file="diff_peaks_FDR_Ent2control_ATLL_th0.01.txt", quote=F, row.names = F, sep="\t")
#annotated_peaks <- annotatePeak(diff_peaks, tssRegion = c(-1000, 1000),TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene) 

annotated_peaks <- annotatePeak(diff_peaks, tssRegion = c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
write.table(annotated_peaks, file="annotated_diff_peaks_FDR_Ent2control_ATLL_th0.01.txt", quote=F, row.names = F, sep="\t")

pdf("DiffBind_diff_peaks_annotated_peaks_Ent2control_ATLL.pdf")
plotAnnoPie(annotated_peaks)
dev.off()
##
cat ("show diff peaks","\n")
dba.show(diff_peaks, bContrasts=TRUE)


diff_peaks <- dba.contrast(dbObj, categories=DBA_CONDITION, minMembers=2)
diff_peaks <- dba.analyze(diff_peaks, method=DBA_DESEQ2)

# check the contrast again
cat ("show diff peaks again","\n")
dba.show(diff_peaks, bContrasts=TRUE)

cat ("print(DBA_CONDITION)","\n")
print(DBA_CONDITION)
cat ("print(DBA_ID)","\n")
print(DBA_ID)

# gen PCA plot
cat ("gen PCA plot","\n")
pdf("DiffBind_diff_peaks_PCA_Ent2control_ATLL.pdf")
dba.plotPCA(diff_peaks, label=DBA_ID,correlations=FALSE)
dev.off()

# gen Volcano plot
cat ("generate the Volcano plot","\n")
pdf("DiffBind_Volcano_diff_peaks_Ent2control_ATLL.pdf")
dba.plotVolcano(diff_peaks,contrast=1)
dev.off()


##Ent2combi comp
#load sample sheet
dbObj <- dba(sampleSheet= paste(samplesheet_path, "sampleSheet_Ent2combi_ATLL.csv", sep="/"))

#########################################################3
#Ent2combi_ATLL_ATLL
# gen the Heatmap
pdf("DiffBind_Heatmap_Ent2combi_ATLL_ATLL.pdf")
par(oma = c(4, 4, 2, 2) +1.5)
dba.plotHeatmap(dbObj)
dev.off()

# cal the diff
## Counts reads for each sample
dbObj <-dba.blacklist(dbObj)
dbObj <- dba.count(dbObj, bParallel=FALSE)  # read in the counts, turn off the Parallel computing
dbObj <-dba.normalize(dbObj)
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION ,minMembers=2)

##terminal mode, run these two command to only use one core 
#library(BiocParallel)
#register(MulticoreParam(1)) 

#check the samples
cat ("OBJ samples","\n")
dbObj$samples


#db/OBJ samplesObj <- dba.contrast(dbObj, contrast=c("Condition","Control","Treatment"), categories=DBA_CONDITION, minMembers=2)  # set up the contrast
#dbObj <- dba.contrast(dbObj, design="~Condition",minMembers=2) 
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION ,minMembers=2) 
cat ("print dbObj")
cat ("\n")
dbObj

cat ("check if contrasts is true","\n")
dba.show(dbObj, bContrasts=TRUE)

cat ("conduct dba.analyze","\n")
#dbObj <- dba.analyze(dbObj, design="~Condition")  # analyze
#force to run the dba.analyze() to get away with the error message below
#Analyze error: Error in pv.DBA(DBA, method, bTagwise = bTagwise, minMembers = 3, bParallel = bParallel)
dbObj <- dba.analyze(dbObj,method=DBA_DESEQ2)

cat ("check contrasts again","\n")
dba.show(dbObj, bContrasts=TRUE)

# diff peaks
diff_peaks <- dba.report(dbObj, th=0.01)  # choose FDR < 0.05
head(diff_peaks)
write.table(diff_peaks, file="diff_peaks_FDR_Ent2combi_ATLL_ATLL_th0.01.txt", quote=F, row.names = F, sep="\t")
#annotated_peaks <- annotatePeak(diff_peaks, tssRegion = c(-1000, 1000),TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene) 

annotated_peaks <- annotatePeak(diff_peaks, tssRegion = c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
write.table(annotated_peaks, file="annotated_diff_peaks_FDR_Ent2combi_ATLL_ATLL_th0.01.txt", quote=F, row.names = F, sep="\t")

pdf("DiffBind_diff_peaks_annotated_peaks_Ent2combi_ATLL_ATLL.pdf")
plotAnnoPie(annotated_peaks)
dev.off()
##
cat ("show diff peaks","\n")
dba.show(diff_peaks, bContrasts=TRUE)


diff_peaks <- dba.contrast(dbObj, categories=DBA_CONDITION, minMembers=2)
diff_peaks <- dba.analyze(diff_peaks, method=DBA_DESEQ2)

# check the contrast again
cat ("show diff peaks again","\n")
dba.show(diff_peaks, bContrasts=TRUE)

cat ("print(DBA_CONDITION)","\n")
print(DBA_CONDITION)
cat ("print(DBA_ID)","\n")
print(DBA_ID)

# gen PCA plot
cat ("gen PCA plot","\n")
pdf("DiffBind_diff_peaks_PCA_Ent2combi_ATLL_ATLL.pdf")
dba.plotPCA(diff_peaks, label=DBA_ID,correlations=FALSE)
dev.off()

# gen Volcano plot
cat ("generate the Volcano plot","\n")
pdf("DiffBind_Volcano_diff_peaks_Ent2combi_ATLL_ATLL.pdf")
dba.plotVolcano(diff_peaks,contrast=1)
dev.off()




