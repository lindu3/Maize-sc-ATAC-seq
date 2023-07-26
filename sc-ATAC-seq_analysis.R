## scripts to analysis scATAC-seq 

library(ArchR)

########################################################################
################ Genome annotation files preparation ###################
########################################################################
# Maize reference genome B73 v4
# annotate genome file
# install.packages("./BSgenome.Zmays.B73.AGPv4.tar.gz", repos = NULL, type="source")
library(BSgenome.Zmays.B73.AGPv4)
B73.v4_genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Zmays.B73.AGPv4)
B73.v4_genomeAnnotation

# annotate gene file

library(ape)
library(GenomicRanges)
gff3 = read.gff('./data/B73.AGPv4.pub.gff3', na.strings = c(".", "?"), GFF3 = TRUE)



#  I split the gff3 up into 'genes', 'mRNA' and 'exons'
genes 		= gff3[gff3$type=="gene",]
mRNA			= gff3[gff3$type=="mRNA",]
exons     = gff3[gff3$type=="exon",]

#  I parse the 9th column to get the gene and transcript IDs
mRNA$transcript_id	= gsub(";.*","",gsub("ID=","",mRNA$attributes))
mRNA$gene_id		= gsub(";.*","",gsub(".*Parent=","",mRNA$attributes))
exons$transcript_id     = gsub(";.*","",gsub(".*Parent=","",exons$attributes))
genes$gene_id           = gsub(";.*","",gsub("^ID=","",genes$attributes))

#  I create a table to convert from gene to transcripts
GeneToTranscripts	= unique(mRNA[,c("transcript_id","gene_id")])

#  I merge this table with the exons table so I can convert from exons to genes
exons		= merge(exons,GeneToTranscripts,by.x="transcript_id",by.y="transcript_id")

#  I create the GenomicRanges objects based on examples they provide
genes.gr		= GRanges(seqnames=genes$seqid,ranges = IRanges(genes$start,genes$end),strand = genes$strand,gene_id=genes$gene_id,symbol=genes$gene_id)
exons.gr		= GRanges(seqnames=exons$seqid,ranges = IRanges(exons$start,exons$end),strand = exons$strand,symbol=exons$gene_id)
tss.gr			= GRanges(seqnames=mRNA$seqid,ranges = IRanges(start=ifelse(mRNA$strand=="+",mRNA$start,mRNA$end),width=1),strand=mRNA$strand,tx_id=c(1:nrow(mRNA)),tx_name=mRNA$transcript_id) 

#  Then I run their createGeneAnnotation function
B73.v4_geneAnnotation 	 	= createGeneAnnotation(TSS=tss.gr,exons=exons.gr,genes=genes.gr)



#############################################################
################ Input data reformatting  ###################
#############################################################

## reformat the scATAC-seq data of GSE155178_ACR_x_cell.binary.sparse.txt.gz
## Before the reformat, using reformat.py to split the data based on different tissues
## then reformat the each tissue data into ArchR readable data format

# Single cell ATAC-seq data downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155178

scatac_tissue_data_path = './data/B73_tissue_Frg/'
tissue_file_list = list.files(scatac_tissue_data_path)
for (tissue in tissue_file_list){
  print(tissue)
  reformatfile = reformatFragmentFiles(
    fragmentFiles = file.path(scatac_tissue_data_path,tissue),
    checkChrPrefix = FALSE
  )
}

##################################################################################
################ all tissue single cell ATAC-seq data Analysis ###################
##################################################################################

## setting and remembering a known seed to facilitate replication of operations requiring randomization.
set.seed(1)

## setting the number of threads to be used for this study
addArchRThreads(threads = 6)


# Load the input data

all_atac_data_input = c('./data/B73_tissue_Frg/GSE155178_Ear_real_frg-Reformat.tsv.gz',
                        './data/B73_tissue_Frg/GSE155178_axillarybud1_real_frg-Reformat.tsv.gz',
                        './data/B73_tissue_Frg/GSE155178_axillarybud2_real_frg-Reformat.tsv.gz',
                        './data/B73_tissue_Frg/GSE155178_crownroot1_real_frg-Reformat.tsv.gz',
                        './data/B73_tissue_Frg/GSE155178_crownroot2_real_frg-Reformat.tsv.gz',
                        './data/B73_tissue_Frg/GSE155178_Leaf2_real_frg-Reformat.tsv.gz',
                        './data/B73_tissue_Frg/GSE155178_Root1_real_frg-Reformat.tsv.gz',
                        './data/B73_tissue_Frg/GSE155178_Root2_real_frg-Reformat.tsv.gz',
                        './data/B73_tissue_Frg/GSE155178_Tassel_real_frg-Reformat.tsv.gz',
                        './data/B73_tissue_Frg/GSE155178_B73Mo17_real_frg-Reformat.tsv.gz')
all_atac_data_input_name = c('Ear','Axillarybud1','Axillarybud2','crownroot1','crownroot2','Leaf','Root1','Root2','Tassel','B73Mo17')


## create the ArrowFiles
B73.v4.ArrowFiles.all_tissue <- createArrowFiles(
  inputFiles = all_atac_data_input,
  sampleNames = all_atac_data_input_name,
  genomeAnnotation = B73.v4_genomeAnnotation,
  geneAnnotation = B73.v4_geneAnnotation,
  minTSS = 1, #Dont set this too high because you can always increase later
  minFrags = 500, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)


## calculate the doublet scores to help with filting out the doublet of single cell studys
doubScores <- addDoubletScores(
  input = B73.v4.ArrowFiles.all_tissue,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)
## creat the ArchR project
B73.v4.all_tissue_proj <- ArchRProject(
  ArrowFiles = B73.v4.ArrowFiles.all_tissue, 
  outputDirectory = "B73_all_tissue",
  genomeAnnotation = B73.v4_genomeAnnotation,
  geneAnnotation = B73.v4_geneAnnotation,
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

## filter out the doublets
B73.v4.all_tissue_proj <- filterDoublets(ArchRProj = B73.v4.all_tissue_proj)

## dimention reduction
B73.v4.all_tissue_proj <- addIterativeLSI(ArchRProj = B73.v4.all_tissue_proj, useMatrix = "TileMatrix", name = "IterativeLSI")

## cluster by Seurat
B73.v4.all_tissue_proj <- addClusters(input = B73.v4.all_tissue_proj, reducedDims = "IterativeLSI")

## add UMAP info
B73.v4.all_tissue_proj <- addUMAP(ArchRProj = B73.v4.all_tissue_proj, reducedDims = "IterativeLSI")
## plot the cluster UMAP results by tissue and clusters seperately
p1 <- plotEmbedding(ArchRProj = B73.v4.all_tissue_proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = B73.v4.all_tissue_proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
## get the cell number count table for each clusters
table(B73.v4.all_tissue_proj$Sample)







########################################################
################ Ear tissue analysis ###################
########################################################


## create the ArrowFiles
ear_B73.v4.ArrowFiles.min <- createArrowFiles(
  inputFiles = './data/B73_tissue_Frg/GSE155178_Ear_real_frg-Reformat.tsv.gz',
  sampleNames = "Ear_500",
  genomeAnnotation = B73.v4_genomeAnnotation,
  geneAnnotation = B73.v4_geneAnnotation,
  minTSS = 1, #Dont set this too high because you can always increase later
  filterFrags = 500, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)


## calculate the doublet scores to help with filting out the doublet of single cell studys
doubScores <- addDoubletScores(
  input = ear_B73.v4.ArrowFiles.min,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

## creat the ArchR project
proj <- ArchRProject(
  ArrowFiles = ear_B73.v4.ArrowFiles.min, 
  outputDirectory = "B73",
  genomeAnnotation = B73.v4_genomeAnnotation,
  geneAnnotation = B73.v4_geneAnnotation,
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

## check the availabel matriuces of this proj
getAvailableMatrices(proj)

## check the distribution of TSS enrichemnt score 
quantile(proj$TSSEnrichment)

## filter out the doublets
proj <- filterDoublets(ArchRProj = proj)

## dimention reduction
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")


## TSS Enrichment plot
TSSEnrichment_plot <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE,
  baseSize = 20
)
TSSEnrichment_plot

## Number of fragments plot

num_of_frag <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE,
  baseSize = 20
)
num_of_frag
ggAlignPlots( TSSEnrichment_plot,num_of_frag, type = "h")

#############################
## cluster by Seurat
Seraut_cluster_proj <- addClusters(input = proj, reducedDims = "IterativeLSI", maxClusters = 15 )

Seraut_cluster_proj <- addUMAP(ArchRProj = Seraut_cluster_proj, reducedDims = "IterativeLSI")

#p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
# plotting
Seurat_cluster_p <- plotEmbedding(ArchRProj = Seraut_cluster_proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")



## plotting
ggAlignPlots( Seurat_cluster_p, type = "h")


table(Seraut_cluster_proj$Clusters)


############################
## cluster by Seurat
## set the min cell number for the clusters (nOurlier) 30
Seraut_cluster_proj_nOutlier30 <- addClusters(input = proj, reducedDims = "IterativeLSI", maxClusters = 7, nOutlier = 30 )

Seraut_cluster_proj_nOutlier30 <- addUMAP(ArchRProj = Seraut_cluster_proj_nOutlier30, reducedDims = "IterativeLSI")

#p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

Seurat_cluster_nOutlier30_p <- plotEmbedding(ArchRProj = Seraut_cluster_proj_nOutlier30, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

## ploting
ggAlignPlots(Seurat_cluster_nOutlier30_p, type = "h")
## get the cell number count table for each clusters
table(Seraut_cluster_proj_nOutlier30$Clusters)

## find ear tissue marker genes
ear_markersGS_30 <- getMarkerFeatures(
  ArchRProj = Seraut_cluster_proj_nOutlier30, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

## get the ear tissue marker gene list
ear_markerList <- getMarkers(ear_markersGS_30, cutOff = "FDR <= 0.05 & Log2FC >= 1")
ear_markerList$C2
## write the marker genes
write.csv(ear_markerList$C7,file="./output/ear_marker_genes/C7.csv")
## annotation by cluster specific genes

markerGenes  <- c(
  "Zm00001d021906", #krn related
  "Zm00001d031254",  # C1
  "Zm00001d021791",  # C2
  "Zm00001d021796",  # C6
  "Zm00001d052890", #UB3
  "Zm00001d031463"
  
  
)
## plot the heatmap based on the marker genes
heatmapGS <- markerHeatmap(
  seMarker = ear_markersGS_30, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")



## Because of the sparsity of scATAC-seq data. 
## We can impute gene scores by smoothing signal across nearby cells.
Seraut_cluster_proj_nOutlier30 <- addImputeWeights(Seraut_cluster_proj_nOutlier30)

p_marker_gene <- plotEmbedding(
  ArchRProj = Seraut_cluster_proj_nOutlier30, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(Seraut_cluster_proj_nOutlier30)
)
## plot the specific gene enrichment map
p_marker_gene$Zm00001d021906
p_marker_gene$Zm00001d031463
p_marker_gene$Zm00001d052890
## Plotting 

p_all <- lapply(p_marker_gene, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 2),p_all))


cluster2_proj <- addGroupCoverages(ArchRProj = Seraut_cluster_proj_nOutlier30, groupBy = "Clusters")


## call ATAC-seq peaks by MACS2
pathToMacs2 <- findMacs2()

macs2_test <- addReproduciblePeakSet(
  ArchRProj = cluster2_proj, 
  groupBy = "Clusters", 
  pathToMacs2 = pathToMacs2,
  genomeSize = 2134373047 # B73.v4
)
## get single cell ATAC-seq peak data
macs2_test
getPeakSet(macs2_test)

ear_proj_peakcalling <- addPeakMatrix(macs2_test)

getAvailableMatrices(ear_proj_peakcalling)
ear_proj_peakcalling
############################################
## get the marker peaks
markersPeaks <- getMarkerFeatures(
  ArchRProj = ear_proj_peakcalling, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaks

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList

markerList$C2


# plot the MA plot of different clusters
pmaC2 <- markerPlot(seMarker = markersPeaks, name = "C2", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
pmaC2



####################################
## plot the Marker features C2 VS C6
markerTest <- getMarkerFeatures(
  ArchRProj = ear_proj_peakcalling, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C2",
  bgdGroups = "C6"
)
# plot the MA plot
pma_C2VC6 <- markerPlot(seMarker = markerTest, name = "C2", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma_C2VC6
# plot the Volcanio plot 
pv <- markerPlot(seMarker = markerTest, name = "C2", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
pv
##################################
## Motif analysis

## downloading chromVARmotif package to do the motif analysis
# devtools::install_github("GreenleafLab/chromVARmotifs")

ear_proj_peakcalling_motif <- addMotifAnnotations(ArchRProj = ear_proj_peakcalling, motifSet = "homer", name = "Motif")


## for Cluster 2
motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = ear_proj_peakcalling_motif,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
motifsUp
## data for plotting
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

##
head(df)
## ggplot for motif ranking
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp

## for Cluster 6
motifsDo <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = ear_proj_peakcalling_motif,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)
## motif p-value df for plotting 
df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

## ggplot for C6 plotting
ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo


####################################
##   Integrate sc RNA-seq data   ###
####################################

# ArchR  accepts both RangedSummarizedExperiment object and unmodified Seurat objects as input 
# to the integration workflow
# ear_single cell RNA-seq data is from doi:10.1016/j.devcel.2020.12.015 in Seurat object format.

ear_RNA_rdata <- load("./scRNA_data/ENAD51_seurat_obj_LinDu.RData")
ear_RNA_rdata = embl_ebi.srt

integrate_RNA_ATAC <- addGeneIntegrationMatrix(
  ArchRProj = Seraut_cluster_proj_nOutlier30, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = ear_RNA_rdata,
  addToArrow = FALSE,
  groupRNA = "inferred_cell_type_authors_labels",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)
# count matrix
cM <- as.matrix(confusionMatrix(integrate_RNA_ATAC$Clusters, integrate_RNA_ATAC$predictedGroup_Un))
cM
# since there are a lot of cells are labeled as NA we get rid of them and lable the rest of the cells
cM <- cM[,-1]
cM
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments


pal <- paletteDiscrete(values = ear_RNA_rdata$inferred_cell_type_authors_labels)
pal
Seraut_cluster_proj_nOutlier30$cellNames
p1 <- plotEmbedding(
  integrate_RNA_ATAC, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un", 
  pal = pal
)
p1

##############################

# plot the ATAC-seq peak tracks based on interested gene Symbol
p <- plotBrowserTrack(
  ArchRProj = ear_proj_peakcalling_motif, 
  groupBy = "Clusters", 
  geneSymbol = 'Zm00001d021906', 
  upstream = 10000,
  downstream = 10000
)

grid::grid.newpage()
grid::grid.draw(p$Zm00001d021906)



# chr7 165,587,952 171,432,419

p_all7 <- plotBrowserTrack(
  ArchRProj = ear_proj_peakcalling_motif, 
  groupBy = "Clusters", 
  geneSymbol = 'Zm00001d021907', 
  upstream = 10000,
  downstream = 8000000
)

grid::grid.newpage()
grid::grid.draw(p_all7$Zm00001d021907)

# chr1 183,922,026 209,130,885

p_all1 <- plotBrowserTrack(
  ArchRProj = ear_proj_peakcalling_motif, 
  groupBy = "Clusters", 
  geneSymbol = 'Zm00001d031254', 
  upstream = 10000,
  downstream = 20000000
)

grid::grid.newpage()
grid::grid.draw(p_all1$Zm00001d031254)

