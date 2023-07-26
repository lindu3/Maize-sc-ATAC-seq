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