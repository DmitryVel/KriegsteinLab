#read normalized log-transformed and filtered data. Cell barcode, metadata and genes files are available at https://autism.cells.ucsc.edu/
data=readMM("normalized_matrix_Log.mtx")
cells=read.table("combined_cells_500_Genes_rm.tsv");
colnames(data)<-cells[,1]
genes=read.table("genes.tsv");
rownames(data)<-genes[,1];
genes.f=as.character(read.table("genes_F.txt")[,1]);
meta=read.table("meta_ext.txt",sep="\t",header=TRUE,row.names=1,check.names=FALSE);
data=data[genes.f,rownames(meta)];

#modify glmer function and set nAGQ to 0
trace(glmer, edit=T)
#use Shift+I to insert "nAGQ = 0" in the beginning of the function's body. Hit ESC, then :wq+ENTER to save.

#now run DEG. I recommend using 512GB of RAM and 8 cores (or scale accordingly)
Sys.setenv(MC_CORES=8);
source("MAST_ASD.R")
TEST.MAST.ASD(data,meta)
#That would run DGE analysis for all clusters with default parameters, which is testing only genes expressed in at least 10% of cells in each given cluster
#If you want to only test specific clusters, you can do this
TEST.MAST.ASD(data,meta,clusters=c("1","2"))

#Once it's done, here is how I filter the data: based on LMM_<cluster_number>.txt file, keep genes with +/-0.14 log fold change (10% absolute fold change). 
#Then, of these genes and based on <cluster_number>_FC.txt file, only keep the genes that have 1) +/-0.14 sample log fold change [last column of this file] 2) Have the same direction of fold change from LMM_<cluster_number>.txt and the sample fold change.
#Sample fold change is derived based on averaged nuclei profiles for each tissue sample and cell type.
#Basically this way you are making sure "bulk" tissue RNA levels would be in line with single-cell RNA measurements.
#I have some scripts to do filtering automatically if that would help. Can try to dig them out.
