#Installation
#biocLite("pasilla")
#biocLite('DESeq')

library('DESeq')
library('pasilla')
##load data set
datafile = system.file( "extdata/pasilla_gene_counts.tsv", package="pasilla" )
pasillaCountTable = read.table( datafile, header=TRUE, row.names=1 )

##remove genes with lower reads number
cutoffReads       <- quantile(as.matrix(pasillaCountTable), 0.5)
pasillaCountTableTrim <- pasillaCountTable[which(rowSums(pasillaCountTable) > cutoffReads ), ]


##a description of the samples

condition <-    c( "untreated", "untreated", "untreated",
                 "untreated", "treated", "treated", "treated" )


cds      <- newCountDataSet( pasillaCountTableTrim, condition )
cds      <- estimateSizeFactors(cds)
head( counts( cds, normalized=TRUE ) )
cds      <- estimateDispersions( cds )

str( fitInfo(cds) )
plotDispEsts(cds)
res        <- nbinomTest( cds, 'untreated', 'treated' )
plotMA(res)


##check P-value distribution
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="") ## check the p-value
resSig = res[ which(res$padj < 0.01), ] ## extract the genes with adjP < 0.1


##extract normalized data
pasillaNormalizedDf <- counts( cds, normalized=TRUE )
heatMapDf           <- pasillaNormalizedDf[(rownames(pasillaNormalizedDf) %in% resSig$id), ]
heatMapDfNor        <- t(apply(heatMapDf, 1, function(x){x/max(x)}))


##plot heatmap

library("RColorBrewer")
library("gplots")
## top 30 genes with highest row means 
#select = order(rowMeans(counts(cdsFull)), decreasing=TRUE)[1:30] 
hmcol =  colorRampPalette(brewer.pal(9, "GnBu"))(100)
## expression heatmap
heatmap.2(heatMapDfNor, col = hmcol, trace="none", margin=c(10, 6))





