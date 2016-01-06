#' @title Generate concordance plot of enrichment results
#' @usage PlotConcord(Out, N=100, stat="set.mean",minsize=10)
#' @param Out output of allez(), easefun() or eacifun()
#' @param N top n gene sets to show
#' @param stat sort by which summary statistics
#' @param minsize the minimal set size to consider
#' @author Ning Leng
#' @return a heatmap showing overlapping rates of top N gene sets.
#' @examples scores <- rnorm(20)


PlotConcord <- function(Out, N=100, stat="set.mean",minsize=10,...){

	Out0 <- Out[[1]][which(Out[[1]]$set.size>minsize),]
	OutO <- Out0[order(Out0[[stat]],decreasing=TRUE),]
	Sets <- OutO[1:N,]
	SetV <- Out[[2]]$set.data[[1]]
	SetIn <- which(SetV%in%names(Sets$set.size))
	SetVUse <- SetV[SetIn]
	GeneV <- Out[[2]]$set.data$symbol
	GeneVUse <- GeneV[SetIn]
	SetSize <- Sets$set.size
	SetNames <- names(SetSize)
	Scores <- Sets$set.mean
	Concord <- sapply(1:N,function(i)sapply(1:N, function(j){
									 S1=GeneVUse[which(SetVUse==SetNames[i])]
									 S2=GeneVUse[which(SetVUse==SetNames[j])]
									 length(intersect(S1,S2))/min(length(S1),length(S2))}))
	library(gplots)
	MaxSize <- max(SetSize)

	rownames(Concord)=colnames(Concord)=paste(SetNames, SetSize)
	heatmap.2(Concord,trace="none",col=greenred, Colv=F, Rowv=F,...)
	Out <- Concord	
}
