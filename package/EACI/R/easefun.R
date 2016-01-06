#' @title Calculate the DAVID EASE score for a single set
#' @usage easefun(FAll, FMk, FIdenIn)
#' @param FAll all genes that are considered in an analysis
#' @param FMk marker list of interest
#' @param FIdenIn genes that are identified as significant; it should be a subset of FAll.
#' @author Ning Leng
#' @return In FMk, only genes in FAll are included in the calculation
#' @examples f1 <- paste0("g",1:100)
#' f2 <- paste0("g",1:50)
#' f3 <- paste0("g",30:55)
#' res <- easefun(f1,f2,f3)
easefun <- function(FAll, FMk, FIdenIn){
	  FIden=intersect(FAll, FIdenIn)

	  Both=max(0,length(intersect(FMk, FIden))-1)
		MkNoIden=length(setdiff(FMk, FIden))
		IdenNoMk=length(setdiff(FIden, FMk))
		Neither=length(FAll)-Both-MkNoIden-IdenNoMk


		ToTest=cbind(c(Both, MkNoIden), c(IdenNoMk, Neither))
		Fet=fisher.test(ToTest, alternative="g")
		Out=Fet$p.value
}

