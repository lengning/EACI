#' @title Run EACI test for gene set enrichment analysis
#' @usage eacitest(scores, lib, idtype=c("ENTREZID", "SYMBOL"),library.loc=NULL,sets = c("GO","KEGG"),
#' locallist = NULL,reduce = NULL,iter=10,minsetsize=10,annotate = TRUE)
#' @param scores numeric vector for gene scores. names should be gene names (IDs)
#' @param lib character string, name of data package corresponding to
#' organism; for example, "org.Hs.eg".
#' @param idtype idtype could be either \code{"ENTREZID"} (default) or 
#'	\code{"SYMBOL"}. It should match names of the scores vector.
#' @param library.loc library location
#' @param sets character string, describing the collection of sets. \code{"GO"} (default) or \code{"KEGG"}
#' @param locallist list contains in-house gene sets of interest. Default is NULL.
#' Each element in the list represents a gene set. Each element should be a character
#' string containing genes' entrez IDs (if \code{idtype="ENTREZID"i}) or
#' gene symbols (if \code{idtype="SYMBOL"}). Element names will be used as 
#' set names. If \code{locallist} is not NULL, \code{universe="local"} will be disabled.
#' @param reduce function used in \code{collapse}. 
#' @param iter number of iterations when running EM
#' @param minsetsize minimal set size to consider; sets with size less than minsetsize will be
#' ignored during the calculation.
#' @param annotate logical, whether to include set names in the output
#' @author Ning Leng
#' @return sorted list of gene sets; summary statistics
#' @examples scores <- rnorm(20)

eacitest <- function (scores,
                   lib,
                   idtype = c("ENTREZID", "SYMBOL"),
                   library.loc=NULL,
                   sets = c("GO","KEGG"),
                   locallist = NULL,
		reduce = NULL,
		iter=10,
		minsetsize=10,
                   annotate = TRUE,
                   ...)
{
  stopifnot(any(!is.na(scores)))

  ## Remove NA's from scores ##
  if (any(is.na(scores))){
    warning("scores containing NA's will be excluded")
    scores <- scores[!is.na(scores)]
  }

  scorenames <- names(scores)
  sets <- match.arg(sets)
	idtype <- match.arg(idtype)

  if (any(duplicated(scorenames))) 
    stop("input IDs must be unique")
  if( !is.numeric(scores) ) 
    stop("scores must be numeric")
    ## Default reduce ##
  uscores <- unique(scores)
  if(is.null(reduce))
   reduce <- if(length(uscores)==2 & identical(uscores,c(0,1))) max else median

  message("Loading necessary libraries...")
  fn_loadPlatformLibraries( Libraries=lib, library.loc=library.loc )

  set_id <- switch(sets, GO="go_id", KEGG="path_id")
  is.org <- substr(lib,1,3)=="org"
  orgpkg <- ifelse(is.org,lib[1],paste(lib[1],"ORGPKG",sep=""))
	idtype_name <- switch(idtype, ENTREZID="gene_id", SYMBOL="symbol")

  ## ANNOTATION ##
  message("Converting annotations to data.frames ...")

    ## Use org info for ENTREZ TO GO/KEGG ID ##
    set2eg <- toTable(getDataEnv(name=ifelse(sets=="GO",
              "GO2ALLEGS","PATH2EG"), lib=orgpkg))
    org.symbol <- toTable(getDataEnv(name="SYMBOL",lib=orgpkg))
    set2eg <- cbind(set2eg, org.symbol[
              match(set2eg$gene_id,org.symbol$gene_id),"symbol",drop=FALSE])
# local list
		if(!is.null(locallist)){
		if(is.null(names(locallist)))names(locallist) = paste0("L",1:length(locallist))
		localunlist=unlist(locallist)
		locallen=length(locallist)
		localeachlen=sapply(locallist,length)
		newlist <- switch(idtype, ENTREZID=data.frame(cbind(gene_id=localunlist,
			go_id=unlist(sapply(1:locallen,function(k)rep(paste0("Local:",names(locallist)[k]),localeachlen[k]),simplify=F)),
			Evidence="local", Ontology="local",symbol=set2eg$symbol[match(localunlist,set2eg$gene_id)]
		),stringsAsFactors=F),
		SYMBOL=data.frame(cbind(gene_id=set2eg$gene_id[match(localunlist,set2eg$symbol)], 
			go_id=unlist(sapply(1:locallen,function(k)rep(paste0("Local:",names(locallist)[k]),localeachlen[k]),simplify=F)),
			Evidence="local", Ontology="local",symbol=localunlist),stringsAsFactors=F))
		set2eg <- rbind(newlist,set2eg)
		if(set_id=="path_id")names(newlist)[2] <- set_id
		}
		set2eg <- switch(idtype, ENTREZID=set2eg[set2eg$gene_id %in% names(scores),],
			SYMBOL=set2eg[set2eg$symbol %in% names(scores),]) # only genes in scores
		
		setsize0 <- tapply(set2eg[[set_id]],set2eg[[set_id]],length)
		setsizeV <- setsize0[set2eg[[set_id]]]
		setlarge <- which(setsizeV > minsetsize)
		set2eg <- set2eg[setlarge, ]

		scores0 <- scores
		scores <- scores0[which(names(scores0)%in%set2eg[[idtype_name]])]

  ## SCORES ##

	EMout <- eacifun(scores, set2eg[[idtype_name]], set2eg[[set_id]], iter=iter)
	mean0 <- EMout[[1]]
	meant <- EMout[[4]]
	sd0 <- EMout[[2]]
	sdt <- EMout[[5]]
	tau0 <- EMout[[3]]
	taut <- EMout[[6]]
	alphat <- EMout[[7]]
	betat <- EMout[[8]]
	pval <- EMout$pval
	meanbg <- EMout$meanbg
	sdbg <- EMout$sdbg
	t2 <- proc.time()[3]

  ## Gene Scores and Annotation ##
  set2eg <- unique(set2eg[,c(set_id, 'gene_id', 'symbol')])
  set.data <- set2eg
  
  set.mean <- mean0 
  set.sd <- sd0 
  set.size <- tapply(set.data[[set_id]],set.data[[set_id]],length)

 
  res <- data.frame(set.mean=set.mean, set.sd=set.sd,
                      set.size=set.size, pval=pval)
  res <- res[order(abs(res$set.mean),decreasing=T),]                    

  aux <- list(set.data = set.data, meantrace=meant, sdtrace=sdt,tautrace=taut,
							alphatrace=alphat, betatrace=betat, meanbg=meanbg, sdbg=sdbg)
  
  if (!annotate) {
    main <- res
  }
  if (annotate) {
    message("Labeling output ...")
    fn_loadSetLibraries( sets=sets )
    if (sets == "GO") {
      gterms <- toTable(GOTERM)
      main <- data.frame(gterms[match(rownames(res),gterms$go_id),
                                c("Term","Ontology")],res)
      rownames(main) <- rownames(res)
    }
    if (sets == "KEGG") {
      kterms <- toTable(KEGGPATHID2NAME)
      main <- data.frame(kterms[match(rownames(res),kterms$path_id),
              "path_name",drop=FALSE], res)
      rownames(main) <- rownames(res)
    }
  }
  out <- list(setscores = main, aux = aux, call = match.call())
  return(out)
}
