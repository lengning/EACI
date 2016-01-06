#' @title EACI model for gene set enrichment analysis
#' @usage eacifun(X, GeneV, SetV, meanIn=NULL, sdIn=NULL, tauIn=NULL, alpha0=4, beta0=1/4,iter=2,
#' pred.dist=T, meanbg=NULL)
#' @param X numeric vector for gene scores.
#' @param GeneV gene names. The order of GeneV should match X.
#' @param SetV gene set names. The order of SetV should match X. If one gene belongs to
#' two sets, X and GeneV should have two duplicated duplicated entries for this gene.
#' In SetV, these two duplicated entries should have different set names.
#' @param meanIn,sdIn,tauIn  hyper parameters mean, sd, tau. Default is NULL. If NULL, they will be 
#' estimated using empirical data
#' @param alpha0,beta0 hyper parameter for prior
#' @param iter number of iterations
#' @param pred.dist posterior predictive distribution? If not, prior predictive distribution will be used.
#' @param meanbg background mean. If it is NULL, it will be estimated using empirical data.
#' @author Ning Leng
#' @return estimated mean sd and tau
#' @examples scores <- rnorm(20)

eacifun <- function(X, GeneV, SetV, meanIn=NULL, sdIn=NULL, tauIn=NULL, alpha0=4, beta0=1/4,iter=2,
									pred.dist=T, meanbg=NULL){
# X gene scores
# GeneV a vector of characters indicating gene name
# SetV a vector of characters indicating sets, order should be matched with geneV
# meanIn sdIn mean and sd for each set. names should be set names
# tauIn tau's for each gene in each set

Len <- c( length(GeneV), length(SetV))	
if(length(unique(Len))>1)stop("Different lengths between GeneV and SetV!")
if(length(setdiff(GeneV, names(X)))>0) stop("Sets in GeneV but not in X names!")
if(length(setdiff(names(X), GeneV))>0) stop("Sets in X names but not in GeneV!")

# reorder - group by set	
names(SetV)=names(GeneV)=paste(GeneV, SetV,sep="_")
GeneOrder=order(GeneV)
SetV2 <- SetV[GeneOrder]
GeneV2 <- GeneV[GeneOrder]
SetOrder <- order(SetV2)
GeneO <- GeneV2[SetOrder]
SetO <- SetV2[SetOrder]
XO <- X[GeneO]


Allsets <- unique(unique(SetO))
Allgenes <- unique(GeneO)
Size <- tapply(SetO, SetO, length)
NumSet <- tapply(GeneO, GeneO, length)

if(length(setdiff(names(meanIn), Allsets))>0) stop("Sets in SetV but not in meanIn names!")
if(length(setdiff(names(sdIn), Allsets))>0) stop("Sets in SetV but not in sdIn names!")

if(is.null(meanIn)){
	NumSetO <- NumSet[GeneO]
	meanIn <- tapply(XO/NumSetO, SetO, mean)}
if(is.null(sdIn))sdIn <- tapply(XO, SetO, function(i)1)


if(is.null(tauIn))tauIn <- unlist(sapply(1:length(Size),function(i){tmp=rep(1/Size[i],Size[i])}))
if(length(tauIn)!=length(GeneV)) stop("length of tauIn should match GeneV!")

t1 <- proc.time()[3]

meanO <- meanIn[SetO]
sdO <- sdIn[SetO]
SizeO <- Size[SetO]

if(!pred.dist)FList <- dnorm(XO, meanO, sdO,log=TRUE )+dgamma(1/sdO^2, shape=alpha0, rate=beta0,log=T)
if(pred.dist){
# posterior predictive
alt <- alpha0 +.5*SizeO
bet <- beta0 + .5*(sdO^2)*SizeO
FList <- dt((XO-meanO)*sqrt(alt)/sqrt(bet), df=2*alt,log=T)
}
#FList <- alpha0*log(beta0) - log(gamma(alpha0)) + log(gamma(alt)) - alt*log(bet) - .5*SizeO*log(2*pi)
MultiList <- FList + log(tauIn)
Nu <- 400
MultiListNuexp <- exp(MultiList+Nu)
#alphaList0 <- tapply(MultiListNuexp, SetO, mean) 
#Max <- max(alphaList0[which(!alphaList0%in%c(NA, Inf))])
#alphaList0[which(alphaList0==Inf)]=Max
#alphaList0[which(is.na(alphaList0))]=0
#alphaList <- alphaList0/sum(alphaList0) # mean of taus; sum alpha = 1
#alphaO <- alphaList[SetO]
#alphaO0 <- alphaList0[SetO]
#taunewList0 <- unlist(tapply(alphaO0, GeneO, function(i)i/sum(i)))
taunewList0 <- unlist(tapply(MultiListNuexp, GeneO, function(i)i/sum(i)))
names(taunewList0)<- unlist(tapply(names(GeneO), GeneO, c))
	

taunewList <- taunewList0[names(SetO)]

genetausum <- tapply(taunewList, GeneO, sum)
genetausumO <- genetausum[GeneO]
meanrawList <- XO * taunewList / genetausumO
meanNew <- tapply(meanrawList, SetO, mean)

MinusList <- (XO-meanNew[SetO])^2
sdrawList <- sqrt(MinusList * taunewList / genetausumO)
sdNew <- tapply(sdrawList, SetO, mean)
meanNew[which(is.na(meanNew))]=0
sdNew[which(is.na(sdNew)| sdNew==0)] <- min(sdNew[which(sdNew>0 & !is.na(sdNew))])

var0 <- sdNew^2
G_mu <- mean(var0)
G_var <- var(var0)
G_a <- G_mu^2 / G_var + 2
G_b <- G_mu * (G_a-1)
if(!is.na(G_a))alpha0 <- G_a
if(!is.na(G_b))beta0 <- G_b

meant = sdt = taut = alphat = betat = NULL
meant <- rbind(meant, meanNew)
sdt <- rbind(sdt, sdNew)
taut <- rbind(taut, taunewList)
alphat <- c(alphat, alpha0)
betat <- c(betat, beta0)
t2 <- proc.time()[3]
  cat(paste("iteration", 1, "done; time ", round(t2-t1,2), "sec \n"))

if(iter==1)Out=list(meanNew=meanNew, sdNew=sdNew, tauNew=taunewList)
else{
for(ii in 2:iter){
t1 <- proc.time()[3]
meanO <- meanNew[SetO]
sdO <- sdNew[SetO]
tauIn <- taunewList
LogtauIn <- log(tauIn)
Minlogtau <- min(LogtauIn[which(LogtauIn!=-Inf)])
LogtauIn[which(LogtauIn==-Inf)] <- Minlogtau


if(!pred.dist)FList <- dnorm(XO, meanO, sdO,log=TRUE )+dgamma(1/sdO^2, shape=alpha0, rate=beta0,log=T)
if(pred.dist){
#posterior predictive
alt <- alpha0 +.5*SizeO
bet <- beta0 + .5*sdO^2*SizeO
FList <- dt((XO-meanO)*sqrt(alt)/sqrt(bet), df=2*alt, log=T)
}
#FList <- alpha0*log(beta0) - log(gamma(alpha0)) + log(gamma(alt)) - alt*log(bet) - .5*SizeO*log(2*pi)
MultiList <- FList + LogtauIn
MultiListNuexp <- exp(MultiList+Nu)
#alphaList0 <- tapply(MultiListNuexp, SetO, mean)

#Max <- max(alphaList0[which(!alphaList0%in%c(NA, Inf))])
#alphaList0[which(alphaList0==Inf)]=Max
#alphaList0[which(is.na(alphaList0))]=0
#alphaO0 <- alphaList0[SetO]
#alphaList <- alphaList0/sum(alphaList0) # mean of taus; sum alpha = 1
#alphaO <- alphaList[SetO]
taunewList0 <- unlist(tapply(MultiListNuexp, GeneO, function(i)i/sum(i)))
names(taunewList0)<- unlist(tapply(names(GeneO), GeneO, c))
taunewList <- taunewList0[names(SetO)]

genetausum <- tapply(taunewList, GeneO, sum)
genetausumO <- genetausum[GeneO]
meanrawList <- XO * taunewList / genetausumO
meanNew <- tapply(meanrawList, SetO, mean)

MinusList <- (XO-meanNew[SetO])^2
sdrawList <- sqrt(MinusList * taunewList / genetausumO)
sdNew <- tapply(sdrawList, SetO, mean)
meanNew[which(is.na(meanNew))]=0
sdNew[which(is.na(sdNew)| sdNew==0)] <- min(sdNew[which(sdNew>0 & !is.na(sdNew))])

var0 <- sdNew^2
G_mu <- mean(var0)
G_var <- var(var0)
G_a <- G_mu^2 / G_var + 2
G_b <- G_mu * (G_a-1)
if(!is.na(G_a))alpha0 <- G_a
if(!is.na(G_b))beta0 <- G_b


meant <- rbind(meant, meanNew)
sdt <- rbind(sdt, sdNew)
taut <- rbind(taut, taunewList)
alphat <- c(alphat, alpha0)
betat <- c(betat, beta0)
t2 <- proc.time()[3]
cat(paste("iteration", ii, "done; time ", round(t2-t1,2), "sec \n"))
}


#PP
#meanO <- meanNew[SetO]
#sdO <- sdNew[SetO]

#if(!pred.dist)P_enrichO <- dnorm(XO, meanO, sdO,log=TRUE )+dgamma(1/sdO^2, shape=alpha0, rate=beta0,log=T)
#if(pred.dist){
#posterior predictive
#	alt <- alpha0 +.5*SizeO
#	bet <- beta0 + .5*sdO^2*SizeO
#	P_enrichO <- dt((XO-meanO)*sqrt(alt)/sqrt(bet), df=2*alt, log=T)
#}
#P_bgO <- dnorm(XO, 0, sdO,log=T)
#P_enrich <-tapply(P_enrichO, SetO, sum)
#P_bg <-tapply(P_bgO, SetO, sum) 

# Z
if(is.null(meanbg))meanbg <- sum(X)/length(XO)
sdbg <- mean(sdNew)
#prb <- pnorm(meanNew, meanbg, sdbg) 
prb <- pnorm(meanNew, meanbg, sd(meanNew)) 
pout <- ifelse(1-prb>prb, prb, 1-prb)*2

MultiList <- FList + LogtauIn


Out=list(meanNew=meanNew, sdNew=sdNew, tauNew=taunewList,
				 meantrace=meant, sdtrace=sdt,tautrace=taut, alphatrace=alphat,
				 betatrace=betat, pval=pout, meanbg=meanbg, sdbg=sdbg)#, PP=c(P_enrich, P_bg))
}

Out
}


