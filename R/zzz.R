#' @rdname getPsi-methods
#' @export
setMethod(
    f= "getPsi",
    signature= signature(x="matrix",y="missing"),
    function(x,y,...){
      getPsi(x=x,y=NULL,...)
    }
)

#' @rdname getPsi-methods
#' @export
setMethod(
    f= "getPsi",
    signature= signature(x="data.frame",y="missing"),
    function(x,y,...){
      getPsi(x=as.matrix(x),y=NULL,...)
    }
)
#' @rdname getPsi-methods
#' @export
setMethod(
    f= "getPsi",
    signature= signature(x="data.frame",y="data.frame"),
    function(x,y,...){
      getPsi(x=as.matrix(x),y=as.matrix(y),...)
    }
)

#' @rdname getPsi-methods
#' @export
setMethod(
    f= "getPsi",
    signature= signature(x="matrix",y="NULL"),
    function(x,y,...){
    	dots<-list(...);
      doCheck = TRUE;
      if (!is.null(dots$doCheck)){doCheck = dots$doCheck;}
      if (doCheck == TRUE) {
        standardChecks <- .doStandardChecks(x=x);
        x <- standardChecks$x
      }
    	psi<-.Call("getPsi202",x)
      return(psi)
    }
)

#' @rdname getPsi-methods
#' @export
setMethod(
    f= "getPsi",
    signature= signature(x="matrix",y="matrix"),
    function(x,y,...){
      dots<-list(...);
      doCheck = TRUE;
      if (!is.null(dots$doCheck)){doCheck = dots$doCheck;}
      if (doCheck == TRUE) {
        standardChecks <- .doStandardChecks(x=x,y=y);
        x <- standardChecks$x
        y <- standardChecks$y
      }

      psi1<-.Call("getPsi202",x)
      psi2<-.Call("getPsi202",y)
      return(c(psi1,psi2))
    }
)



getVar<-function(mat,...){
    dots<-list(...)
    doCheck = TRUE
    if (!is.null(dots$doCheck)){doCheck = dots$doCheck}
    if (doCheck == TRUE) {
        standardChecks <- .doStandardChecks(x=x,y=NULL);
        x <-standardChecks$x
    }
    bn<-rowSums(!is.na(x))
    t<-sum(bn)
    omega<-getOmega(bn)
    cij<-0
    for (i in c(1:(length(bn)-1))){
    for (j in c((i+1):length(bn))){
        cij<-cij+2*(bn[i]-1)*bn[i]*(bn[j]-1)*bn[j]
    }}
    cii<-0
    for (i in c(1:length(bn))){
        cii<-cii+(bn[i]-1)*bn[i]*(bn[i]+3)*(t-bn[i])
    }
    return((cii-cij)/(45*omega^2)*(t+1))
}

.bootstrapCI<-function(x,y=NULL){ 
    ##Obtain confidence interval by bootstrapping
	if (exists(".Random.seed", .GlobalEnv)) {
        oldseed <- .GlobalEnv$.Random.seed
    } else {
        oldseed <- NULL
    }

	set.seed( as.integer(options("nopaco.seed")$nopaco.seed) )
    psi<-.Call(
        "bootstrapCI",
        x,
        y,
        as.integer(options("nopaco.nDraws.CI")$nopaco.nDraws.CI),
        as.integer(options("nopaco.nCPU")$nopaco.nCPU)
    )


	if (!is.null(oldseed)) {
        .GlobalEnv$.Random.seed <- oldseed
    } else {
        rm(".Random.seed", envir = .GlobalEnv)
    }
    return(psi)
}

exactProb<-function(bn, ...){
    dotList<-list(...)
    skipTests<-dotList[["skipTests"]]
    validDotArguments<-c("skipTests")
    if (!all(names(dotList)%in%validDotArguments)){
        stop(paste("Unknown arguments given to exactProb function. Valid are:",
            paste(validDotArguments,collapse=", "))
        )
    }
    if (is.null(skipTests)){skipTests<-FALSE}
    storage.mode(bn) <- "integer"
    res<-.Call(
        'exactDistr202',
        bn,
        as.logical(skipTests))

    return(res)
}

#' @aliases rfromPsi
#' @title Convertion between Pearson correlation and the non paramtric concordance coefficient
#'
#' @description Convertion between Pearson correlation and the non paramtric concordance coefficient
#'
#' @param psi a (vector of) non paramtric concordance coefficient(s)
#'
#' @return A (vector of) corresponding Pearson correlation coefficient(s).
#'
#' @details
#' The convertion is performed following the relationship described by Rothery (1979).
#' \code{2*cos(pi*(1-psi))-1}
#'
#' @references Rothery, P. 'A nonparametric measure of intraclass correlation', Biometrika, 66, 3, 629-639  (1979).
#'
#'
#' @family concordance functions
#'
#' @docType methods
#' @rdname convertionMethods
#'
#' @examples
#'
#' #Generate a matrix without concordance
#' matRandom <- matrix(rnorm(30),10,3)
#' result<-concordance.test(matRandom)
#' getPsi(result) #concordance coefficient
#' result$ci      #95% confidence interval
#'
#' #Corresonding Pearson correlation
#' rfromPsi(getPsi(result))
#' rfromPsi(result$ci)
#'
#' @export
rfromPsi<-function(psi){
    psi<-sapply(sapply(psi,max,0.5),min,1)
    2*cos(pi*(1-psi))-1
}

#' @aliases psifromR
#'
#' @param r a (vector of) Pearson correlation coefficient(s)
#' @docType methods
#' @rdname convertionMethods
#'
#' @examples
#' #Plot the relation between Pearson correlation and the nonparamatric concordance coefficient.
#' r<-seq(-1,1,0.01)
#' psi<-psifromR(r)
#' plot(r,psi,type='l',xlab="Pearson correlation", ylab="nonparametric concordance")
#'
#' @export
psifromR<-function(r){
    r<-sapply(sapply(r,max,-1),min,1)
    1-acos((r+1)/2)/pi
}

.doStandardChecks<-function(x,y=NULL){
		#Checks if nopaco.nCPU<=64 : MAXIMUM_WAIT_OBJECTS must not be > 64, 
		#Checks if nopaco.nDraws.CI>=1000
		
    #Adds numeric rownames if rownames are absent
    #Adds numeric colnames if colnames are absent
    #Stops in case of non numeric values (NA's are allowed)

    ##If y is given:
    ###Stops if rownames between x and y are not similar
    ###Warns if colnames between x and y are not similar

    ##Excludes observations with less than two replicates in either x or y
    #Sets any NA or NaN explicitly to NaN which is used in external C code.
    
    if (options("nopaco.nCPU")> 64 ) stop("options('nopaco.nCPU') must be <= 64")
    if (options("nopaco.nCPU")<1   ) stop("options('nopaco.nCPU') must be >0")
    if (options("nopaco.nDraws.CI")<1000   ) stop("options('nopaco.Draws.CI') must be >=1000")
    
    if(is.data.frame(x)){
        if (!all(unlist(lapply(x,class))=="numeric")) stop("x contains non numeric values")
    }
    if(is.data.frame(y)){
        if (!all(unlist(lapply(y,class))=="numeric")) stop("y contains non numeric values")
    }
    x<-as.matrix(x)
    y<- if (is.null(y)) {NULL } else { as.matrix(y) }

    if( is.null(rownames(x))) rownames(x)<-seq_len(nrow(x))
    if( is.null(colnames(x))) colnames(x)<-seq_len(ncol(x))
    x[is.na(x)]<-NaN ##change e.g. NaN in NA

    #Check compatability between x and y
    if (!is.null(y) ){
		  if(is.null(rownames(y))) rownames(y)<-c(1:nrow(y))
      if(is.null(colnames(y))) colnames(y)<-c(1:ncol(y))
      y[is.na(y)]<-NaN ##change e.g. NaN in NA

  
      if ( all(dim(x)==dim(y)) ){
          rnames1<-rownames(x)
          rnames2<-rownames(y)
          if (any(rnames1!=rnames2) ) {stop("Rownames of matrices are not in the same order!")}
      }
      else { stop("Matrix dimensions are not comparable!") }

      #Check similarity between matrices column names
      cnames1<-colnames(x)
      cnames2<-colnames(y)
      l<-length(intersect(cnames1,cnames2))
     	if (l!=length(cnames1) | l!=length(cnames2) ){
      	if (ncol(x)!=ncol(y)) stop("Different numbers of columns between matrices")
      	warning("Colnames of matrices are not comparable. Going to assume that corresponding columns are in the same order!")
     	} else {
          y<-y[,cnames1,drop=FALSE];
      }
      
      if (any(is.na(x)!=is.na(y))){
          warning("Excluding some measurements to make matrices comparable!")
          y[!is.na(x)]<-NaN
          x[!is.na(y)]<-NaN
      }
    }
    
    ##Exclude observations with less than two replicates in either x or y
    lessThanTwo<-rowSums(!is.na(x))<2
    if (!is.null(y) ) lessThanTwo<-lessThanTwo | rowSums(!is.na(y))<2
    if (any(lessThanTwo)) {
     	warning(paste0("Exluding ", sum(lessThanTwo), " measurements with <2 replicates. Remaining number of observations: ", sum(lessThanTwo==FALSE)))
	    x<-x[which(lessThanTwo==FALSE),,drop=FALSE]
      if (!is.null(y) ) y<-y[which(lessThanTwo==FALSE),,drop=FALSE]
    }

	  ##Check for ties
    numOfTies.x<-0
    numOfTies.y<-0
    t1<-table(rank(x)) ##First rank, otherwise numbers close to machine precision will be seen as ties
    numOfTies.x<-sum(t1[which(t1>1)])
		if (!is.null(y) ) {
		  t2<-table(rank(y)) ##First rank, otherwise number close to machine precision will be seen as ties
  	  numOfTies.y<-sum(t2[which(t2>1)])
	  }

    return(list(x=x,y=y,ties.x=numOfTies.x,ties.y=numOfTies.y))
}

getOmega<-function(bn){
    sum(bn*(bn-1)*(sum(bn)-bn))
}


#' @importFrom stats qbeta
.minPsi<-function(bn){
    bn<-sort(bn,decreasing=TRUE)
    maxB<-max(bn)
    n<-length(bn)
    Q<-R<-matrix(c(1:(maxB*n)),nrow=n)
    for (i in seq_along(bn)){
        R[i,-c(0:bn[i])]<-NA
    }
    Q[]<-rank(R)
    Q[is.na(R)]<-NA
    getPsi(Q)
}

#' @importFrom stats qbeta
.confEstimatorBeta<-function(x,mu,p,lower.tail,targetValue){
    b<-x
    a<-mu*b/(1-mu)
    (qbeta(p,a,b,lower.tail=lower.tail)-targetValue)^2
}


#' @importFrom stats optimize pbeta qbeta
.rBetaApproximation<-function(x, alpha){
		p<-ci.lower<-NA_real_
    psi<-getPsi(x)
    if (!is.na(psi)) {
		  v<-getVar(x,doCheck=FALSE) ##Get variance given H0(psi=2/3)
		  bn<-rowSums(!is.na(x))
		  meanB<-sum(bn*(bn*(bn-1)))/sum(bn*(bn-1))


		  zeta<-(2/3 - sqrt(meanB+1)/(9/2*(meanB-1)^1.5))
		  iii<-2/getOmega(bn)
		  alphaPar<-(4*(meanB+1)/(81*(meanB-1)^3)-((8/9)*(meanB+1)^1.5)/(81*(meanB-1)^4.5))/v-sqrt(4*(meanB+1)/(81*(meanB-1)^3))
		  betaPar<-alphaPar*(9*(meanB-1)^1.5/(2*sqrt(meanB+1))-1)
		  p<-pbeta(psi-zeta-iii,shape1=alphaPar,shape2=betaPar,lower.tail=FALSE)
		  est<-optimize(interval=c(1,1e8),mu=psi-zeta-iii,p=p,lower.tail=TRUE,targetValue=2/3-zeta-iii, f=.confEstimatorBeta)
		  betaEst<-est$minimum
		  alphaEst<-(psi-zeta-iii)*betaEst/(1-(psi-zeta-iii))
		  ci.lower=max(.minPsi(bn),qbeta(alpha,alphaEst,betaEst,lower.tail=TRUE)+zeta+iii)
	  }
    list(p = p, ci.lower = ci.lower, ci.upper=1, ci.method="Rbeta")
}

.betaApproximation<-function(x, alpha){
		p<-ci.lower<-NA_real_
    psi<-getPsi(x)
	  if (!is.na(psi)) {
		  v<-getVar(x,doCheck=FALSE) 
		  bn<-rowSums(!is.na(x))
		  m<-2/3

		  if (v>=(m*(1-m))){stop("Cannot approximate by a beta distribution due to overdispersion of the estimates")}
		  alphaPar<-m*(m*(1-m)/v-1)
		  betaPar<-(1-m)*(m*(1-m)/v-1)

		  p<-pbeta(psi,shape1=alphaPar,shape2=betaPar,lower.tail=FALSE)
		  est<-optimize(interval=c(1,1e8),mu=psi,p=p,lower.tail=TRUE,targetValue=m, f=.confEstimatorBeta)
		  betaEst<-est$minimum
		  alphaEst<-psi*betaEst/(1-psi)
      ci.lower = max(.minPsi(bn), qbeta(alpha,alphaEst,betaEst,lower.tail=TRUE) )
		  }
	  list(p = p, ci.lower = ci.lower, ci.upper=1, ci.method="beta")
}

.normalApproximation<-function(x, alpha){
		p<-ci.lower<-NA_real_
	  psi<-getPsi(x)
	  if (!is.na(psi)) {
		  v<-getVar(x,doCheck=FALSE) 
		  ci.lower<-ci.upper<-NA

		  p<-pnorm(q=psi,mean=2/3,sd=sqrt(v),lower.tail=FALSE)
		    if (p>(alpha/5)) { ##otherwise perform estimate by bootstrap (see below)
		      ci.lower<-qnorm(alpha,mean=psi,sd=sqrt(v),lower.tail=TRUE)

		  }
		}
    list(p = p, ci.lower = ci.lower, ci.upper=1, ci.method="normal")
}

getAgreeMat<-function(mat){
    am<-matrix(2/3,ncol=ncol(mat),nrow=ncol(mat))
    n<-matrix(1,ncol=ncol(mat),nrow=ncol(mat))
    for (i in c(1:(ncol(mat)-1))){
    for (j in c((i+1):(ncol(mat)))){
		#cat("i=",i," j=",j,"\n");
        n[i,j]<-n[j,i]<-sum(rowSums(!is.na(mat[,c(i,j)]))==2)
        if (n[i,j]>1){
			tmpMat<-mat[rowSums(!is.na(mat[,c(i,j)]))==2,c(i,j)]
            am[i,j]<-am[j,i]<-getPsi(tmpMat)
        }
    }}
    diag(am)<-1
    diag(n)<-apply(mat,2,function(x){sum(!is.na(x))})
    list(mat=am,n=n)
}
