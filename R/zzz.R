#' @rdname getPsi-methods
#' @export
setMethod( 
    f= "getPsi", 
    signature= signature(x="matrix",y="missing"),
    function(x,y,...){
        getPsi(x=x,y=matrix(ncol=0,nrow=0))
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
        getPsi(x=x,y=matrix(ncol=0,nrow=0),doCheck=doCheck)
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
        standardChecks <- .doStandardChecks(x=x,y=y,doCheck=doCheck);
        mat <-standardChecks$x
        y<-standardChecks$y
        if (length(y)==0){
            psi<-.Call("getPsi202",mat)
            return(psi)
        } else{
            psi1<-.Call("getPsi202",mat)
            psi2<-.Call("getPsi202",y)
            return(c(psi1,psi2))
        }
    }
)

getVar<-function(mat,...){
    dots<-list(...)
    doCheck = TRUE
    if (!is.null(dots$doCheck)){doCheck = dots$doCheck}
    standardChecks <- .doStandardChecks(x=mat,y=NULL,doCheck=doCheck)
    mat <-standardChecks$x
    bn<-rowSums(is.finite(mat))
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

.sampleVariance<-function(choleskimat,missingx,missingy=NULL){
    ##Obtain variance for psi by choleskimat
    missingx<-is.finite(missingx)
    if (!is.null(missingy)){
        missingy<-is.finite(missingy)
    }
	oldseed <- .Random.seed
	set.seed( as.integer(options("concordance.seed")$concordance.seed) )
    psi<-.Call(
        "samplePsi",
        choleskimat,
        missingx,
        missingy,
        as.integer(options("concordance.nDraws")$concordance.nDraws),
        as.integer(options("concordance.nCPU")$concordance.nCPU)
    )
	.Random.seed <- oldseed
    if (is.null(missingy)){
        return(psi[,1,drop=FALSE])
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
        as.logical(skipTests),
        as.logical(options("concordance.verbose")$concordance.verbose))
    
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

.doStandardChecks<-function(x,y=NULL,doCheck=TRUE){
    #Adds numeric rownames if rownames are absent
    #Adds numeric colnames if colnames are absent
    #Stops in case of non numeric values (NA's are allowed)
    #Warns in case of ties
    ##If y is given:
    ###Stops if rownames between x and y are not similar
    ###Warns if colnames between x and y are not similar
    
    ##Excludes observations with less than two replicates in either x or y
    #Sets any non-finite values (e.g. NA's) explicitly to infinite which is used in external C code.
    if (doCheck == TRUE) {
        if(is.data.frame(x)){
            if (!all(unlist(lapply(x,class))=="numeric")) stop("x contains non numeric values")
        }
        if(is.data.frame(y)){
            if (!all(unlist(lapply(y,class))=="numeric")) stop("y contains non numeric values")
        }
        x<-as.matrix(x)
        if (length(y)==0){y<-NULL;}
        else {y<-as.matrix(y);}
        if(is.null(rownames(x))){rownames(x)<-c(1:nrow(x))}
        if(is.null(colnames(x))){colnames(x)<-c(1:ncol(x))}
        if (!is.null(y)){
            if(is.null(rownames(y))){rownames(y)<-c(1:nrow(y))}
            if(is.null(colnames(y))){colnames(y)<-c(1:ncol(y))}
        }
        x[!is.finite(x)]<-NA
        ##Check for ties
        t1<-table(rank(x)) ##First rank, otherwise number close to machine precision will be seen as ties
        t1.s<-sum(t1[which(t1>1)])
        if (t1.s>0) {warning("Found ",t1.s," ties in x")}
        #Check y compatability with x 
        if (!is.null(y) ){
            y[!is.finite(y)]<-NA
            ##Check for ties
            t2<-table(rank(y))
            t2.s<-sum(t2[which(t2>1)])
            if (t2.s>0) {warning("Found ",t2.s," ties in y")}
            if(all(dim(x)==dim(y))){
                rnames1<-rownames(x)
                rnames2<-rownames(y)
                if (any(rnames1!=rnames2) ) {stop("Rownames of matrices are not in the same order!")}
            }
            else { stop("Matrix dimensions are not comparable!")}
            #Check similarity columnnames matrices
            cnames1<-colnames(x)
            cnames2<-colnames(y)
            l<-length(intersect(cnames1,cnames2))
            if (l!=length(cnames1) | l!=length(cnames2) ){warning("Colnames of matrices are not comparable!")}
            else if (l==ncol(x) & l==ncol(y)){
                y<-y[,cnames1,drop=FALSE];
            }
            if (any(is.finite(x)!=is.finite(y))){
                warning("Excluding some measurements to make matrices comparable!")
                y[!is.finite(x)]<-NA
                x[!is.finite(y)]<-NA
            }
        }
        ##Exclude observations with less than two replicates
        if(0){###BEGIN: TEST
			isOK1<-rowSums(is.finite(x))>1
		    if (!all(isOK1)){
		        x<-x[isOK1,,drop=FALSE]
		        if (!is.null(y)){
		            y<-y[isOK1,,drop=FALSE]
		        }
		    }
		    isOK2<-TRUE;
		    if (!is.null(y)){
		        isOK2<-rowSums(is.finite(y))>1
		        if (!all(isOK2)){
		            y<-y[isOK2,,drop=FALSE]
		            if (!is.null(x)){
		                x<-x[isOK2,,drop=FALSE]
		            }
		        }
		    }
		    if ( (sum(!isOK1)+sum(!isOK2))>0){
		        warning(paste("Excluding observations with less than two known replicates (n=",sum(!isOK1)+sum(!isOK2),")",sep=""))
		    }    
		}###END: TEST
    }
    x[!is.finite(x)]<-Inf
    if (!is.null(y)) y[!is.finite(y)]<-Inf
    return(list(x=x,y=y))
}

getOmega<-function(bn){
    sum(bn*(bn-1)*(sum(bn)-bn))
}

getAgreeMat<-function(mat){
    am<-matrix(2/3,ncol=ncol(mat),nrow=ncol(mat))
    n<-matrix(1,ncol=ncol(mat),nrow=ncol(mat))
    for (i in c(1:(ncol(mat)-1))){
    for (j in c((i+1):(ncol(mat)))){
        n[i,j]<-n[j,i]<-sum(rowSums(is.finite(mat[,c(i,j)]))==2)
        if (n[i,j]>1){
            am[i,j]<-am[j,i]<-getPsi(mat[rowSums(is.finite(mat[,c(i,j)]))==2,c(i,j)])
        }
    }}
    diag(am)<-1
    diag(n)<-apply(mat,2,function(x){sum(!is.na(x))})
    list(mat=am,n=n)
}
