singleMatMethods<-c("exact","Rbeta","beta","normal","bootstrap") ##First must be "exact"
twoMatMethods<-c("bootstrap")

setClassUnion("numericOrNull", c("numeric", "NULL"))

############################################
#' Class ConcordanceTest
#'
#' Class \code{ConcordanceTest} stores results from a concordance test.
#'
#' @name ConcordanceTest-class
#' @rdname ConcordanceTest-class
#' @description  This class stores results obtained from a concordance test.
#'
#' @slot pvalue The pvalue
#' @slot psi1 The concordance in matrix \code{x}
#' @slot psi2 The concordance in matrix \code{y}
#' @slot method The method used to obtain the pvalue
#' @slot alternative The alternative hypothesis
#' @slot ci.lower The lower confidence boudary
#' @slot ci.upper The upper confidence boudary
#' @slot ci.method The method used to obtain the confidence interval
#' @slot alpha The significance level
#' @slot call The call made to the \code{\link{concordance.test}} function
#' @exportClass ConcordanceTest
setClass("ConcordanceTest",
    representation = representation(
        pvalue = "numeric",
        psi1 = "numeric",
        psi2 = "numericOrNull",
        method="character",
        alternative="character",
        ci.lower="numeric",
        ci.upper="numeric",
        ci.method="character",
        alpha="numeric",
        call="call"
    )
)


setMethod("initialize",
    signature  = signature(.Object = "ConcordanceTest"),
    definition = function(.Object, ...){
        dotList <- list(...)
        if (!is.numeric(dotList[["pvalue"]]) | length(dotList[["pvalue"]])!=1 ){
            stop("pvalue must be a numeric value of length 1")
        }
        if (!is.na(dotList[["pvalue"]])){
	        if( dotList[["pvalue"]]<0 | dotList[["pvalue"]]>1  ) {
            stop("pvalue must be 0 <= p <= 1")
          }
        }
        if (!is.numeric(dotList[["alpha"]]) | length(dotList[["alpha"]])!=1 ){
            stop("alpha must be a numeric value of length 1")
        }
        if (dotList[["alpha"]]<0 | dotList[["alpha"]]>1  ) {
            stop("alpha must be 0 < p < 1")
        }
        if (!is.numeric(dotList[["psi1"]]) | length(dotList[["psi1"]])!=1 ){
            stop("psi1 must be a numeric value of length 1. Found psi1 = ",dotList[["psi1"]])
        }
				if (!is.na(dotList[["psi1"]])){
		      if (dotList[["psi1"]]<0 | dotList[["psi1"]]>1  ) {
		          stop("psi1 must be 0<=psi1<=1")
		      }
	      }
        if (!is.null(dotList[["psi2"]])){
	        if (!is.numeric(dotList[["psi2"]]) | length(dotList[["psi2"]])!=1 ){
	            stop("psi2 must be a numeric value of length 1")
	        }
					if (!is.na(dotList[["psi2"]])){
	          if (dotList[["psi2"]]<0 | dotList[["psi2"]]>1  ) {
	              stop("psi2 must be 0<=psi2<=1")
	          }
          }
        }
        if (!is.numeric(dotList[["ci.lower"]]) | length(dotList[["ci.lower"]])!=1 ){
            stop("ci.lower must be a numeric value of length 1")
        }
        if (!is.numeric(dotList[["ci.upper"]]) | length(dotList[["ci.upper"]])!=1 ){
            stop("ci.upper must be a numeric value of length 1")
        }
        if (!is.character(dotList[["method"]]) | length(dotList[["method"]])!=1 ){
            stop("method must be a character value of length 1")
        }
        if (! dotList[["method"]] %in% c(singleMatMethods,twoMatMethods) ) {
            stop("method must be one of ",paste(c(singleMatMethods,twoMatMethods),collapse=','))
        }
        alt<-c('less','two.sided','greater')
        if (! dotList[["alternative"]] %in% alt ){
            stop("method must be one of ",paste(alt,collapse=','))
        }
        .Object@pvalue=dotList[["pvalue"]]
        .Object@psi1=dotList[["psi1"]]
        .Object@psi2=dotList[["psi2"]]
        .Object@method=dotList[["method"]]
        .Object@alternative=dotList[["alternative"]]
        .Object@ci.lower<-dotList[["ci.lower"]]
        .Object@ci.upper=dotList[["ci.upper"]]
        .Object@ci.method=dotList[["ci.method"]]
        .Object@alpha=dotList[["alpha"]]
        .Object@call=dotList[["call"]]
        return(.Object)
    }
)


setMethod("show",
    signature  = signature(object = "ConcordanceTest"),
    definition = function(object){
        if (length(getPsi(object))==2){
            if (object$alternative=="greater"){
                cat("One sided concordance test for difference between two concordances\nAlternative hypothesis: psi1 > psi2 (i.e. delta < 0)\n")
            }
            else if (object$alternative=="less"){
                cat("One sided concordance test for difference between two concordances\nAlternative hypothesis: psi1 < psi2 (i.e. delta > 0)\n")
            }
            else if (object$alternative=="two.sided"){
                cat("Two sided concordance test for difference between two concordances\nAlternative hypothesis: |psi1 - psi2| > 0\n")
            }
        } else {
            if (object$alternative=="greater"){
                cat("One sided concordance test\nAlternative hypothesis: psi > 2/3\n")
            }
            else if (object$alternative=="less"){
                cat("One sided concordance test\nAlternative hypothesis: psi < 2/3\n")
            }
            else if (object$alternative=="two.sided"){
                cat("Two sided concordance test\nAlternative hypothesis: |psi| > 2/3\n")
            }
        }
        cat("Estimates are obtained by ",paste(paste(names(object$method),": ",object$method,sep=""),collapse=" and "),"\n")
        coefs<-as.data.frame(signif(coef(object),3))
        print(coefs)
            
    }
)




##############################################################
#' Extract argument names from a ConcordanceTest object
#'
#' \code{names} extracts argument names from a \code{\link{ConcordanceTest-class}} object
#' @param x An object of \code{\link{ConcordanceTest-class}}
#'
#' @return A character vector
#'
#' @docType methods
#' @rdname names-methods
#' @aliases names,ConcordanceTest-method
#'
#' @examples
#'
#' matRandom <- matrix(rnorm(3*20),20,3)
#' testResult <- concordance.test(matRandom)
#' names(testResult)
#' @export
setMethod( "names", 
    signature  = signature("ConcordanceTest"),
    definition = function(x){
        c("pvalue","ci","psi","alternative","call","alpha","method")
    }
)

##############################################################
#' Extract argument values from a ConcordanceTest object
#'
#' Extracts argument values from a \code{\link{ConcordanceTest-class}} object
#'
#' @param x An object of \code{\link{ConcordanceTest-class}}
#' @param name The argument to get the value of
#'
#'
#' @return The value of the requested argument
#'
#' @docType methods
#' @rdname get-methods
#' @aliases get,ConcordanceTest-method
#'
#' @examples
#'
#' matRandom <- matrix(rnorm(3*20),20,3)
#' testResult <- concordance.test(matRandom)
#' names(testResult)
#' testResult$psi
#' @export
setMethod( "$", "ConcordanceTest", function(x, name ){
    name<-match.arg(name,names(x))
    if (name=="pvalue") return(x@pvalue)
    else if (name=="ci") {
        ci<-c(x@ci.lower,x@ci.upper)
        names(ci)<-c(paste("lower ",signif((1-x$alpha)*100,2),"%CI",sep=""),paste("upper ",signif((1-x$alpha)*100,2),"%CI",sep=""))
        return(ci)
    } else if (name=="psi") return(getPsi(x))
    else if (name=="alternative") return(x@alternative)
    else if (name=="call") return(x@call)
    else if (name=="alpha") return(x@alpha)
    else if (name=="method") {
        m<-c(x@method,x@ci.method)
        names(m)<-c("pvalue","confidence interval")
        return(m)
    }
    return(NULL)
} 
)

#' @rdname getPsi-methods
#' @export
setMethod( 
    f="getPsi",
    signature(x="ConcordanceTest",y="missing"),
    function(x){
        psi<-c(x@psi1,x@psi2)
        return(psi)
    }
)



##############################################################
#' Extract test results from the results of a concordance.test
#'
#' \code{coef} extract the test results from the results of a concordance.test 
#' @name coef
#' @param object An object of \code{\link{ConcordanceTest-class}}
#' @param ...  Not used
#'
#' @return A matrix
#'
#' @family concordance functions
#'
#' @export
#' @docType methods
#' @rdname coef-methods
#' @aliases coef,ConcordanceTest-method
#'
#' @examples
#'
#' matRandom <- matrix(rnorm(3*20),20,3)
#' testResult <- concordance.test(matRandom)
#' getPsi(testResult)
#' coef(testResult)
#'
#' @importFrom stats coef
#' @export
setMethod("coef",
    signature(object = "ConcordanceTest"),
    function (object, ...) 
    {
        if (length(object$psi)==2) {
            coefMat<-matrix(NA,ncol=6,nrow=1,dimnames=list("",c("psi1","psi2","delta",names(object$ci),"p")))
            coefMat[,1:2]<-object$psi
            coefMat[,3]<- -1*diff(object$psi)
            coefMat[,4:5]<-object$ci
            coefMat[,6]<-object$pvalue
        } else {
            coefMat<-matrix(NA,ncol=4,nrow=1,dimnames=list("",c("psi",names(object$ci),"p")))
            coefMat[,1]<-object$psi
            coefMat[,2:3]<-object$ci
            coefMat[,4]<-object$pvalue
        }
        return(coefMat)
    }
)

#' @aliases concordance.test nopaco
#' @title Perform a nonparametric concordance test.
#'
#' @description \code{concordance.test} performs a test for a random concordance (if a single matrix is given) or tests for equal concordance between two matrices.
#'
#' @param x a numeric matrix, subjects in the rows, repeated measurements in the columns
#' @param y (optional) a numeric matrix of equal size as argument \code{x}
#' @param alternative  "less", "greater" or "two.sided". Only used when \code{y} is given.
#' @param alpha significance level (default = 0.05)
#' @param ... see details
#'
#' @return An object of \code{\link{ConcordanceTest-class}}
#'
#' @details
#'  \itemize{
#' \item Testing the deviation from random concordance: if only one matrix is given (i.e. argument \code{x}), its concordance will be tested against the alternative hypothesis of finding a higher concordance under random sampling conditions. For small matrices (depending on number of replicate measurements) an exact method will be used to determine a p-value. In case of larger matrices, where the exact approach in not feasible, either the revised-beta approach (default), a beta approximation, or a normal approximation is used. To enforce the use of either one method, the \code{method} argument can be used with value "exact","Rbeta","beta" or "normal".
#' \item Testing for a difference between concordances: if both arguments \code{x} and \code{y} have been given, the equality of concordances of both matrices is tested. The default alternative hypothesis is 'two.sided'. Both matrices must be of equal size and have corresponding missing entries (\code{NA} values). In case of missing data in one matrix, the same entries in the other matrix will also be set to missing. This method involves bootstrapping to estimate variance and confidence estimates, with adjustable options via the global R options "nopaco.nCPU" (default=2), "nopaco.nDraws.CI" (default=1e4) and "nopaco.seed" (default = 1). 
#' }
#' Unbalanced data due to randomly missing data or an unequal number of repeated measurements per subject is allowed. In that case, missing or unknown values must be set to \code{NA}.

#'
#'  

#'
#' @family concordance functions
#'
#' @references P.Rothery (1979) Biometrika 66(3):629-639
#'
#' @docType methods
#' @rdname concordance.test
#'
#' @examples
#' require(MASS) ##to use the 'mvrnorm' function
#' 
#' #Generate a matrix without concordance
#' # for three observers in two samples
#' matRandom <- matrix(rnorm(3*20),20,3)
#' concordance.test(matRandom) 
#'
#' #Generate a corresponding matrix with strong concordance
#' sigma<-matrix(0.8,3,3)
#' diag(sigma)<-1
#' matConcordant <- mvrnorm(20,mu=rep(0,3),Sigma=sigma)
#' concordance.test(matConcordant)
#'
#' #Test concordances between the two matrices
#' aTest <- concordance.test(matConcordant, matRandom)
#' 
#' getPsi(aTest)
#' coef(aTest)
#'
#' @importFrom Matrix nearPD
#' @importFrom methods new
#' @importFrom stats cor optimize optim pbeta pnorm qbeta qnorm quantile var
#' @export
concordance.test <-
function(x,y=NULL,alternative=NULL,alpha=0.05,...){
    dotList<-list(...)
    validDotArguments<-c("method","force.ci.bootstrap")
    if (!all(names(dotList)%in%validDotArguments)){
        stop(paste("Unknown arguments given to the 'concordance.test' function. Valid are:",paste(validDotArguments,collapse=", ")))
    }
    method<-dotList[["method"]]
    force.ci.bootstrap<-dotList[["force.ci.bootstrap"]]
    

    
    if (is.null(force.ci.bootstrap)){force.ci.bootstrap<-TRUE;} ##Default: Estimate CI by bootstrapping
    if (!is.logical(force.ci.bootstrap)){stop("force.ci.bootstrap must be logical");}

    forceExact<-FALSE
    if (is.null(method)){
        method<-twoMatMethods
        if (is.null(y)){
            method<-singleMatMethods
        }  
    } else if (method=="exact") forceExact<-TRUE
    
    if (is.null(y)){
        method<-match.arg(method,singleMatMethods)
    } else {
        method<-match.arg(method,twoMatMethods)
    }
    if (is.null(y)){
        alternative=match.arg(alternative,c('greater','less','two.sided')) #Null hypothesis psi<=2/3
    } else {
        alternative=match.arg(alternative,c('two.sided','greater','less')) #Null hypothesis |psi1-psi2|==0
    }
    if (is.null(y)){
        if (alternative !="greater") warning("Using alternative='greater' which is the only valid one for testing 'absence of concordance'.")
        alternative='greater'
    }

     
    cl <- match.call()
    standardChecks <- .doStandardChecks(x=x,y=y)
		numOfTies<-standardChecks$ties.x+standardChecks$ties.y
    x <- standardChecks$x
    y <- standardChecks$y
   
		if (sum(!is.na(x))>0){
			if ( (numOfTies/sum(!is.na(x)))>0.05) {
				warning("Large proportion of tied values may affect test outcome.");
			}
		}
		
    obs<-getPsi(x=x,y=y,doCheck=FALSE)
    psi1<-obs[1]
    psi2<- if(length(obs)==2) {obs[2]} else {NULL}
    
    p<-ci.lower<-ci.upper<-NA_real_

    if (is.null(y)){
      b<-ncol(x)
      if ("exact"%in%method){
          exactDist<-list()
          bn<- rowSums(!is.na(x))
          sumN<-length(bn)
          exactDist<-try({
              exactProb(bn,skipTests=forceExact)
          },silent=TRUE)

          if (!inherits(exactDist,"try-error")){
              p<-max(0,min(1,sum(exactDist[ ((exactDist[,1])>=obs),2])))
              ci.upper<-max(exactDist[,1])
          }
          else if (!is.null(forceExact)){
              if(forceExact==TRUE) stop(exactDist)
              else {method<- singleMatMethods[2]}
          } else { stop("Runtime error: Something went wrong in estimating the exactDistribution.")  }
      }
      
      if (!"exact"%in%method){
	      if (method[1]=="Rbeta"){
            res<-.rBetaApproximation(x, alpha); 
        } else if (method[1]=="beta"){
            res<-.betaApproximation(x, alpha); 
        } else if (method[1]=="normal"){ 
            res<-.normalApproximation(x, alpha); 
        } else if (method[1] == "bootstrap"){
          stop("bootstrap is not available for testing a single concordance coefficient")
      	} else { stop("Unknown method given!") }
      	
				p<-res$p
	      ci.method<-res$ci.method
	      ci.lower<-res$ci.lower
	      ci.upper<-res$ci.upper
      } #END Single concordance methods
               
      if (force.ci.bootstrap == TRUE | is.na(ci.lower) | is.na(ci.upper)) {
           
          samplePsi<-.bootstrapCI(x=x,y=NULL)
          bs<-samplePsi[,1]
          ci.method<-"bootstrap"
          ci.lower<-quantile(bs,alpha,na.rm=TRUE)	
          ci.upper<-1
      }

    } else if (!is.null(y)){
			cs<-colSums(!is.na(x));
			includeReps<-cs>2
			if (any(includeReps==FALSE)){
				warning("Ignoring replicates with less than three values in order to stabilize variance estimation.")
			}
			
			x<-x[,includeReps, drop=FALSE]
			y<-y[,includeReps, drop=FALSE]
			
      b<-ncol(x)
		  bn <- rowSums(!is.na(x))

			observed.delta  <- -1*diff(obs) #psi1-psi2
			
      method<-"bootstrap"

      
      bootstrapped.delta<--1*apply(.bootstrapCI(x,y),1,diff)
      
      obsvar<-var(bootstrapped.delta)
      p<-NA_real_
      
      if (alternative=='less') { #Ha:psi1<psi2
        p<-pnorm(observed.delta,mean=0,sd=sqrt(obsvar) ,lower.tail=TRUE)
        if ((p>(alpha/5)) & force.ci.bootstrap==FALSE) { ##otherwise perform estimate by bootstrap (see below)
          ci.method<-"normal"
          ci.lower<-.minPsi(bn)-1 
          ci.upper<-qnorm(alpha,mean=observed.delta,sd=sqrt(obsvar),lower.tail=FALSE)
        }
      } else if (alternative=='greater') { #Ha:psi1>psi2
        p<-pnorm(observed.delta,mean=0,sd=sqrt(obsvar) ,lower.tail=FALSE)  

        if ((p>(alpha/5)) & force.ci.bootstrap==FALSE) { ##otherwise perform estimate by bootstrap (see below)
          ci.method<-"normal"
          ci.lower<-qnorm(alpha,mean=observed.delta,sd=sqrt(obsvar),lower.tail=TRUE)
          ci.upper<--1*(.minPsi(bn)-1 )
          }
      } else if (alternative=='two.sided')  {  #Ha:psi1-psi2 != 0
        p<-2*pnorm(abs(observed.delta),mean=0,sd=sqrt(obsvar) ,lower.tail=FALSE)

        if ((p>(alpha/5)) & force.ci.bootstrap==FALSE) { ##otherwise perform estimate by bootstrap (see below)
          ci.method<-"normal"
          ci.lower<-qnorm(alpha/2,mean=observed.delta,sd=sqrt(obsvar),lower.tail=TRUE)
          ci.upper<-qnorm(alpha/2,mean=observed.delta,sd=sqrt(obsvar),lower.tail=FALSE)
        }
      }
      if (force.ci.bootstrap == TRUE | is.na(ci.lower) | is.na(ci.upper)) {
          ci.method<-"bootstrap"
          if (alternative=='less') { #Ha:psi1<psi2
              ci.lower<-.minPsi(bn)-1 
              ci.upper<-quantile(bootstrapped.delta,1-alpha)
          } else if (alternative=='greater') { #Ha:psi1>psi2
              ci.lower<-quantile(bootstrapped.delta,alpha)
              ci.upper<--1*(.minPsi(bn)-1 )
          } else if (alternative=='two.sided') { #Ha:psi1-psi2 != 0
              ci.lower<-quantile(bootstrapped.delta,alpha/2)
              ci.upper<-quantile(bootstrapped.delta,1-alpha/2)
          }
      }
    }
    result<- new("ConcordanceTest",
        pvalue = p,
        psi1 = psi1,
        psi2 = psi2,
        method=method,
        alternative=alternative,
        ci.lower=ci.lower,
        ci.upper=ci.upper,
        ci.method=ci.method,
        alpha=alpha,
        call=cl
    )
    result
}
