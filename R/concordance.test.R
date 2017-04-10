singleMatMethods<-c("exact","altBeta","beta","normal","simulation")

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
        psi2 = "numeric",
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
        if (dotList[["pvalue"]]<0 | dotList[["pvalue"]]>1  ) {
            stop("pvalue must be 0 <= p <= 1")
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
        if (dotList[["psi1"]]<0 | dotList[["psi1"]]>1  ) {
            stop("psi1 must be 0<=psi1<=1")
        }
        if (!is.numeric(dotList[["psi2"]]) | length(dotList[["psi2"]])!=1 ){
            stop("psi2 must be a numeric value of length 1")
        }
        if (!is.na(dotList[["psi2"]])){
            if (dotList[["psi2"]]<0 | dotList[["psi2"]]>1  ) {
                stop("psi2 must be 0<=psi2<=1")
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
        if (! dotList[["method"]] %in% singleMatMethods ) {
            stop("method must be one of ",paste(singleMatMethods,collapse=','))
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
                cat("One sided concordance test for difference between two concordances\nAlternative hypothesis: psi1 > psi2\n")
            }
            else if (object$alternative=="less"){
                cat("One sided concordance test for difference between two concordances\nAlternative hypothesis: psi1 < 2/3\n")
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
        return(psi[!is.na(psi)])
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
            coefMat[,3]<-diff(object$psi)
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
#' \item Testing the deviation from random concordance: if only one matrix is given (i.e. argument \code{x}), its concordance will be tested for a deviation from a random concordance. The default alternative hypothesis is 'greater', thus higher concordance than would have been observed by chance. For small matrices (depending on number of replicate measurements) an exact method will be used to determine to p-value. In case of larger matrices either, the alternative beta approach (default), a beta approximation or a normal approximation is used. To enforce the use of either one method, the \code{method} argument can be used with value "exact","alt.beta","beta" or "normal".
#' \item Testing for a difference between concordances: if both arguments \code{x} and \code{y} have been given, the equality of concordances of both matrices is tested. The default alternative hypothesis is 'two.sided'. Both matrices must be of equal size and have corresponding missing entries. Eventual missing data in one matrix will also be set missing in the other.
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
#' require(MASS) ##to use mvrnorm function
#' 
#' #Generate a matrix without concordance
#' matRandom <- matrix(rnorm(3*20),20,3)
#' concordance.test(matRandom) 
#'
#' #Generate a matrix with strong concordance
#' sigma<-matrix(0.8,3,3)
#' diag(sigma)<-1
#' matConcordant <- mvrnorm(20,mu=rep(0,3),Sigma=sigma)
#' concordance.test(matConcordant)
#'
#' #Test concordances between matrices
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
    validDotArguments<-c("method","calcCI")
    if (!all(names(dotList)%in%validDotArguments)){
        stop(paste("Unknown arguments given to the 'concordance.test' function. Valid are:",paste(validDotArguments,collapse=", ")))
    }
    method<-dotList[["method"]]
    calcCI<-dotList[["calcCI"]]
    if (is.null(calcCI)){calcCI<-TRUE;}
    if (!is.logical(calcCI)){stop("calcCI must be logical");}
        
    forceExact<-FALSE
    if (is.null(method)){
        if (is.null(y)){
            method<-singleMatMethods
        } else method<-"simulation"
    } else if (method=="exact") forceExact<-TRUE
    if (is.null(y)){
        method<-match.arg(method,singleMatMethods)
    } else {
        if (method!="simulation") { stop("The method must be set to 'simulation' when testing the difference between concordances.") }
    }
    if (is.null(alternative) & is.null(y)){
        alternative=match.arg(alternative,c('greater','less','two.sided')) #Null hypothesis psi<=2/3
    } else {
        alternative=match.arg(alternative,c('two.sided','greater','less')) #Null hypothesis |psi1-psi2|==0
    }
    if (is.null(y)){
        if (alternative !="greater") warning("Only alternative='greater' is valid for testing absence of concordance.")
        alternative='greater'
    }
    ci.method <- method
    cl <- match.call()
    standardChecks <- .doStandardChecks(x=x,y=y,doCheck=TRUE)
	numOfTies<-standardChecks$ties
    x <- standardChecks$x
    y <- standardChecks$y

	if ( (numOfTies/sum(is.finite(x)))>0.05) {
		warning("Large proportion of tied values may affect test outcome.");
	}
    obs<-getPsi(x=x,y=y,doCheck=FALSE)
    psi1<-obs[1]
    psi2<-obs[2]
    p<-ci.lower<-ci.upper<-as.numeric(NA)

    .minPsi<-function(bn){
        bn<-sort(bn,decreasing=TRUE)
        maxB<-max(bn)
        n<-length(bn)
        Q<-R<-matrix(c(1:(maxB*n)),nrow=n)
        for (i in seq_along(bn)){
            R[i,-c(0:bn[i])]<-NA
        }
        Q[]<-rank(R)#,ties.method='random')
        Q[is.na(R)]<-NA
        getPsi(Q)
    }

    .confEstimatorBeta<-function(x,mu,p,lower.tail,targetValue){
        b<-x
        a<-mu*b/(1-mu)
        (qbeta(p,a,b,lower.tail=lower.tail)-targetValue)^2
    }
   
    b<-ncol(x)
    if (is.null(y)){
        if (method=="exact"){
            exactDist<-list()
            bn<- rowSums(is.finite(x))
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
            } else { stop("Should not be reached!") }
        }
        if (method != "exact"){        
            obsvar<-getVar(x,doCheck=FALSE) ##Get variance given H0(psi=2/3)
            bn<-rowSums(is.finite(x))
            meanB<-sum(bn*(bn*(bn-1)))/sum(bn*(bn-1)) #mean(bn)
            n<-length(bn)
            if (method=="altBeta"){
                x<-(2/3 - sqrt(meanB+1)/(9/2*(meanB-1)^1.5))
                iii<-2/getOmega(bn)
                alphaPar<-(4*(meanB+1)/(81*(meanB-1)^3)-((8/9)*(meanB+1)^1.5)/(81*(meanB-1)^4.5))/obsvar-sqrt(4*(meanB+1)/(81*(meanB-1)^3))
                betaPar<-alphaPar*(9*(meanB-1)^1.5/(2*sqrt(meanB+1))-1)
             
                    p<-pbeta(obs-x-iii,shape1=alphaPar,shape2=betaPar,lower.tail=FALSE)
                    est<-optimize(interval=c(1,1e8),mu=obs-x-iii,p=p,lower.tail=TRUE,targetValue=2/3-x-iii, f=.confEstimatorBeta)
                    betaEst<-est$minimum
                    alphaEst<-(obs-x-iii)*betaEst/(1-(obs-x-iii))
                    ci.lower=max(.minPsi(bn),qbeta(alpha,alphaEst,betaEst,lower.tail=TRUE)+x+iii)
                    ci.upper=1
                
            } else if (method=="beta"){
                m<-2/3
                v<-obsvar
                if (v>=(m*(1-m))){stop("Cannot approximate by beta due to overdispersion of estimates")}
                alphaPar<-m*(m*(1-m)/v-1)
                betaPar<-(1-m)*(m*(1-m)/v-1)
               
                    p<-pbeta(obs,shape1=alphaPar,shape2=betaPar,lower.tail=FALSE)
                    est<-optimize(interval=c(1,1e8),mu=obs,p=p,lower.tail=TRUE,targetValue=2/3, f=.confEstimatorBeta)
                    betaEst<-est$minimum
                    alphaEst<-obs*betaEst/(1-obs)
                    ci.upper=1
                    ci.lower=max(.minPsi(bn),qbeta(alpha,alphaEst,betaEst,lower.tail=TRUE))
                
            } else if (method=="normal"){ 
               
                    p<-pnorm(q=obs,mean=2/3,sd=sqrt(obsvar),lower.tail=FALSE)
                    if (p>(alpha/10)) { ##otherwise perform estimate by simulation (see below)
                        ci.lower<-qnorm(alpha,mean=obs,sd=sqrt(obsvar),lower.tail=TRUE)
                           ci.upper<-1
                    }
                
            } else if (method == "simulation"){ stop("simulation is not available for testing a single concordance coefficient")
            } else {    stop("Unknown method given!") }
        }#END Exact
        if (calcCI){
            if (is.na(ci.lower)|is.na(ci.upper)){
			    rmat1<-x
                rmat1[]<-rank(x)#,ties.method='random')
                rmat1[!is.finite(x)]<-NA
                am<-getAgreeMat(rmat1) 
                cormat<-am$mat
                cormat[]<-rfromPsi(am$mat)
                if (!all(eigen(cormat)$value>0)){
                    cormat<-as.matrix(nearPD(cormat,keepDiag=TRUE,maxit=10000)$mat)
                }
                SigmaEV <- eigen(cormat)
                mychol<-t(SigmaEV$vectors %*% diag(sqrt(SigmaEV$values)))
                bs<-.sampleVariance(mychol,x)
               	ci.method<-"simulation"
                if (is.na(ci.lower)) ci.lower<-quantile(bs,alpha)
                if (is.na(ci.upper)) ci.upper<-1
            }
        }
    } else if (!is.null(y)){
        bn<- rowSums(is.finite(x))
        method<-"simulation"
        ci.method<-"normal"
        rmat1<-x
        rmat1[]<-rank(x)
        rmat2<-y
        rmat2[]<-rank(y)
        rmat1[!is.finite(x)]<-NA
        rmat2[!is.finite(y)]<-NA
        mergemat<-cbind(rmat1,rmat2)
        nmat<-matrix(NA,ncol(mergemat),ncol(mergemat))
        mergemat[]<-rank(mergemat) 
        mergemat[!is.finite(cbind(rmat1,rmat2))]<-NA
		am<-getAgreeMat(mergemat)  
        cormat<-am$mat
        cormat[]<-rfromPsi(am$mat)
        
        b<-ncol(x)
        r1<-(sum(cormat[1:b,1:b]+cormat[-c(1:b),-c(1:b)])/2 - b)/(b*b-b) #Average observed correlation
        r2<-(sum(cormat[1:b,-c(1:b)]+cormat[-c(1:b),1:b])/2)/(b*b) #Average observed correlation
        par<-c(r1,r2)
        cormatH0<-matrix(par[2],2*b,2*b)
        cormatH0[1:(b),1:(b)]<-par[1]
        cormatH0[-c(1:(b)),-c(1:(b))]<-par[1]
        diag(cormatH0)<-1
        if (!all(eigen(cormatH0)$value>0)){ cormatH0<-as.matrix(nearPD(cormatH0,keepDiag=TRUE,maxit=10000)$mat)}
        SigmaEV <- eigen(cormatH0)
        mychol<-t(SigmaEV$vectors %*% diag(sqrt(SigmaEV$values)))
        
		j<-.sampleVariance(mychol,x,y)
        delta<-diff(sort(obs))
        meanN<-mean(am$n)
        obsvar<-var(j[,1]-j[,2])
        p<-NA
        if (alternative=='less') {
            p<-pnorm(delta,mean=0,sd=sqrt(obsvar) ,lower.tail=TRUE)
            if (p>(alpha/10)) { ##otherwise perform estimate by simulation (see below)
                ci.lower<-.minPsi(bn)-1 
                ci.upper<-qnorm(alpha,mean=delta,sd=sqrt(obsvar),lower.tail=FALSE)
            }
        } else if (alternative=='greater') {
            p<-pnorm(delta,mean=0,sd=sqrt(obsvar) ,lower.tail=FALSE)  
            if (p>(alpha/10)) { ##otherwise perform estimate by simulation (see below)
                ci.lower<-qnorm(alpha,mean=delta,sd=sqrt(obsvar),lower.tail=TRUE)
                ci.upper<--1*(.minPsi(bn)-1 )
            }
        } else if (alternative=='two.sided')  {
            p<-2*pnorm(abs(delta),mean=0,sd=sqrt(obsvar) ,lower.tail=FALSE)
             if (p>(alpha/10)) { ##otherwise perform estimate by simulation (see below)
                ci.lower<-qnorm(alpha/2,mean=delta,sd=sqrt(obsvar),lower.tail=TRUE)
                ci.upper<-qnorm(alpha/2,mean=delta,sd=sqrt(obsvar),lower.tail=FALSE)
            }
        }
        if (is.na(ci.lower)|is.na(ci.upper)){
            SigmaEV <- eigen( as.matrix(nearPD(cormat,keepDiag=TRUE,maxit=10000)$mat))
            mychol<-t(SigmaEV$vectors %*% diag(sqrt(SigmaEV$values)))
            samplePsi<-.sampleVariance(mychol,x,y)
            bs<-samplePsi[,2]-samplePsi[,1]
            ci.method<-"simulation"
            if (alternative=='less') {
                if (is.na(ci.lower)) ci.lower<-.minPsi(bn)
                if (is.na(ci.upper)) ci.upper<-quantile(bs,1-alpha)
            } else if (alternative=='greater') {
                if (is.na(ci.lower)) ci.lower<-quantile(bs,alpha)
                if (is.na(ci.upper)) ci.upper<-1
            } else if (alternative=='two.sided') {
                if (is.na(ci.lower)) ci.lower<-quantile(bs,alpha/2)
                if (is.na(ci.upper)) ci.upper<-quantile(bs,1-alpha/2)
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
