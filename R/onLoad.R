#' @importFrom parallel detectCores
.onLoad <- function(libname, pkgname) {
    options("concordance.nCPU" = parallel::detectCores(),"concordance.nDraws"=1e4,"concordance.seed"=1,"concordance.verbose"=TRUE)
    library.dynam(pkgname,pkgname,libname)
}

.onunLoad<-function(libname, pkgname){
    library.dynam.unload(pkgname,pkgname,libname)
}

#' @useDynLib nopaco, .registration = TRUE 
NULL
