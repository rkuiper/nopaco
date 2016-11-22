#' @importFrom parallel detectCores
.onLoad <- function(libname, pkgname) {
    options("concordance.nCPU" = parallel::detectCores(),"concordance.nDraws"=1e4,"concordance.seed"=1,"concordance.verbose"=TRUE)
    library.dynam("nopaco",pkgname,libname)
}

.onunLoad<-function(libname, pkgname){
    library.dynam.unload("nopaco",pkgname,libname)
}

#' @useDynLib nopaco
NULL
