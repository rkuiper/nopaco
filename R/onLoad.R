#' @importFrom parallel detectCores
.onLoad <- function(libname, pkgname) {
    options(
    "nopaco.nCPU"    = 2, #CRAN policy: "If running a package uses multiple threads/cores it must never use more than two simultaneously: the check farm is a shared resource and will typically be running many checks simultaneously. "
    "nopaco.nDraws.CI"  = 1e4,
    "nopaco.seed"    = 1,
    "nopaco.verbose" = TRUE)
    library.dynam(pkgname,pkgname,libname)
}

.onunLoad<-function(libname, pkgname){
    library.dynam.unload(pkgname,pkgname,libname)
}

#' @useDynLib nopaco, .registration = TRUE 
NULL
