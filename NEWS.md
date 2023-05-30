
# NEWS file for the nopaco package

## Version 1.0.7

* The default nopaco.nCPU setting previously set by parallel::detectCores(), is now maximized to a value of 2. If higher values are required, it can be set manually. 
* Count ties in both matrices (if given), instead of only matrix 'x'.
* Can now handle infinite values in the input matrices. Previously Inf and -Inf were set to NA.
* Changed the direction of 'delta' in the difference test.
* The options `concordance.nCPU`, `concordance.seed` and `concordance.verbose`, are replaced by `nopaco.nCPU`, `nopaco.seed` and `nopaco.verbose` respectively. The option `concordance.nDraws` is replaced by `nopaco.nDraws.CI` (for CI estimation).
* In the *concordance.test* function, method = `simulation` has been renamed as `bootstrap`.

## Version 1.0.6

* Fixed memory protection inconsitentsies in c++ code that were introduced in the last bugfix

## Version 1.0.5

* In order to resolve the linking error on Solaris, the c++ headers are reordered in line with the description in 'Writing R Extensions': R headers after system headers.
* Fixed memory protection inconsitentsies in c++ code
* Added ORCID to DESCRIPTION file
* Updated CITATIONS file

## Version 1.0.4

* Avoids non-stable CI estimates by excluding replicates with less than 3 values in the concordance test for differences.
* Now obtains confidence intervals by bootstrapping instead of assuming asympototic normality.
* Debuged code: Confidence intervals were not correct if (psi1 > psi2) in the concordance.test function
* Removed the Latex package headers for language and font encodings: 'english', 'utf8x', and 'T1'

## Version 1.0.3

* Added a vignette explaining the basic use of the package. The vignette works out an expample which data is given in the newly added data scores.rda (also requiring newly added data documentation)
* A major change was made by replacing the Google dense-hash (and corresponding configuration script) by an c++11 std unordered_map in order to solve the poor portability of the package. 
* According to rcheck tool by Thomas Kalibera (see github.com/kalibera/rchk) the samplePsi.cpp file did contain some unprotected variables which are now protected.
* Concordance.test could sometimes give p values slightly larger (in the order of machine precision) than 1. This is resolved.
* Implementation of tied values is handled in this version.
* The .samplePsi function did not correctly release the random.seed within the global.env causing random number generators to be non random anymore. 
* Fixed a warning (Found no call to: ‘R_useDynamicSymbols’) by adding a call to  R_useDynamicSymbols(info, TRUE); in register.cpp

## Version 0.99.8
    
* The previous version still did not compile on solaris (and possible also OSX). I found the cause in the configuration file which was not enforcing c++11. This is fixed in this version.


## Version 0.99.7

* The previous version contained memory leaks. These have been corrected by:
   	replacing a delete by free (in getPsi.cpp:68)
   	freeing pointers pStateLimits and pState (in exactdistribution202.cpp:296-297)
   	No more leaks are found by the valgrind tool.
* After submission to CRAN, Kurt Hornik reported a failure to install on Debian Linux. Although the error could not be reproduced, it is assumed to be related to the c++11 dependency. The current version includes 'CXX_STD = CXX11' in the Makevars file in order to force compilation and linking with a C++11 compiler.
	* The previous nopaco package did not compile on some systems because of non supported c++ headers. These headers were used in the external code (sparsehash @google). The current version includes a configure script that came with the sparsehash code hopefully solving the issue.
    * remove some non-used files shipped with the external sparsehash code.
    * in addition to  the 'cph' role for google, I added 'aut' in the description file.
    * more specifically identified the files to which the copyright applies

## Version 0.99.4

* First version submitted to CRAN

