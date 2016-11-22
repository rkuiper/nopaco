#ifndef __EXACTDISTR202_H__
#define __EXACTDISTR202_H__



#define STRICT_R_HEADERS 

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>



extern "C"{
SEXP exactDistr202(SEXP bn,SEXP skipTests, SEXP verbose);
}

#endif /*__EXACTDISTR202_H__*/
