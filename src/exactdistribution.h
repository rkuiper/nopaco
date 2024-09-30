//----------------------------------------------------------------
// Name        : exactdistribution.h
// Author      : Remco Hoogenboezem
// Version     : 1.0.1
// Copyright   :
// Description : concordance
//----------------------------------------------------------------

#ifndef __EXACTDISTR202_H__
#define __EXACTDISTR202_H__

#define STRICT_R_HEADERS
#define R_NO_REMAP
//----------------------------------------------------------------
#include <R.h>
#include <Rinternals.h>
//----------------------------------------------------------------

extern "C"{
SEXP exactDistr202(SEXP bn,SEXP skipTests);
}

#endif /*__EXACTDISTR202_H__*/
