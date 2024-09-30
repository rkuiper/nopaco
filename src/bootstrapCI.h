//----------------------------------------------------------------
// Name        : bootstrapCI.h
// Author      : Rowan Kuiper
// Version     : 1.0.1
// Copyright   :
// Description : concordance
//----------------------------------------------------------------

#ifndef __BOOTSTRAPCI_H__
#define __BOOTSTRAPCI_H__


#define STRICT_R_HEADERS
#define R_NO_REMAP
//----------------------------------------------------------------
#include <R.h>
#include <Rinternals.h>
//----------------------------------------------------------------


extern "C" {
	SEXP bootstrapCI(SEXP MAT1,SEXP MAT2,SEXP rP_nDraws,SEXP rP_nCPU);
}

#endif /*__BOOTSTRAPCI_H__*/
