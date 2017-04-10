//----------------------------------------------------------------
// Name        : register.cpp
// Author      : Rowan Kuiper
// Version     : 1.0.1
// Copyright   :
// Description : concordance
//----------------------------------------------------------------

#include "getPsi.h"
#include "samplePsi.h"
#include "exactdistribution.h"
#include <R_ext/Rdynload.h>
//----------------------------------------------------------------

static const R_CallMethodDef callMethods[] = {
{"getPsi202",  (DL_FUNC) &getPsi202, 1},
{"exactDistr202",(DL_FUNC) &exactDistr202,2},
{"samplePsi",(DL_FUNC) &samplePsi,5},
{NULL, NULL, 0}
};
//----------------------------------------------------------------


extern "C"{
    void R_init_nopaco(DllInfo *info)
    {
        R_registerRoutines(info,NULL,callMethods,NULL, NULL);
		R_useDynamicSymbols(info, TRUE);
    }
}
//----------------------------------------------------------------

