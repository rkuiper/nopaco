//----------------------------------------------------------------
// Name        : register.cpp
// Author      : Rowan Kuiper
// Version     : 1.0.1
// Copyright   :
// Description : concordance
//----------------------------------------------------------------

#include "getPsi.h"
#include "bootstrapCI.h"
#include "exactdistribution.h"
#include <R_ext/Rdynload.h>
//----------------------------------------------------------------

static const R_CallMethodDef callMethods[] = {
{"getPsi202",  (DL_FUNC) &getPsi202, 1},
{"exactDistr202",(DL_FUNC) &exactDistr202,2},
{"bootstrapCI",(DL_FUNC) &bootstrapCI,4},
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

