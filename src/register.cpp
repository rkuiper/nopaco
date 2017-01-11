#include "getPsi.h"
#include "samplePsi.h"
#include "exactdistribution202.h"
#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[] = {
{"getPsi202",  (DL_FUNC) &getPsi202, 1},
{"exactDistr202",(DL_FUNC) &exactDistr202,3},
{"samplePsi",(DL_FUNC) &samplePsi,5},
NULL
};


extern "C"{
    void R_init_nopaco(DllInfo *info)
    {
        R_registerRoutines(info,NULL,callMethods,NULL, NULL);
    }
}
