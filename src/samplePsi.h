#ifndef __SAMPLEPSI_H__
#define __SAMPLEPSI_H__

#define STRICT_R_HEADERS 

extern "C" {
	SEXP samplePsi(SEXP rP_choleski,SEXP rP_mat1Missing,SEXP rP_mat2Missing,SEXP rP_nDraws,SEXP rP_nCPU);
}

#endif /*__SAMPLEPSI_H__*/
