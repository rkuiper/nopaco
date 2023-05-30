//----------------------------------------------------------------
// Name        : getPsi.h
// Author      : Rowan Kuiper
// Version     : 1.0.1
// Copyright   :
// Description : concordance
//----------------------------------------------------------------

#ifndef __GETPSI_H__
#define __GETPSI_H__

#define STRICT_R_HEADERS
//----------------------------------------------------------------
#include <R.h>
#include <Rinternals.h>
//----------------------------------------------------------------



double getPsi(double* MAT1, unsigned int n, unsigned int maxB);
//----------------------------------------------------------------

extern "C" {
SEXP getPsi202(SEXP MAT1);
}
//----------------------------------------------------------------


class DataClass
{
    protected:
		unsigned int seed;
		unsigned int nrow, ncol;
		double * sMAT;
		double * rMAT;
		double * qMAT;

		unsigned int* BN;
		unsigned long T;
		unsigned long omega;

		void R2Q( void );
		void S2R( void );
		void orderPerSubject( void );
		void BN_from_S( void );



	 public:
		~DataClass();

		double calculatePSI(void );

		void preprocess( void );

		DataClass(double* pmat1, unsigned int n, unsigned int maxB);
		DataClass(const DataClass &obj);


};
//----------------------------------------------------------------

#endif /*__GETPSI_H__*/
