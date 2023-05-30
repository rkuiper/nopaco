//----------------------------------------------------------------
// Name        : getPsi.cpp
// Author      : Rowan Kuiper
// Version     : 1.0.1
// Copyright   :
// Description : concordance
//----------------------------------------------------------------

#include <iostream> //DEBUG ONLY!!!
#include <vector>
#include <map>
#include <time.h>
#include <algorithm>
#include <deque>
#include <limits>
#include <stdlib.h>
#include "getPsi.h"

#include <R.h>
#include <Rmath.h>
//----------------------------------------------------------------




//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------


typedef struct
{
    long    lPosition;
    double  dValue;

}lKEY_dVALUE_PAIR;
//----------------------------------------------------------------

typedef struct
{
    long  lPosition;
    long	lValue;

}lKEY_lVALUE_PAIR;
//----------------------------------------------------------------



struct PREDICATE_K_V_PAIR
{
   bool operator() (lKEY_dVALUE_PAIR pair1,lKEY_dVALUE_PAIR pair2)
    {
	     return pair1.dValue<pair2.dValue;
    }
	bool operator() (lKEY_lVALUE_PAIR pair1,lKEY_lVALUE_PAIR pair2)
    {
		return pair1.lValue<pair2.lValue;
    }
} predicate_k_v_pair;
//----------------------------------------------------------------



template <typename T>
inline T  transpose (T mat, int const& nrow,int const& ncol)
{
	int n = nrow*ncol;
	int idx1=0,idx2 = 0;
	int row=0,col=0;
	T placeholder = static_cast<T>(malloc(n * sizeof(mat[0]))) ;
	for (idx1 = 0; idx1 < n;idx1++){
		row = (idx1)%nrow+1;
		col = (idx1)/nrow+1;
		idx2 = (row-1)*ncol+col-1;
		placeholder[idx2] = mat[idx1];
	}
	for (idx1 = 0; idx1 < n;idx1++){
		mat[idx1] = placeholder[idx1];
	}
	free(placeholder);
	return mat;
}
//----------------------------------------------------------------




DataClass::~DataClass(){
	if (this->BN!=NULL){ free(this->BN);}

	if (this->sMAT!=NULL){ free(this->sMAT);}
	if (this->qMAT!=NULL){ free(this->qMAT);}
	if (this->rMAT!=NULL){ free(this->rMAT);}

	this->BN=NULL;
	this->sMAT=NULL;
	this->qMAT=NULL;
	this->rMAT=NULL;
}
//----------------------------------------------------------------

double DataClass::calculatePSI(void ){
	double *QMAT = this->qMAT;
	unsigned int* bn = this->BN ;
	unsigned int n = this->ncol;
	unsigned int bmax = this->nrow;

	unsigned long k,j;

	//std::cerr << "nrow = " << n <<"; ncol = "<< bmax<<endl;
	//Calculate psi
	long theta;
	double* pdTempPsi = (double*)malloc(sizeof(double) * n);
	for (j = 0; j < n; j++){pdTempPsi[j] = 0;}

	for (j = 0; j < n; j++){

		//std::cerr<<"\n";
		for (k = 1; k < bn[j]; k++){
			theta = -2 * k*(k-bn[j]) ; // -2*k*(k-b)
			pdTempPsi[j]+= theta * (QMAT[j*bmax+(k-1)]);
			//std::cerr<< QMAT[j*bmax+(k-1)]<<", ";
		}

	}
	double psi = 0.0;


	for (j = 0; j < n; j++){psi += pdTempPsi[j];}
	psi /= this->omega;

	free(pdTempPsi);
	return(1-psi);
}
//----------------------------------------------------------------


void DataClass::R2Q(){
	//Convert R matrix to Q matrix (i.e. Q_1 = R_2-R_1-1):
	//std::cerr << "Q:\n";
	unsigned int i=0,j=0;
	for (j = 0; j < this->ncol; j++){
		for (i = 0; i < (this->BN[j]-1); i++){
			this->qMAT[j*this->nrow+i] = this->rMAT[j*this->nrow+i+1]-this->rMAT[j*this->nrow+i]-1;
			//std::cerr << "\t" << this->qMAT[j*this->nrow+i];
		}
		//std::cerr<< "\n";
	}
}
//----------------------------------------------------------------




void DataClass::S2R(){
	//Convert to int
	unsigned long count = 0, i = 0, j = 0, k = 0, v = 0, nNAN = 0;
	double l =0.0;

  //a vector (for each element in the matrix) of key value pairs of doubles
	lKEY_dVALUE_PAIR* pKVP = (lKEY_dVALUE_PAIR*)malloc(sizeof(lKEY_dVALUE_PAIR) * (this->ncol* this->nrow)); 

	//First sort all values in the whole matrix	to obtain rank values
	//std::cerr<<"S:\n";
	for (j = 0; j < this->ncol; j++){
	//std::cerr<<"\n";
		for (i = 0; i < this->nrow; i++){
			//std::cerr<< *(sMAT+j* this->nrow+i) <<",";
			pKVP[count].lPosition = count; //Set each key value pair
			pKVP[count++].dValue = *(sMAT+j* this->nrow+i);
			if (isnan(*(sMAT+j* this->nrow+i) ) ) {nNAN++;}
		}
	}

	sort(pKVP,pKVP+(this->nrow*this->ncol),predicate_k_v_pair);

	//Next loop  through the sorted values except the NaNs which are at the end
	i = 0;
	unsigned long from = 0, to = 0, from2 = 0, to2 = 0;


	for (i=0; i < (this->nrow*this->ncol - nNAN);i++){
		//And detect stretches of tied values

		from = i;to=i;
		while (i < (this->nrow*this->ncol - nNAN-1)){
		 	if (pKVP[i].dValue!=pKVP[i+1].dValue){break;}
			i++; //While tied value
		}
		to = i;


		//in case of no tied values:
		if ((to-from)==0) {
			//cout << "\nAssessing NON-tied value : " << pKVP[from].dValue << "(From: " << from << "   to:"<<to << ")"<<endl;
			this->rMAT[pKVP[from].lPosition] = 1 + from ;  //Assign the ranks based on the sorted key value pairs
		} else { //otherwise, if case of tied value:

			double mu = (2+to+from)/2.0;

			//cout << "\nAssessing tied value : " << pKVP[from].dValue << "(From: " << from << "   to:"<<to << ")"<<endl;
			//cout <<"\t mu=" << mu << endl;

			//Make sure that the tied values are sorted by subject.

			for (j = from; j <= to; j++){
				pKVP[j].dValue = pKVP[j].lPosition ;
			}
			sort(pKVP+from,pKVP+to+1,predicate_k_v_pair);

//			for (j = from; j <= to; j++){
//			   cout << "Value: " << pKVP[j].dValue <<"\torder: " << j << "\tposition :"<< pKVP[j].lPosition << "\tsubject :"<< (long)pKVP[j].lPosition/this->nrow << "\n";
//			}

			long curSubj;j=from;
			//Loop through them
			while(j<=to){
				//And determine the begin and end of each subject.
				curSubj = (long) (pKVP[j].lPosition / this->nrow);
				from2 = j;
				while ( j<=to ) {
					if (curSubj == ((long) (pKVP[j].lPosition / this->nrow))){ j++;}
					else {break;}
				}
				to2 = j;

				//Based on the number inside and ouside the subject, calculate and set the expected ranking
				k= to2-from2;
				v = 1+to-from-k;
				//cout << "\tFrom2: " << from2 << "   to2:"<<to2 << endl;
				//cout << "\t"<<k << " tie(s) within subject" << curSubj << " and " << v << "ouside subject"<< endl;
				count = 0;
				for (l = 1-(k+1)/2.0; l<=(k -(k+1)/2.0); l++){
					//cout <<"\t\t\t"<<count <<": adding " <<mu+l*(1.0+(v)/(k+1.0))<<" to position "<<pKVP[from2+count].lPosition<<endl;
					this->rMAT[pKVP[from2+count].lPosition] = mu+l*(1.0+(v)/(k+1.0)) ;
					count++;

				}
			}
		}
	}

	for (i = (this->nrow*this->ncol - nNAN); i < (this->nrow*this->ncol ); i++){
		//this->rMAT[pKVP[i].lPosition] = std::numeric_limits<double>::infinity();  //Assign the missing values as inifinies.
        this->rMAT[pKVP[i].lPosition] = NAN;  //Assign the missing values as inifinies.
	}
 /*std::cerr<<"R:\n";
	for (j = 0; j < (this->ncol ); j++){
		std::cerr<<"\ni";
		for (i = 0; i < (this->nrow ); i++){
			std::cerr<<this->rMAT[ j* this->nrow + i] <<", ";
		}
	}*/
	
	free(pKVP);
}
//----------------------------------------------------------------



void DataClass::orderPerSubject(){
	//pmat1 is a pointer to double matrix, In the rows are the observers (microarrays), in the columns the subjects (patients)
	unsigned int j;
	for (j = 0; j < this->ncol; j++){
		sort(this->sMAT+j* this->nrow, this->sMAT+(j+1)* this->nrow); //sort column j
		//std::cerr<<"\ncol " << j <<": ";
		//for (int k = 0; k < this->nrow; k++) {
		// std::cerr<< *(this->sMAT+j*this->nrow+k)<<",";
	  //}
		//std::cerr<<"\n";
	}
}
//----------------------------------------------------------------

void DataClass::BN_from_S(){
	unsigned long i,j;
	this->T = 0;
	this->omega = 0;
	//Determine BN (= vector of length b indicating the number of known subjects with 1..b measurements)  and T (=total number of pivots) and omega (=b*(b-1)*(t-b)
    if (this->nrow>0){
        for (j = 0; j < this->ncol; j++){ //loop subjects
		    i = this->nrow-1;
		    //while ((*(this->sMAT + j * this->nrow + i))==std::numeric_limits<double>::infinity()){ i--;}
		    while (isnan(*(this->sMAT + j * this->nrow + i))){ i--;}
		    this->BN[j] = i+1; //Set the number of real values for patient j
		    (this->T)+=this->BN[j]; //And the number of pivots
	    }
	    //Now that T is known , we can determine omega
	    for (j = 0; j < this->ncol; j++){
		    this->omega += this->BN[j]*(this->BN[j]-1) * (this->T-this->BN[j]);
	    }
    }
}

//----------------------------------------------------------------


void DataClass::preprocess( void ){
	//Order each column
	this->orderPerSubject();
	//Determine BN (= vector of length b indicating the number of known patients with 1..b measurments)  and T (=total number of pivots irrespective of b: t_b +b = T);
	this->BN_from_S();
	
	//Convert to int
	this->S2R();

	//Determine BN (= vector of length b indicating the number of known patients with 1..b measurments)  and T (=total number of pivots irrespective of b: t_b +b = T);
	//this->BN_from_R();

	this->R2Q();

}
//----------------------------------------------------------------


DataClass::DataClass(double* pmat1, unsigned int n, unsigned int maxB){
	//Assume subject in the rows and observations in the columns
	//Transpose the matrix such that observation are in the rows and subjects in the columns
	transpose(pmat1,n,maxB);

	this->nrow = maxB;
	this->ncol = n;

	this->BN = (unsigned int*)malloc(sizeof(unsigned int) * this->ncol);

	this->sMAT = (double*)malloc(sizeof(double) * (this->ncol*this->nrow));
	this->qMAT = (double*)malloc(sizeof(double) * (this->ncol*this->nrow));
	this->rMAT = (double*)malloc(sizeof(double) * (this->ncol*this->nrow));

	for (unsigned int i = 0; i < (this->ncol*this->nrow); i++){		this->sMAT[i] = pmat1[i];}
	this->T =0;
	this->omega = 0;
	this->seed =1;
	this->preprocess( );
}
//----------------------------------------------------------------


DataClass::DataClass (const DataClass &obj) {
	this->nrow = obj.nrow;
	this->ncol = obj.ncol;

	this->T = obj.T;
	this->omega = obj.omega;
	this->seed = obj.seed;

	this->BN = NULL;
	this->qMAT = NULL;
	this->rMAT = NULL;

	if (obj.BN != NULL) {
		this->BN = (unsigned int*)malloc(sizeof(unsigned int) * this->ncol);
		for (unsigned int i = 0; i < this->ncol; i++){
			this->BN[i] = obj.BN[i];
		}
	}

	if (obj.qMAT != NULL) {
		this->qMAT = (double*)malloc(sizeof(double) * (this->ncol*this->nrow));
		for (unsigned int i = 0; i < (this->ncol*this->nrow); i++){
			this->qMAT [i] = obj.qMAT[i];
		}
	}

	if (obj.rMAT != NULL) {
		this->rMAT = (double*)malloc(sizeof(double) * (this->ncol*this->nrow));
		for (unsigned int i = 0; i < (this->ncol*this->nrow); i++){
			this->rMAT [i] = obj.rMAT[i];
		}
	}
	if (obj.sMAT != NULL) {
		this->sMAT = (double*)malloc(sizeof(double) * (this->ncol*this->nrow));
		for (unsigned int i = 0; i < (this->ncol*this->nrow); i++){
			this->sMAT [i] = obj.sMAT[i];
		}
	}
}
//----------------------------------------------------------------



double getPsi(double* MAT1, unsigned int n, unsigned int maxB){
	DataClass dc(MAT1, n, maxB);
	return( dc.calculatePSI());

}
//----------------------------------------------------------------


extern "C"{
	SEXP getPsi202(SEXP MAT1){
		SEXP MAT1_copy = PROTECT(duplicate(MAT1));
		SEXP dim = getAttrib( MAT1_copy, R_DimSymbol ) ;
		int nrow = INTEGER(dim)[0];
		int ncol = INTEGER(dim)[1];

		
		SEXP psi = allocVector(REALSXP,1);
		*REAL(psi) = getPsi(REAL(MAT1_copy),nrow,ncol);
		UNPROTECT(1);
		return(psi);
	}


}
//----------------------------------------------------------------
