//----------------------------------------------------------------
// Name        : bootstrapCI.cpp
// Author      : Rowan Kuiper
// Version     : 1.0.1
// Copyright   :
// Description : concordance
//----------------------------------------------------------------

/*
Bootstrap input matrices to find a condidence interval for psi
*/

#if defined _WIN64 || defined _WIN32
	#include <windows.h>
#else
	#include <pthread.h>
	#include <unistd.h>
	#define Sleep(x) usleep((x)*1000)
#endif
#include <deque>
#include <vector>
#include <limits>
#include <iostream>
#include <algorithm>    

#include "getPsi.h"
#include "bootstrapCI.h"

//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------

void * ThreadFunc_bootstrapCI(void * pUserData);
//----------------------------------------------------------------

//----------------------------------------------------------------
//Job class
//----------------------------------------------------------------

class CJob
{
    private:

    public:
		//Whatever you need inside the worker thread for this job
		unsigned int id;
		unsigned int nDraws;
		unsigned int offset;

		CJob(){}

};

//----------------------------------------------------------------
//Jobs class
//----------------------------------------------------------------

class CJobs
{
    private:
    public:
		double* pMAT1;
		double* pMAT2;
		
		int maxB1;
		int maxB2;
		int n1;
		int n2;
		int maxn;
		double* pd_result1;
		double* pd_result2;
		int counter;



		#if defined _WIN64 || defined _WIN32
			CRITICAL_SECTION criticalSection; //for windows
		#else
	        pthread_mutex_t mutex;
		#endif

        deque<CJob>     queue;

		//Whatever you need inside the worker threads for all jobs


		CJobs (double*& pMAT1, double*& pMAT2, int& maxB1,int& maxB2, int& n1,int& n2,  double* pd_result1, double* pd_result2)
		{
			#if defined _WIN64 || defined _WIN32
				InitializeCriticalSection(&criticalSection); //for windows
			#else
	            pthread_mutex_init(&mutex,NULL);
			#endif

			this->pMAT1 = pMAT1;
			this->pMAT2 = pMAT2;
			this->maxB1 = maxB1;
			this->maxB2 = maxB2;
			this->n1 = n1;
			this->n2 = n2;
			this->pd_result1 = pd_result1;
			this->pd_result2 = pd_result2;
			this->counter =0;
			this->maxn = max(n1,n2);
        }

        ~CJobs(void)
        {
			#if defined _WIN64 || defined _WIN32
				DeleteCriticalSection(&criticalSection); //for windows
			#else
				pthread_mutex_destroy(&mutex);
			#endif
        }
};

//----------------------------------------------------------------
//Worker thread
//----------------------------------------------------------------

//For windows
#if defined _WIN64 || defined _WIN32
long unsigned int WINAPI ThreadFuncWin_bootstrapCI(void * pUserData)
{
	ThreadFunc_bootstrapCI(pUserData);
	return 0;
}
#endif
//----------------------------------------------------------------

void * ThreadFunc_bootstrapCI(void * pUserData)
{
	CJob  job;
    CJobs * pJobs;

	pJobs=(CJobs*)pUserData;

	int *indexSelections = (int*)malloc(sizeof(int) *pJobs->maxn);
	
    double *dataMatrix1 = (double*)malloc(sizeof(double) * (pJobs->maxB1)*pJobs->n1);
	double *dataMatrix2 = (double*)malloc(sizeof(double) * (pJobs->maxB2)*pJobs->n2);


	for(;;)
    {
        //----------------------------------------------------------------
        //Get job from queue
        //----------------------------------------------------------------

		#if defined _WIN64 || defined _WIN32
			EnterCriticalSection(&pJobs->criticalSection);    //for windows
		#else
		    pthread_mutex_lock(&pJobs->mutex);
		#endif
        if(pJobs->queue.size()==0)
        {

			#if defined _WIN64 || defined _WIN32
				LeaveCriticalSection(&pJobs->criticalSection);    //for windows
			#else
	            pthread_mutex_unlock(&pJobs->mutex);
			#endif
			free(indexSelections);
			free(dataMatrix1);
			free(dataMatrix2); 

			return (void*)pJobs;
        }

        job=pJobs->queue.front();
        pJobs->queue.pop_front();

		#if defined _WIN64 || defined _WIN32
			LeaveCriticalSection(&pJobs->criticalSection);    //for windows
		#else
	        pthread_mutex_unlock(&pJobs->mutex);
		#endif


        //----------------------------------------------------------------
		//Do work on job
        //----------------------------------------------------------------
		unsigned int cycle;
		for (cycle = 0; cycle< job.nDraws; cycle++){
			//if (bSigintReceived) { break;}

			#if defined _WIN64 || defined _WIN32
			EnterCriticalSection(&pJobs->criticalSection);    //for windows
			#else
				pthread_mutex_lock(&pJobs->mutex);
			#endif
            //Create the indices to use for bootstrap			
            GetRNGstate();
            for (int i = 0; i < pJobs->n1; i++){
		        indexSelections[i] = R_unif_index(pJobs->n1);
	        }   
	        PutRNGstate();
			#if defined _WIN64 || defined _WIN32
				LeaveCriticalSection(&pJobs->criticalSection);    //for windows
			#else
			    pthread_mutex_unlock(&pJobs->mutex);
			#endif

            //Fill the bootrapped matrices
            
            for (int i = 0; i < pJobs->n1; i++){
            for (int j = 0; j < pJobs->maxB1; j++){

                dataMatrix1[pJobs->n1*j+i] = pJobs->pMAT1[pJobs->n1*j+indexSelections[i]];
            }}

            for (int i = 0; i < pJobs->n2; i++){
            for (int j = 0; j < pJobs->maxB2; j++){

                dataMatrix2[pJobs->n2*j+i] = pJobs->pMAT2[pJobs->n2*j+indexSelections[i]];
            }}


            pJobs->pd_result1[job.offset +cycle] = getPsi(dataMatrix1, pJobs->n1, pJobs->maxB1);
            pJobs->pd_result2[job.offset +cycle] = getPsi(dataMatrix2, pJobs->n2, pJobs->maxB2);
            
            #if defined _WIN64 || defined _WIN32
				EnterCriticalSection(&pJobs->criticalSection);    //for windows
			#else
			    pthread_mutex_lock(&pJobs->mutex);
			#endif
				pJobs->counter ++;
			#if defined _WIN64 || defined _WIN32
				LeaveCriticalSection(&pJobs->criticalSection);    //for windows
			#else
			    pthread_mutex_unlock(&pJobs->mutex);
			#endif
            

 
		}

	}
}
//----------------------------------------------------------------


void startMultithreadedSampling(double* pMAT1,double* pMAT2,unsigned int nDraws,int maxB1,int maxB2, int n1, int n2,unsigned int nCPU,double* pd_result1, double* pd_result2) {

	unsigned long i;
	unsigned int iThread;

	CJobs jobs(pMAT1, pMAT2 , maxB1,maxB2, n1,n2, pd_result1, pd_result2);
	#if defined _WIN64 || defined _WIN32
		HANDLE *	pThreads;
		pThreads=(HANDLE*)malloc(nCPU*sizeof(HANDLE));
	#else
		pthread_t *     pThreads;
		pThreads=(pthread_t*)malloc(nCPU*sizeof(pthread_t));

	#endif

	//----------------------------------------------------------------
	//Prepare jobs
	//----------------------------------------------------------------
	CJob* job = (CJob*)malloc(nCPU*sizeof(CJob));
	unsigned int offset = 0;
	for (i = 0; i < nCPU; i++){
		job[i].id=i;
		job[i].nDraws = (nDraws-offset)/(nCPU-i);
		job[i].offset = offset;
		offset +=  	job[i].nDraws;
		jobs.queue.push_back(job[i]);
	}
	//----------------------------------------------------------------
	//Create threads
	//----------------------------------------------------------------


	for(iThread=0;iThread<nCPU;iThread++)
	{
		#if defined _WIN64 || defined _WIN32
			pThreads[iThread]=CreateThread(NULL,0, ThreadFuncWin_bootstrapCI,(void*)&jobs,0, NULL); //for windows
		#else
			pthread_create(&pThreads[iThread],NULL, ThreadFunc_bootstrapCI, (void*)&jobs);
		#endif
	}
	//----------------------------------------------------------------
	//Wait for threads to finish
	//----------------------------------------------------------------
	for(iThread=0;iThread<(nCPU);iThread++)
	{
		#if defined _WIN64 || defined _WIN32
		while(WaitForMultipleObjects(nCPU,pThreads,TRUE,100)){
		} //for windows
		#else
		while (pthread_join(pThreads[iThread],0)) {
			Sleep(100);
		}
		#endif

	}

	free(job);
	free(pThreads);
}
//----------------------------------------------------------------

 

extern "C" {
    SEXP bootstrapCI(SEXP MAT1,SEXP MAT2,SEXP rP_nDraws,SEXP rP_nCPU){
		//MAT1 matrix1
		//MAT2 matrix2 or NULL
		//rP_bn: A vector of numbers of observations per subject
		//rP_nDraws: An integer for the number of matrices to draw

        int nrow1 = 0,nrow2 = 0,ncol1 = 0,ncol2 = 0;
   
        SEXP dim1, dim2;

		rP_nDraws = PROTECT(Rf_coerceVector(rP_nDraws, INTSXP));
		rP_nCPU = PROTECT(Rf_coerceVector(rP_nCPU, INTSXP));
    	MAT1 = PROTECT(Rf_coerceVector(MAT1, REALSXP));

		/*check input parameters*/
        PROTECT(dim1 = Rf_getAttrib( MAT1, R_DimSymbol ) );
		nrow1 = INTEGER(dim1)[0];
		ncol1 = INTEGER(dim1)[1];
  
        if (MAT2!=R_NilValue){
    		MAT2 = PROTECT(Rf_coerceVector(MAT2, REALSXP));
            dim2 = PROTECT(Rf_getAttrib( MAT2, R_DimSymbol ));
		    nrow2 = INTEGER(dim2)[0];
		    ncol2 = INTEGER(dim2)[1];

            if ((nrow1!=nrow2) | (ncol1!=ncol2))  {Rf_error("Dimensions of both matrices must be the same.");}
        }
        
        
		int nCPU = *INTEGER(rP_nCPU);
		if (nCPU > 64 ) {Rf_error("nCPU must be < 65.");}
        SEXP output1;
		PROTECT(output1 = Rf_allocMatrix(REALSXP, *INTEGER(rP_nDraws),2));
		
        if (nrow2>0){            
			startMultithreadedSampling(REAL(MAT1),REAL(MAT2),*INTEGER(rP_nDraws),ncol1,ncol2, nrow1, nrow2, nCPU, REAL(output1),REAL(output1)+*INTEGER(rP_nDraws));
		} else {
            startMultithreadedSampling(REAL(MAT1),NULL,*INTEGER(rP_nDraws),ncol1,ncol2,  nrow1, nrow2, nCPU, REAL(output1),REAL(output1)+*INTEGER(rP_nDraws));
		}
        UNPROTECT(5);
        if (nrow2>0){ UNPROTECT(2); }
		return output1;
        
    
	}

}
//----------------------------------------------------------------
