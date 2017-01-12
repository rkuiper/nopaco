//----------------------------------------------------------------
// Name        : main.cpp
// Author      : Remco Hoogenboezem
// Version     : 0.99
// Copyright   :
// Description : concordance
//----------------------------------------------------------------


#include "sparsehash/dense_hash_map"
#include <stdint.h>
#include <iostream>
#include <map>
#include "exactdistribution_base.h"
#include "exactdistribution202.h"

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

//----------------------------------------------------------------

//----------------------------------------------------------------
using namespace std;
using google::dense_hash_map;
//----------------------------------------------------------------

map<long,double>::iterator            iPrevNode;
map<long,double>::iterator            iNode;
map<long,double>::iterator            iLeft;
map<long,double>::iterator            iRight;

pair<map<long,double>::iterator,bool> pResult;



double estimateProgress(unsigned int b, unsigned int n,unsigned int level ){
    //A very raw progress estimator
    double f = 1/(b*n*1.0);    
    double lvl  = (double)level;
    if (level == 0){lvl=R_NegInf;}
    else if (level == b*n){lvl=R_PosInf;}
    else { lvl = (level*f)*(level*f);}


    double mean = 0.5;
    double sd = sqrt(f*f*b*n*n/12.0);
    return( 100*pnorm(lvl, mean,sd, 1, 0));
    
}

//----------------------------------------------------------------
char * Key202(void)
{
    int     iElement;
    char *  pKey;
    
    pKey=(char*)malloc(iKeyLen);
    
    for(iElement=0;iElement<nStateElements;iElement+=2)
    {
        pKey[iElement>>1]=char(pState[iElement]<<4)|char(pState[iElement+1]);
    }
    
    return pKey;
}
//----------------------------------------------------------------
void State(const char * pKey)
{
    int iElement;
    
    for(iElement=0;iElement<nStateElements;iElement+=2)
    {
        pState[iElement]=int(pKey[iElement>>1])>>4;
        pState[iElement+1]=int(pKey[iElement>>1])&15;
    }
}
//----------------------------------------------------------------
void Node(int iElement,map<long,double> * pPrevNode)
{
    double                  dProbability;
    long                  dConcordance;

    char *                  pKey;
    map<long,double> *    pNode;
    
    dProbability=Probability(iElement);
    
    pState[iElement]++;
    
    dConcordance=Concordance(iElement);
   
    
    pKey=Key202();
    
	
    if((*pLevel).count(pKey)==0) //If first encounter of this key
    {
        pNode=&(*pLevel)[pKey];
        
        for(iPrevNode=pPrevNode->begin();iPrevNode!=pPrevNode->end();iPrevNode++)
        {
			

            (*pNode)[iPrevNode->first+dConcordance]=iPrevNode->second+dProbability;
        }
    }
    
    else //If key already seen
    {
        pNode=&(*pLevel)[pKey]; free(pKey);
        
        for(iPrevNode=pPrevNode->begin();iPrevNode!=pPrevNode->end();iPrevNode++)
        {
            pResult=pNode->insert(pair<long,double>(iPrevNode->first+dConcordance,iPrevNode->second+dProbability));
            iNode=pResult.first;
            
            if(pResult.second==false)
            {
                iNode->second=LogSum(iNode->second,iPrevNode->second+dProbability);
                continue;
            }
            
            if(iNode!=pNode->begin())
            {
                iLeft=iNode; iLeft--;

                if(fabs(iNode->first-iLeft->first)<EPSILON)
                {
                    iLeft->second=LogSum(iLeft->second,iNode->second);
                    pNode->erase(iNode);
                    continue;
                }
            }

            iRight=pResult.first; iRight++;

            if(iRight!=pNode->end())
            {
                if(fabs(iNode->first-iRight->first)<EPSILON)
                {
                    iRight->second=LogSum(iRight->second,iNode->second);
                    pNode->erase(iNode);
                    continue;
                }
            }
        }
    }
 

    pState[iElement]--;

}
//----------------------------------------------------------------
extern "C"{
SEXP exactDistr202(SEXP bn,SEXP skipTests, SEXP verbose)
{
    int         nElements;
    int         iElement;

    int         iLevel;
    
    int         dStateLimit;
    double      dLogSum;
    double      dScale;
    
    double *    pConcordance;
    double *    pProbability;
    SEXP          plhs;
    
    int *prhs = INTEGER(bn);
    int nprhs = length(bn);

    dense_hash_map<const char *,map<long,double>,CHashFcn,CEqualKey>              levels[2];
    dense_hash_map<const char *,map<long,double>,CHashFcn,CEqualKey> *            pPrevLevel;        
    dense_hash_map<const char *,map<long,double>,CHashFcn,CEqualKey>::iterator    iPrevLevel;

    map<long,double> *    pPrevNode;
    
   
    
    
    //-------------------------------------------------------------------------
    //Initialize
    //-------------------------------------------------------------------------

    nStateElements= nprhs;    
    iKeyLen=(nStateElements>>1)+(nStateElements&1);
	

    pStateLimits=(int*)calloc(nStateElements+1,sizeof(int));
    pState=(int*)calloc(nStateElements+(nStateElements&1),sizeof(int));

    
    iMaxLevel=0;
    int imaxB = 0;
    
	int nSmallerThan1 = 0;
	int nBetween1and7 = 0;
	int nGreaterThan7 = 0;
    for(iElement=0;iElement<nStateElements;iElement++) //from 1 to n
    {
     	dStateLimit= prhs[iElement]; //statelimit = max(b)
         
		iMaxLevel+=pStateLimits[iElement]=dStateLimit;
		if (dStateLimit >imaxB ){imaxB =dStateLimit; }
		
		if(dStateLimit>7) nGreaterThan7++;
       	else if(dStateLimit>1) nBetween1and7++;
		else if(dStateLimit<1) nSmallerThan1++;
    }


	if (nGreaterThan7>0 || nBetween1and7==0 || nSmallerThan1>0) 
	{
	    error("Number of replicates in all subject cannot exceed 7. In addition there must be at least one subject that has >1 replicates.")  ;         
	} 


	if (LOGICAL(skipTests)[0]==false){
	    if( (imaxB>6) &  (nStateElements>10))    error("Exact: Number of subjects is limited to 10 in case of 7 replicate measurements")  ;         
	    if( (imaxB>5) &  (nStateElements>13))    error("Exact: Number of subjects is limited to 13 in case of 6 replicate measurements")  ;         
	    if( (imaxB>4) &  (nStateElements>17))    error("Exact: Number of subjects is limited to 17 in case of 5 replicate measurements")  ;         
	    if( (imaxB>3) &  (nStateElements>25))    error("Exact: Number of subjects is limited to 25 in case of 4 replicate measurements")  ;         
	    if( (imaxB>2) &  (nStateElements>50))    error("Exact: Number of subjects is limited to 50 in case of 3 replicate measurements")  ;  
	    if( (imaxB>1) &  (nStateElements>150))    error("Exact: Number of subjects is limited to 150 in case of 2 replicate measurements")  ;         
	}    

    sort(pStateLimits,pStateLimits+nStateElements);

    unsigned long omega = 0;
    for(iElement=nStateElements-1;iElement>=0;iElement--)
    {
        omega+=pStateLimits[iElement]*(pStateLimits[iElement]-1)*(iMaxLevel-pStateLimits[iElement]);
    }
    
    
    //-------------------------------------------------------------------------
    //Lets go
    //-------------------------------------------------------------------------
    
    levels[0].set_empty_key(NULL);
    levels[1].set_empty_key(NULL);    
    
    levels[0][Key202()][0.0]=0.0;
    
    pPrevLevel=&levels[0];
    pLevel=&levels[1];


    for(iLevel=1;iLevel<=iMaxLevel;iLevel++)
    {
		#if !defined _WIN64 && !defined _WIN32        
        if (LOGICAL(verbose)[0]==true){
            Rprintf("\r Exact distribution  %.1f %           ",estimateProgress(imaxB, nStateElements,iLevel ));        
        }
        #endif
        R_CheckUserInterrupt();
        for(iPrevLevel=pPrevLevel->begin();iPrevLevel!=pPrevLevel->end();iPrevLevel++)
        {
            State(iPrevLevel->first);

            if(pState[0]<pStateLimits[0])
            {
                Node(0,&iPrevLevel->second);
            }
            
            for(iElement=1;iElement<nStateElements;iElement++)
            {
                if(pStateLimits[iElement-1]==pStateLimits[iElement])
                {
                    if(pState[iElement-1]>pState[iElement])
                    {
                        Node(iElement,&iPrevLevel->second);
                    }
                }
              
                else
                {
                    if(pState[iElement]<pStateLimits[iElement])
                    {
                        Node(iElement,&iPrevLevel->second);
                    }
                }
            }
            
            free((void*)iPrevLevel->first);
        }
        
        pPrevLevel->clear();
        
        pPrevLevel=&levels[iLevel&1];
        pLevel=&levels[1-(iLevel&1)];
    }
    
    //-------------------------------------------------------------------------
    //Normalize and create output
    //-------------------------------------------------------------------------
    
    pPrevNode=&(pPrevLevel->begin()->second);
    nElements=pPrevNode->size();


    PROTECT(plhs = allocMatrix(REALSXP, nElements, 2));
    pConcordance=REAL(plhs);
    pProbability=REAL(plhs)+nElements;

    
    dLogSum=log(0.0);
    
    dScale=1.0/omega;
	 
    for(iElement=0,iPrevNode=pPrevNode->begin();iElement<nElements;iElement++,iPrevNode++)
    {
	    pConcordance[iElement]=1- iPrevNode->first*dScale;
        dLogSum=LogSum(dLogSum,(pProbability[iElement]=iPrevNode->second));
    }

    free((void*)pPrevLevel->begin()->first);
    pPrevLevel->clear();
    
    for(iElement=0;iElement<nElements;iElement++)
    {
        pProbability[iElement]=exp(pProbability[iElement]-dLogSum);
    }
    
    UNPROTECT(1);
    #if !defined _WIN64 || !defined _WIN32        
    if (LOGICAL(verbose)[0]==true){    
        Rprintf("\r                                                 ");
        Rprintf("\r");
        R_FlushConsole ();
    }
    #endif
    return plhs;
    //-------------------------------------------------------------------------
    //Done
    //-------------------------------------------------------------------------
}}
//----------------------------------------------------------------
