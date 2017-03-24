//----------------------------------------------------------------
// Name        : exactdistribution.cpp
// Author      : Remco Hoogenboezem
// Version     : 1.0.1
// Copyright   :
// Description : concordance
//----------------------------------------------------------------

#define EPSILON 1.0e-10
#define SEED    437401448   //Some random number that Matlab gave me
//----------------------------------------------------------------

#include "exactdistribution.h"

#include <map>
#include <unordered_map>
#include <algorithm>
#include <string>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------

int         nSubjects;
int *       pMaxBperSubject;
int *       pBperSubject;
int		    iMaxLevel;
int         iKeyLen;

map<long,double>::iterator            iPrevNode;
map<long,double>::iterator            iNode;
map<long,double>::iterator            iLeft;
map<long,double>::iterator            iRight;

pair<map<long,double>::iterator,bool> pResult;


double LogSum(double dLogX,double dLogY)
{
    if(dLogX>=dLogY)
    {
        return dLogX+log(1.0+exp(dLogY-dLogX));
    }

    return dLogY+log(1.0+exp(dLogX-dLogY));
}
//----------------------------------------------------------------
long Concordance(int iPivot)
{
    int     iElement;
    double  dStateElement;
    long    dSum;

    dSum=0;
    for(iElement=0;iElement<nSubjects;iElement++)
    {
		if (iElement!=iPivot){
       		dStateElement=double(pBperSubject[iElement]);
			dSum+=2*(dStateElement*pMaxBperSubject[iElement]-dStateElement*dStateElement);
		}
    }
    return dSum;
}
//----------------------------------------------------------------
double Probability(int iElement)
{
    int iSum;

    iSum=pMaxBperSubject[iElement]-pBperSubject[iElement];

    for(iElement++;pMaxBperSubject[iElement-1]==pMaxBperSubject[iElement] && pBperSubject[iElement-1]==pBperSubject[iElement];iElement++)
    {
        iSum+=pMaxBperSubject[iElement]-pBperSubject[iElement];
    }

    return log(double(iSum));
}
//----------------------------------------------------------------
void State(const string pKeystr)
{
    int iElement;
    for(iElement=0;iElement<nSubjects;iElement+=2)
    {
        pBperSubject[iElement]=(int(pKeystr[iElement>>1])>>4)&15;
        pBperSubject[iElement+1]=int(pKeystr[iElement>>1])&15;
    }
}
//----------------------------------------------------------------

const string Key202(void)
{
    int     iElement;
    char* pKey = (char*)malloc(iKeyLen);

    for(iElement=0;iElement<nSubjects;iElement+=2)
    {
               pKey[iElement>>1]=char(pBperSubject[iElement]<<4)|char(pBperSubject[iElement+1]);
    }
    const std::string s(pKey,iKeyLen);
    free(pKey);

    return s;
}

//----------------------------------------------------------------
void Node(int iElement,map<long,double> * pPrevNode,unordered_map< string, map<long,double> >* pLevel)
{
    //Creates a new node from the prevNode by incrementing element iElement
    double              dProbability;
    long                dConcordance;
    map<long,double> *  pNode;

    dProbability=Probability(iElement);

    pBperSubject[iElement]++;

    dConcordance=Concordance(iElement);


    const string pKeystr=Key202();//Obtain the key for the current pBperSubject state
    if((*pLevel).count(pKeystr)==0) //If first encounter of this key
    {
        pNode=&(*pLevel)[pKeystr];

        for(iPrevNode=pPrevNode->begin();iPrevNode!=pPrevNode->end();iPrevNode++)
        {
            (*pNode)[iPrevNode->first+dConcordance]=iPrevNode->second+dProbability;
        }
    }

    else //If key already seen
    {
        pNode=&(*pLevel)[pKeystr]; //R! free(pKey);

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

                if(fabs((double)(iNode->first-iLeft->first))<EPSILON)
                {
                    iLeft->second=LogSum(iLeft->second,iNode->second);
                    pNode->erase(iNode);
                    continue;
                }
            }

            iRight=pResult.first; iRight++;

            if(iRight!=pNode->end())
            {
                if(fabs((double)(iNode->first-iRight->first))<EPSILON)
                {
                    iRight->second=LogSum(iRight->second,iNode->second);
                    pNode->erase(iNode);
                    continue;
                }
            }
        }
    }


    pBperSubject[iElement]--;

}
//----------------------------------------------------------------
extern "C"{
SEXP exactDistr202(SEXP bn,SEXP skipTests)
{
    int         nElements;
    int         iElement;

    int         iLevel;

    int         b;
    double      dLogSum;
    double      dScale;

    double *    pConcordance;
    double *    pProbability;
    SEXP        plhs;

    unordered_map< string ,map<long,double> >              levels[2];
    unordered_map< string, map<long,double> > *            pLevel;
    unordered_map< string ,map<long,double> > *            pPrevLevel;
    unordered_map< string ,map<long,double> >::iterator    iPrevLevel;

    map<long,double> *    pPrevNode;

    //-------------------------------------------------------------------------
    //Initialize
    //-------------------------------------------------------------------------

    nSubjects= length(bn);
	iKeyLen=(nSubjects>>1)+(nSubjects&1); //=ceiling(n/2)

    pMaxBperSubject=(int*)calloc(nSubjects+1,sizeof(int));
    pBperSubject=(int*)calloc(nSubjects+(nSubjects&1),sizeof(int));

    iMaxLevel=0;
    int imaxB = 0;

	int nSmallerThan1 = 0;
	int nBetween1and15 = 0;
	int nGreaterThan15 = 0;
    for(iElement=0;iElement<nSubjects;iElement++) //from 1 to n
    {
     	b= (INTEGER(bn))[iElement];

		iMaxLevel+=pMaxBperSubject[iElement]=b;
		if (b >imaxB ){imaxB =b; }

		if(b>15) nGreaterThan15++;
       	else if(b>1) nBetween1and15++;
		else if(b<1) nSmallerThan1++;
    }

	if (nGreaterThan15>0 || nBetween1and15==0 || nSmallerThan1>0)
	{
	    error("Number of replicates in all subject cannot exceed 7. In addition there must be at least one subject that has >1 replicates.")  ;
	}

	if (LOGICAL(skipTests)[0]==false){
	    if( (imaxB>6) &  (nSubjects>10))    error("Exact: Number of subjects is limited to 10 in case of 7 replicate measurements")  ;
	    if( (imaxB>5) &  (nSubjects>13))    error("Exact: Number of subjects is limited to 13 in case of 6 replicate measurements")  ;
	    if( (imaxB>4) &  (nSubjects>17))    error("Exact: Number of subjects is limited to 17 in case of 5 replicate measurements")  ;
	    if( (imaxB>3) &  (nSubjects>25))    error("Exact: Number of subjects is limited to 25 in case of 4 replicate measurements")  ;
	    if( (imaxB>2) &  (nSubjects>50))    error("Exact: Number of subjects is limited to 50 in case of 3 replicate measurements")  ;
	    if( (imaxB>1) &  (nSubjects>150))    error("Exact: Number of subjects is limited to 150 in case of 2 replicate measurements")  ;
	}

    sort(pMaxBperSubject,pMaxBperSubject+nSubjects);

    unsigned long omega = 0;
    for(iElement=nSubjects-1;iElement>=0;iElement--)
    {
        omega+=pMaxBperSubject[iElement]*(pMaxBperSubject[iElement]-1)*(iMaxLevel-pMaxBperSubject[iElement]);
    }

    //-------------------------------------------------------------------------
    //Lets go
    //-------------------------------------------------------------------------
    levels[0][Key202()][0]=0.0; //initialize the first node

    pPrevLevel=&levels[0];
    pLevel=&levels[1];

    for(iLevel=1;iLevel<=iMaxLevel;iLevel++)
    {
        R_CheckUserInterrupt();
        for(iPrevLevel=pPrevLevel->begin();iPrevLevel!=pPrevLevel->end();iPrevLevel++)
        {
            State(iPrevLevel->first); //Set the state to the one observed in the given prev node

            if(pBperSubject[0]<pMaxBperSubject[0])
            {
                Node(0,&iPrevLevel->second,pLevel);
            }

            for(iElement=1;iElement<nSubjects;iElement++)
            {
                if(pMaxBperSubject[iElement-1]==pMaxBperSubject[iElement])
                {
                    if(pBperSubject[iElement-1]>pBperSubject[iElement])
                    {
                        Node(iElement,&iPrevLevel->second,pLevel);
                    }
                } else {
                    if(pBperSubject[iElement]<pMaxBperSubject[iElement])
                    {
                        Node(iElement,&iPrevLevel->second,pLevel);
                    }
                }
            }
        }

        pPrevLevel->clear();

        pPrevLevel=&levels[iLevel&1];
        pLevel=&levels[1-(iLevel&1)];
    }

	free(pMaxBperSubject);
    free(pBperSubject);


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

    pPrevLevel->clear();

    for(iElement=0;iElement<nElements;iElement++)
    {
        pProbability[iElement]=exp(pProbability[iElement]-dLogSum);
    }

    UNPROTECT(1);

    return plhs;
    //-------------------------------------------------------------------------
    //Done
    //-------------------------------------------------------------------------
}}
//----------------------------------------------------------------
