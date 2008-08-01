//---------------------------------------------------------------------------

/**
 * Changes:
 * - 20070801 (OME): 
 *   + Change sqrtl() to sqrt()
 *   + Change powl() to pow()
 */
// Changes:
//#pragma hdrstop

#include "LQDAdvanced.h"

/** Remove the oldest point added that has not
  * yet been removed. Do nothing if there is no
  * such point.
  */
void LQDAdvanced::removePunkt()
{
     if (anzPunkte < 1) {
        return;
     }
     
     X.pop_front();
     Y.pop_front();
     --anzPunkte;
}

/**
 * The method searchOptimalPoint performs a search for an optimal
 * point. It iteratively reduces a search interval containing local
 * solutions.  
 */
  void LQDAdvanced::searchOptimalPoint(int method,double power)
  {
    switch (method)
    {
      case 1:
      {
        double f;
        double divisor=1.0;
        double factor;
        bool lastDecision=true;
        bool endGeometric=false;
        if ((lqdDetermineStart==2) && (lqdDescendMethod>1)) decideLQD(fmax,1); // there would not be a line to descend otherwise
        if (lqdDescendMethod>1) fmax=min(fmax,descendToLocalSolution()+buffer);
//        while(fmax-fmin>eps)
//        while(fmax>(1+eps)*fmin)
        while((1-eps)*fmax>fmin)
        {
//          if (lastDecision) f=(fmin+(fmax-fmin)/divisor); else f=(fmax-(fmax-fmin)/divisor);
// ---------------------------------------------------------------Geometric mean---------------------------------------------------------------
          if ((anzPunkte*(1-eps)*fmax)>fmin) endGeometric=true; // just log n steps to the optimal point
          if (endGeometric) f=sqrt(fmin*fmax);
          else
          {
	    // factor=1/powl(2,divisor);
	    // OME: Better?
	    // factor = 1/exp2(divisor);		 
	    factor = 1/pow(2, divisor);

	    // OME: Changed powl() to pow():
            if (!lastDecision) 
	      f = (pow(fmin,factor)*
		   pow(fmax, (1-factor))); 
	    else 
	      f = (pow(fmin,(1-factor)) *
		   pow(fmax,factor));
            if ((f==fmin) || (f==fmax) || (f<fmin) || (f>fmax)) 
	      f = sqrt(fmin*fmax);
          }

/*        if (!lastDecision)
          {
            f=(powl((long double)fmin,factor)*powl((long double)fmax,(1-factor)));
            if ((f==fmin) || (f==fmax) || (f<fmin) || (f>fmax) || (endGeometric)) f=sqrt(fmin*fmax);
          }
          else // if (fmin>0.0)
          {
            f=(powl((long double)fmin,(1-factor))*powl((long double)fmax,factor));
            if ((f==fmin) || (f==fmax) || (f<fmin) || (f>fmax) || (endGeometric)) f=sqrt(fmin*fmax);
          }

          else
          {
            f=(fmin+factor*(fmax-fmin));
            if ((f==fmin) || (f==fmax) || (f<fmin) || (f>fmax)) f=(fmin+(fmax-fmin)/2.0);
          }*/
          // no way to be more precise
          if ((f==fmin) || (f==fmax) || (f<fmin) || (f>fmax)) {
//            printf("Abbruch! mit fmin=%e f=%e fmax=%e factor=%e divisor=%e wegen %d %d %d %d\n",fmin,f,fmax,factor,divisor,(f==fmin),(f==fmax),(f<fmin),(f>fmax));
            fmin=fmax;
          }
          else
          {
            if (decideLQD(f,1)==true)
            {
              fmax=f;
              if (lqdDescendMethod==3) fmax=min(f,descendToLocalSolution()+buffer);
              if ((lqdDescendMethod==4) && ((anzPunkte*(1-eps)*fmax)>fmin)) fmax=min(f,descendToLocalSolution()+buffer);
              if (lastDecision) divisor += power; else divisor = 1.0;
              lastDecision=true;
            }
            else
            {
              fmin=f;
              if (!lastDecision) divisor += power; else divisor = 1.0;
              lastDecision=false;
            }
            if (fmin>fmax) fmin=fmax;
          }
        }
        if ((lqdDescendMethod==2) || (lqdDescendMethod==4)) fmax=min(fmax,descendToLocalSolution()+buffer);
        break;
      }
      case 2:
      {
        double f;
        decideLQD(fmax,2);
        //crossings=1ui64;					// auskommentiert von Sebastian Ruthe
        crossings=(u__int64)1;
        //while(crossings!=0ui64)			// auskommentiert von Sebastian Ruthe
        while(crossings!=(u__int64)0)
        {
          //crossings=0ui64;				// auskommentiert von Sebastian Ruthe
          crossings=(u__int64)0;
          perm.resize(fmaxCuts->cutp.size());
          permutedCuts.resize(fmaxCuts->cutp.size());
          temp.resize(fmaxCuts->cutp.size());
          inversionTable.resize(fmaxCuts->cutp.size());
          sort(fmaxCuts->cutp.begin(),fmaxCuts->cutp.end(),lessFMax(*this));
          const int cutpSize = fmaxCuts->cutp.size();
          for (int i=0; i < cutpSize; i++) {
            perm[fmaxCuts->cutp[i].p->origin]=i;
            inversionTable[i]=0;
          }
          for (int i=0; i<cutpSize; i++) {
            permutedCuts[i]=perm[fminCuts->cutp[i].p->origin];
          }
          mergeSort(0,(fminCuts->cutp.size()-1));
          //if (crossings!=0ui64) {			// auskommentiert von Sebastian Ruthe
          if (crossings!=(u__int64)0) {
            //u__int64 selectCrossing=randomRange64(crossings-1ui64)+1ui64;
            u__int64 selectCrossing=randomRange64(crossings-(u__int64)1)+ (u__int64)1;
            f=computeNthCrossing(selectCrossing);
            if (decideLQD(f,2)==true) fmax=f; else fmin=f;
          }
        }
        break;
      }
    }
    _fmax_output=fmax;
  }

bool LQDAdvanced::lessBW::operator()(const LQDAdvanced::p_cutAndInfo &firstOp,const LQDAdvanced::p_cutAndInfo &secondOp) const
{
	return ((firstOp.p->value < secondOp.p->value) || ((firstOp.p->value==secondOp.p->value) && (((firstOp.p->origin>secondOp.p->origin)&&(firstOp.p->ascending)&&(!secondOp.p->ascending))||((firstOp.p->origin<secondOp.p->origin)&&((firstOp.p->ascending)||(!secondOp.p->ascending))))));
};



bool LQDAdvanced::lessFMin::operator()(const LQDAdvanced::p_cutAndInfo &firstOp, const LQDAdvanced::p_cutAndInfo &secondOp) const
{
  bool less=true;
  double cutY;
  double slopeDiff = computeLQDAdvanced.transformedInput[firstOp.p->origin].slope
                     - computeLQDAdvanced.transformedInput[secondOp.p->origin].slope;
  double interceptDiff = computeLQDAdvanced.transformedInput[firstOp.p->origin].intercept
                      - computeLQDAdvanced.transformedInput[secondOp.p->origin].intercept;
  if (slopeDiff!=0)
  {
    cutY=computeLQDAdvanced.computeYCut(*firstOp.p,*secondOp.p);
    if (computeLQDAdvanced.transformedInput[firstOp.p->origin].slope*computeLQDAdvanced.transformedInput[secondOp.p->origin].slope>0) less=!less;
    if (cutY<=computeLQDAdvanced.fmin) less=!less;
    if (slopeDiff<0) less=!less;
  }
  else
  {
    if (interceptDiff == 0) return firstOp.p->origin < secondOp.p->origin;
    if (interceptDiff < 0) less = !less;
    if (computeLQDAdvanced.transformedInput[firstOp.p->origin].slope < 0) less = !less;
  }
  return less;
}

bool LQDAdvanced::lessFMid::operator()(const LQDAdvanced::p_cutAndInfo &firstOp,const LQDAdvanced::p_cutAndInfo &secondOp) const
{
	bool less=true;
  	double cutY;
  	double slopeDiff=computeLQDAdvanced.transformedInput[firstOp.p->origin].slope-computeLQDAdvanced.transformedInput[secondOp.p->origin].slope;
  	double interceptDiff=computeLQDAdvanced.transformedInput[firstOp.p->origin].intercept-computeLQDAdvanced.transformedInput[secondOp.p->origin].intercept;
  	if (slopeDiff!=0)
    {
      	cutY=computeLQDAdvanced.computeYCut(*firstOp.p,*secondOp.p);
      	if (computeLQDAdvanced.transformedInput[firstOp.p->origin].slope*computeLQDAdvanced.transformedInput[secondOp.p->origin].slope>0) less=!less;
      	if (cutY<=computeLQDAdvanced.fmid) less=!less;
      	if (slopeDiff<0) less=!less;
  	}
  	else
  	{
    	if (interceptDiff==0) return ((firstOp.p->origin<secondOp.p->origin));
    	if (interceptDiff<0) less=!less;
    	if (computeLQDAdvanced.transformedInput[firstOp.p->origin].slope<0) less=!less;
  	}
  	return less;

}

bool LQDAdvanced::lessFMax::operator()(const LQDAdvanced::p_cutAndInfo &firstOp,const LQDAdvanced::p_cutAndInfo &secondOp) const
  {
  	bool less=true;
  	double cutY;
  	double slopeDiff=computeLQDAdvanced.transformedInput[firstOp.p->origin].slope-computeLQDAdvanced.transformedInput[secondOp.p->origin].slope;
  	double interceptDiff=computeLQDAdvanced.transformedInput[firstOp.p->origin].intercept-computeLQDAdvanced.transformedInput[secondOp.p->origin].intercept;
  	if (slopeDiff!=0)
  	{
    	cutY=computeLQDAdvanced.computeYCut(*firstOp.p,*secondOp.p);
    	if (computeLQDAdvanced.transformedInput[firstOp.p->origin].slope*computeLQDAdvanced.transformedInput[secondOp.p->origin].slope>0) less=!less;
    	if (cutY<computeLQDAdvanced.fmax) less=!less;
    	if (slopeDiff<0) less=!less;
  	}
  	else
  	{
    	if (interceptDiff==0) return (((firstOp.p->origin<secondOp.p->origin)));
    	if (interceptDiff<0) less=!less;
    	if (computeLQDAdvanced.transformedInput[firstOp.p->origin].slope<0) less=!less;
  	}
  	return less;
  };

//---------------------------------------------------------------------------

//#pragma package(smart_init)
