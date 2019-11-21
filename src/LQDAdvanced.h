//---------------------------------------------------------------------------

#ifndef LQDAdvancedH
#define LQDAdvancedH

#include "RMquick.h"
#include <vector>
#include <queue>
#include <list>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <functional>


//#include "values.h"
//values.h is Windows-only. Using uint64_t from stdint.h in place of unsigned __int64:
#include <stdint.h>
typedef uint64_t u__int64;
//Replaced all occurrences of "unsigned __int64" with "u__int64"

#include <float.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>

#include <R_ext/Print.h> // for Rprintf
#include <time.h> // for time

using namespace std; 

class lessFMin;
class lessFMax;
class lessBW;
class lessFMid;

/**
 * Least Quartile Difference Regression
 */
class LQDAdvanced
{	
   //Previously defined in robust.cpp and cmd.cpp
   int _noOfSols;

   //Previously defined in robust.cpp and cmd.cpp
   double _fmax_output;
	
  struct cutAndInfo
  {
    double value;
    unsigned int ascending:1;
    unsigned int origin:31;
    
   	bool operator<(const cutAndInfo &b) const {
  		return (value < b.value)
                || ((value==b.value) && ( ((origin>b.origin) && !ascending && b.ascending)
                                        ||((origin<b.origin)&&(!ascending || b.ascending))));
  	}
  	
    bool operator>(const cutAndInfo &b) const {
    	return !(*this < b);
    }
  };
    friend class lessFMin;
    friend class lessFMid;
    friend class lessFMax;
    friend class lessBW;
  
  struct p_cutAndInfo
  {
    cutAndInfo *p;
    
    bool operator<(const p_cutAndInfo &b) const {
    	return (*(this->p) < *b.p);
    }
    
    bool operator>(const p_cutAndInfo &b) const {
         //XXX! Does this give reasonable results when equality holds?
    	return (!(*this < b));
    }
  };
    
  class lessFMax
  {
    private:
    const LQDAdvanced& computeLQDAdvanced;        
  	public:
    lessFMax(const LQDAdvanced& parent) : computeLQDAdvanced(parent) {}
  	bool operator()(const LQDAdvanced::p_cutAndInfo &firstOp, const LQDAdvanced::p_cutAndInfo &secondOp) const;
  };

  class lessFMid
  {
    private:
    const LQDAdvanced& computeLQDAdvanced;        
  	public:
    lessFMid(const LQDAdvanced& parent) : computeLQDAdvanced(parent) {}
  	bool operator()(const LQDAdvanced::p_cutAndInfo &firstOp,const LQDAdvanced::p_cutAndInfo &secondOp) const;
  };


  class lessBW
  {
    private:
    const LQDAdvanced& computeLQDAdvanced;        
  	public:
    lessBW(const LQDAdvanced& parent) : computeLQDAdvanced(parent) {}
 	bool operator()(const LQDAdvanced::p_cutAndInfo &firstOp,const LQDAdvanced::p_cutAndInfo &secondOp) const;
  };

	
  class lessFMin
  {
  //XXX! Doesn't seem to be used anywhere.
    private:
    const LQDAdvanced& computeLQDAdvanced;        
  	public:
    lessFMin(const LQDAdvanced& parent) : computeLQDAdvanced(parent) {}
  	bool operator()(const LQDAdvanced::p_cutAndInfo &firstOp,const LQDAdvanced::p_cutAndInfo &secondOp) const;
  };	

  struct vectorCutAndInfo
  {
    vector<cutAndInfo> cuts;
    vector<p_cutAndInfo> cutp;
  };

  struct line
  {
    double slope;
    double intercept;
  };

  vector<line> transformedInput;

  //XXX! list<T>::size() can be O(n)!
  //Thus we also track the number of points
  //in X and Y ourselves.
  
  //Number of points in X and Y.
  int anzPunkte;
  
  //XXX! Why not use a circular buffer (as in the present
  //implementation of MedianFilter) rather than lists?
  //It's never the case that insertions/deletions are made into the
  //middle of the sequences of points, is it?
  list<double> X;
  list<double> Y;

  vectorCutAndInfo *fmaxCuts;
  vectorCutAndInfo *fminCuts;
  vectorCutAndInfo *fmidCuts;
//  vectorCutAndInfo *xCuts;

  vector<unsigned int> perm;
  vector<unsigned int> permutedCuts;
  vector<unsigned int> inversionTable;
  vector<unsigned int> temp;
  u__int64 crossings;
  double fmax,fmin,fmid;
  int changeIndex;
  int opt;
  int maxLevel;
  int wLen;
  int h;
  double eps;
  double buffer;
  double searchExponent;
  
  /**
   * 1 if lqdApprox == true, 2 otherwise.
   */
  int searchMethod;
  
  char lqdDetermineStart;
  char lqdDescendMethod;
  
  /**
   * Number of entries used in transformedInput.
   */
  int transformedInputSize;
  
  int h_over_2;

  double best;

  double computeExact ()
  {
    vector<line>::iterator lineItI,lineItJ;
    lineItI=transformedInput.begin();
    double cutX;
    double cutY;
    double optX;
    double optY=std::numeric_limits<double>::max();
    _noOfSols=1;
    while(lineItI != transformedInput.end())
    {
      lineItJ=lineItI;
      lineItJ++;
      while(lineItJ != transformedInput.end())
      {
        if (((*lineItI).slope-(*lineItJ).slope)!=0) {
          cutX=((*lineItJ).intercept-(*lineItI).intercept)/((*lineItI).slope-(*lineItJ).slope);
          cutY=((*lineItI).slope * (*lineItJ).intercept-(*lineItJ).slope * (*lineItI).intercept)/((*lineItI).slope-(*lineItJ).slope);
          if (cutY<=optY && cutY>=0) {
            vector<double> cutsY;
            cutsY.resize(transformedInputSize);
            vector<line>::iterator lineIterator;
            lineIterator=transformedInput.begin();
            int count=0;
            while (lineIterator!=transformedInput.end())
            {
              cutsY[count]=(*lineIterator).slope * cutX + (*lineIterator).intercept;
              count++;
              lineIterator++;
            }
            vector<double>::iterator nth = cutsY.begin() + min(transformedInputSize / 2 + h_over_2 - 1, transformedInputSize - 1);
            nth_element(cutsY.begin(),nth,cutsY.end());
            if (*nth<=optY) {
              // W?hle linkeste L?sung
              if (*nth==optY)
              {
                _noOfSols++;
                if (cutX>optX) optX=cutX;
              }
              else
              {
                _noOfSols=1;
                optY=*nth;
                optX=cutX;
              }
            }
          }
        }
        lineItJ++;
      }
      lineItI++;
    }
    _fmax_output=optY;
    return -optX;
  }

  /**
   * The method transformInput transforms the Input to modified dual space.
   */
  void transformInput ()
  {
    transformedInputSize = anzPunkte * (anzPunkte - 1);
	transformedInput.resize(transformedInputSize);
    list<double>::iterator Xi,Xj,Yi,Yj;
    Xi=X.begin();
    Yi=Y.begin();
    int count=0;
    while(Xi != X.end())
    {
      Xj=Xi;
      Yj=Yi;

      Xj++;
      Yj++;
	  //Since Xj is incremented before entering the loop, we execute the inner loop body
	  //X.end() - Xi - 1 times, that is anzPunkte - i - 1 times, where i ranges from
	  //0 to anzPunkte - 1, for each execution of the outer loop. Thus the inner loop
	  //body is executed (anzPunkte - 1 + 1) * (anzPunkte - 1) / 2 times and so
	  //count is incremented anzPunkte * (anzPunkte - 1) times.
      while(Xj != X.end())
      {
        transformedInput[count].slope=(*Xi-*Xj);
        transformedInput[count].intercept=(*Yi-*Yj);
        count++;
        transformedInput[count].slope=-transformedInput[count-1].slope;
        transformedInput[count].intercept=-transformedInput[count-1].intercept;
        count++;
        Xj++;Yj++;
      }
      Xi++;Yi++;
    }
  }

  /*int randomInteger(int low_bound , int up_bound){
  	up_bound -= low_bound;
  	double f = 1.0f;
  	while (f == 1.0f){
  		f = (double (rand())) / ((double)(RAND_MAX));
  	}
  	return low_bound + (int)(f * (up_bound+1));
  }*/
  
  int randomInteger(int low_bound, int up_bound) {
      GetRNGstate();
      int f = (int) runif(low_bound, up_bound);
      PutRNGstate();
      return f;
  }

  /**
   * The method determineStartPoint determines a Start Point for the regression.
   * method1: start at height eps and iteratively double the height until a local solution is found
   * method2: determine a starting point by computing the h over 2 + n over 2 level on a vertical line
   */
  void determineStartPoint(int method)
  {
    switch (method)
    {
      case 1:
      {
        fmax=eps;
        while(decideLQD(fmax,1)==false)
        {
          fmax=fmax*2.0;
        }
        fmin=fmax/2.0;
        break;
      }
      case 2:
      {
        int zufall1;
        int zufall2;
        
        // ******* Sebastian Ruthe ********
        //srand( (unsigned int)time(NULL) ); // besser w?re es nat?rlich dies nur einmal zu tun
        
        // ********************************
               
        do
        {
          zufall1 = randomInteger(0, transformedInputSize - 1);
          zufall2 = randomInteger(0, transformedInputSize - 1);
        }
        while(transformedInput[zufall1].slope-transformedInput[zufall2].slope==0);

        double midX = (transformedInput[zufall2].intercept-transformedInput[zufall1].intercept)/(transformedInput[zufall1].slope-transformedInput[zufall2].slope);
        vector<double> cutY;
        cutY.resize(transformedInputSize);
        vector<line>::iterator lineIterator;
        lineIterator=transformedInput.begin();
        int count=0;
        while (lineIterator!=transformedInput.end())
        {
          cutY[count]=(*lineIterator).slope * midX + (*lineIterator).intercept;
          count++;
          ++lineIterator;
        }
        vector<double>::iterator nth = cutY.begin() + min(transformedInputSize / 2 + h_over_2 - 1, transformedInputSize - 1);
        nth_element(cutY.begin(),nth,cutY.end());
        fmax=*nth;
        fmin=DBL_MIN;
        break;
      }
    }
  }

  /**
   * The method searchOptimalPoint performs a search for an optimal point. It iteratively reduces a search interval containing local solutions.
   */
  void searchOptimalPoint(int method,double factor);

  void mergeSort(int left, int right)
  {
    int mid;
    if (right > left)
    {
      mid = (right + left) / 2;
      mergeSort(left, mid);
      mergeSort(mid+1, right);
      merge(left, mid+1, right);
    }
  }

  void merge(int left, int mid, int right)
  {
    int i, left_end, num_elements, tmp_pos;
    left_end = mid - 1;
    tmp_pos = left;
    num_elements = right - left + 1;

    while ((left <= left_end) && (mid <= right))
    {
      if (permutedCuts[left] <= permutedCuts[mid])

      {
        temp[tmp_pos] = permutedCuts[left];
        inversionTable[permutedCuts[left]]=inversionTable[permutedCuts[left]]+mid-left_end-1;
        crossings=(u__int64)(crossings+(u__int64)mid-(u__int64)left_end-(u__int64)1); // darf man das wegmachen ??
        tmp_pos = tmp_pos + 1;
        left = left +1;
      }
      else
      {
        temp[tmp_pos] = permutedCuts[mid];
        inversionTable[permutedCuts[mid]]=inversionTable[permutedCuts[mid]]+left_end-left+1;
        //crossings=(u__int64)(crossings+(u__int64)left_end-(u__int64)left+1ui64); // auskommentiert von Sebastian Ruthe
        
        crossings=(u__int64)(crossings+(u__int64)left_end-(u__int64)left+(u__int64)1); // ui64 entfernt
        tmp_pos = tmp_pos + 1;
        mid = mid + 1;
      }
    }

    while (left <= left_end)
    {
      temp[tmp_pos] = permutedCuts[left];
      inversionTable[permutedCuts[left]]=inversionTable[permutedCuts[left]]+mid-left_end-1;
      //crossings=(u__int64)(crossings+(u__int64)mid-(u__int64)left_end-1ui64);  // auskommentiert von Sebastian Ruthe
      crossings=(u__int64)(crossings+(u__int64)mid-(u__int64)left_end-(u__int64)1);
      left = left + 1;
      tmp_pos = tmp_pos + 1;
    }
    while (mid <= right)
    {
      temp[tmp_pos] = permutedCuts[mid];
      mid = mid + 1;
      tmp_pos = tmp_pos + 1;
    }
    for (i=0; i <= num_elements; i++)
    {
      permutedCuts[right]=temp[right];
      right = right - 1;
    }
  }

  u__int64 randomRange64(u__int64 to)
  {
    u__int64 random;
    u__int64 toTemp;
//    int bits=(int)ceil(log((double)to)/log(2.0));
    do
    {
      toTemp=to;
      random=0;
      int shifted=0;
      //while (toTemp>0xffffui64)						// auskommentiert von Sebastian Ruthe
      while (toTemp>(u__int64)0xffff)
      {
        //random*=0x10000ui64;                        // auskommentiert von Sebastian Ruthe
        random*=(u__int64)0x10000;
//        printf("Vor  Shift %f\n",(float)toTemp);
//        to>>16;
        // toTemp/=0x10000ui64;							// auskommentiert von Sebastian Ruthe 
        toTemp/=(u__int64)0x10000; 
//        printf("Nach Shift %f\n",(float)toTemp);
        shifted+=16;
//        printf("Random vor Addition %f\n",(float)random);
        //random+=(u__int64)RandomRange(0,0xffff);    // auskommentiert von Sebastian Ruthe
        random+=(u__int64)randomInteger(0,0xffff);
//        printf("Random nach Addition %f\n",(float)random);
      }
      //random+=(((u__int64)RandomRange(0,(int)toTemp) )*(u__int64)pow(2,shifted));   // auskommentiert von Sebastian Ruthe
      random+=(((u__int64)randomInteger(0,(int)toTemp) )*(u__int64)pow(2.0,shifted));
      //XXX! Can we safely use the fact that 1 << shifted == 2^shifted here?
      //1 << shifted should be _much_ faster, is exact, and would not necessitate
      //type conversions.
      
//      printf("Random nach letzter Addition %f\n",(float)random);
//      if (random>to) printf("Gr??er!!!");
    }
    while (random>to);
    return random;
  }

  double computeNthCrossing(u__int64 n)
  {
    //u__int64 reachedCrossing=0ui64;			// auskommentiert von Sebastian Ruthe
    u__int64 reachedCrossing=0;
    unsigned int i=0;
    while (reachedCrossing<n)
    {
      reachedCrossing+=(u__int64)inversionTable[i];
      i++;
    }
    i--;
    reachedCrossing-=(u__int64)inversionTable[i];

    //unsigned int searchedCrossing=int((u__int64)n-(u__int64)reachedCrossing-1ui64);
    unsigned int searchedCrossing=int((u__int64)n-(u__int64)reachedCrossing-(u__int64)1);
    vector<double> possibleCrossings;
    double cutY;
    for (unsigned int j=0;j<fmaxCuts->cutp.size();j++)
    {
//      printf("Problem mit %d %d",i,j);
      cutY=computeYCut((*fmaxCuts->cutp[i].p),(*fmaxCuts->cutp[j].p));
      if ((cutY>fmin) && (cutY<fmax)) possibleCrossings.push_back(cutY);
    }
      vector<double>::iterator nth=possibleCrossings.begin()+searchedCrossing;
      nth_element(possibleCrossings.begin(),nth,possibleCrossings.end());
    if (possibleCrossings.size()==0) {
      return fmax;
    }
    else return (*nth);
  }

  double computeYCut(cutAndInfo a,cutAndInfo b) const
  {
    double cutY=std::numeric_limits<double>::max();
    double slopeDiff=transformedInput[a.origin].slope-transformedInput[b.origin].slope;
    if (slopeDiff!=0) cutY=(transformedInput[a.origin].slope * transformedInput[b.origin].intercept - transformedInput[b.origin].slope * transformedInput[a.origin].intercept)/(slopeDiff);
    return cutY;
  }

  double computeYCut(double slope1, double intercept1, double slope2, double intercept2) const
  {
    double cutY=std::numeric_limits<double>::max();
    double slopeDiff=slope1-slope2;
    if (slopeDiff!=0) cutY=(slope1*intercept2 - slope2*intercept1)/(slopeDiff);
//    if (cutY==0.0) printf("Sl:%f In:%f Sl:%f In:%f\n",slope1,intercept1,slope2,intercept2);
    return cutY;
  }

  vector<cutAndInfo> computeCuts(double slope,double intercept,bool xPos)
  {
    vector<cutAndInfo> cuts(transformedInputSize);
    vector<line>::iterator lineIterator = transformedInput.begin();
    int countAbove=0;
    int countBelow=transformedInputSize - 1;
    int index;
    int count=0;
    while (lineIterator!=transformedInput.end())
    {
      if ((slope-(*lineIterator).slope)!=0)
      {
        if ((*lineIterator).slope<slope)
        {
          index=countAbove;
          cuts[index].ascending=false;
          countAbove++;
        }
        else
        {
          index=countBelow;
          cuts[index].ascending=true;
          countBelow--;
        }
        if (xPos) cuts[index].value=((*lineIterator).intercept - intercept)/(slope-(*lineIterator).slope); else cuts[index].value=computeYCut(slope,intercept,(*lineIterator).slope,(*lineIterator).intercept);
      }
      // Else for parallel lines
      else
      {
        if ((intercept-(*lineIterator).intercept)>0)
        {
          index=countAbove;
          cuts[index].ascending=true;
          countAbove++;
        }
        else
        {
          index=countBelow;
          cuts[index].ascending=false;
          countBelow--;
        }
        if (xPos||(intercept==(*lineIterator).intercept)) cuts[index].value=-std::numeric_limits<double>::max(); else cuts[index].value=std::numeric_limits<double>::max();
      }
      cuts[index].origin=count;
      count++;
      lineIterator++;
    }
    changeIndex = countBelow; // last descending
    return cuts;
  }

  /**
   * The method descendToLocalSolution searches the lowest solution under fmax on a line defined by a and b.
   */
  double descendToLocalSolution()
  {
    vector<cutAndInfo> cuts;
    vector<p_cutAndInfo> cutp;
    int actLevel=0;
    int sameLines=0;
    double localSolution=fmax;
    cuts = computeCuts(transformedInput[opt].slope,transformedInput[opt].intercept,false);
    cutp.resize(cuts.size());
    const int cutpSize = cutp.size();
    for (int i = 0; i < cutpSize; ++i)
    {
      cutp[i].p=&cuts[i];
    }
    p_cutAndInfo p_compareTo;
    cutAndInfo compareTo;
    p_compareTo.p=&compareTo;
    compareTo.value=fmax;
    compareTo.ascending=true;
    compareTo.origin=std::numeric_limits<short>::max();
    vector<p_cutAndInfo>::iterator sortEnd=partition(cutp.begin(),cutp.end(),std::bind(less<p_cutAndInfo>(),std::placeholders::_1,p_compareTo));
    vector<p_cutAndInfo>::iterator sameLine=cutp.begin();
    while (sameLine!=sortEnd)
    {
      if (((*sameLine).p->value==-std::numeric_limits<double>::max()) || (((*sameLine).p->value==0.0) && (transformedInput[(*sameLine).p->origin].slope<=0.0) && (transformedInput[(*sameLine).p->origin].slope>transformedInput[opt].slope))) sameLines++;
      sameLine++;
    }
    actLevel+=sameLines;
    compareTo.value=0;
    compareTo.ascending=true;
    compareTo.origin=std::numeric_limits<short>::max();
    sortEnd=partition(cutp.begin(),sortEnd,std::bind(greater<p_cutAndInfo>(),std::placeholders::_1,p_compareTo));
    sort(cutp.begin(),sortEnd,lessBW(*this));
    vector<p_cutAndInfo>::iterator sortBegin=cutp.begin();
    while(sortBegin!=sortEnd)
    {
      if ((*sortBegin).p->ascending) actLevel++;
      if (actLevel>=h_over_2)
      {
        localSolution=(*(sortBegin)).p->value;
        best=-(localSolution-transformedInput[opt].intercept)/transformedInput[opt].slope;
        sortBegin=sortEnd-1;
      }
      if (!(*sortBegin).p->ascending) actLevel--;
      sortBegin++;
    }
    return localSolution;
  }

  /**
   * The method decideLQD decides for a given height if there is a local optimum.
   */
  bool decideLQD(double f,int method)
  {
    eps_steps++;
    double actBest=0;
    int countSolutions=0;
    switch (method)
    {
      case 1:
      {
        fmidCuts->cuts = computeCuts(0.0, f, true);
        break;
      }
      case 2:
      {
        fmidCuts->cuts.resize(transformedInputSize);
        for (int k=0; k < transformedInputSize; k++)
        {
          fmidCuts->cuts[k].origin=k;
          fmidCuts->cuts[k].ascending=true;
          if ((transformedInput[k].slope<0) || ((transformedInput[k].slope==0.0) && (transformedInput[k].intercept<f))) fmidCuts->cuts[k].ascending=false;
        }
        break;
      }
    }
    fmidCuts->cutp.resize(transformedInputSize);

    for (int i = 0; i < transformedInputSize; ++i)
    {
      fmidCuts->cutp[i].p=&fmidCuts->cuts[i];
    }

    switch (method)
    {
      case 1:
      {
	    if (h_over_2 - 1 > static_cast<int>(fmidCuts->cutp.size()) - h_over_2) {
			throw std::logic_error("ho2thLCut > ho2thRCut");
		}
	  
        vector<p_cutAndInfo>::iterator ho2thLCut = fmidCuts->cutp.begin()+h_over_2-1;
        vector<p_cutAndInfo>::iterator ho2thRCut = fmidCuts->cutp.end()-h_over_2;
        nth_element(fmidCuts->cutp.begin(),ho2thLCut,fmidCuts->cutp.end());
        nth_element(ho2thLCut+1,ho2thRCut,fmidCuts->cutp.end());
        sort(ho2thLCut+1,ho2thRCut); // der zweite Zeiger zeigt HINTER das letzte zu sortierende Element
        break;
      }
      case 2:
      {
        fmid=f;
        sort(fmidCuts->cutp.begin(),fmidCuts->cutp.end(),lessFMid(*this));
        break;
      }
    }

    vector<p_cutAndInfo>::iterator cutIt=fmidCuts->cutp.begin();
    int level = transformedInputSize / 2 + 1;
    maxLevel = level;
    /*double lastcut; Note SA: commented out because of set but not used warning*/
    while (cutIt!=fmidCuts->cutp.end())
    {
      if (!(*cutIt).p->ascending) {
        level++;
      }
      if (level > maxLevel)
      {
        opt = (*(cutIt)).p->origin;
        maxLevel=level;
      }
      if (level>=(transformedInputSize / 2 + h_over_2 + 1))
      {
        countSolutions++;
        switch (method)
        {
          case 1:
          {
            actBest+=(*(cutIt)).p->value;
            break;
          }
          case 2:
          {
            actBest+=(transformedInput[(*(cutIt)).p->origin].intercept - f)/(-transformedInput[(*(cutIt)).p->origin].slope);
            break;
          }
        }
      }
      if ((*(cutIt)).p->ascending) {
        level--;
      }
      /*lastcut=(*cutIt).p->value;*/
      ++cutIt;
    }
    if (countSolutions!=0) best=-actBest/countSolutions;
    if (maxLevel>=(transformedInputSize / 2 + h_over_2 + 1))
    {
      swap(fmidCuts,fmaxCuts);
      return true;
    }
    else
    {
      swap(fmidCuts,fminCuts);
      return false;
    }
  }


  /**
   * The method computeSlope computes the slope of the LQD regression line.
   */
  double computeSlope(void)
  {
    eps_steps=0;
    transformInput();
    if (decideLQD(0.0, searchMethod))
    {
	  //This branch is used in the default code
      fmax=0.0;
      fmin=0.0;
    }
    else if (decideLQD(DBL_MIN,searchMethod))
    {
      fmax=DBL_MIN;
      fmin=DBL_MIN;
    }
    else
    {
      determineStartPoint(lqdDetermineStart);
      searchOptimalPoint(searchMethod,searchExponent);
    }
    _fmax_output=fmax;
    return best;
  };

public:
       
  LQDAdvanced(void) : anzPunkte(0),
  //XXX! Do the following three objects need to be dynamically allocated?
   fmaxCuts(new vectorCutAndInfo),
   fminCuts(new vectorCutAndInfo),
   fmidCuts(new vectorCutAndInfo),
   crossings(0),
   fmax(0.0), fmin(0.0), fmid(0.0),
   changeIndex(0),
   opt(0),
   maxLevel(0),
   wLen(0),
   h(0),
   eps(0.0),
   buffer(0.0),
   searchExponent(0.0),
   searchMethod(0),
   lqdDetermineStart(0),
   lqdDescendMethod(0),
   transformedInputSize(0),
   h_over_2(0),
   best(0.0)
  {
  }

  //XXX! Used nowhere
  int eps_steps;

  void init(int windowLength,int _h, double _eps, double _searchExponent,char _lqdDetermineStart,bool lqdApprox,char _lqdDescendMethod,double _buffer)
  {       
    wLen=windowLength;
    h=_h;
    eps=_eps;
    searchExponent=_searchExponent;
    lqdDetermineStart=_lqdDetermineStart;
    lqdDescendMethod=_lqdDescendMethod;
    if (lqdApprox) searchMethod=1; else searchMethod=2;
//    printf("Initialisiert mit eps=%e exp=%e startMethode=%d Abstiegsmethode=%d Suchmethode=%d\n",eps,searchExponent,lqdDetermineStart,lqdDescendMethod,searchMethod);
    buffer=_buffer;
	transformedInput.resize(wLen*(wLen-1));
    h_over_2=(h*(h-1))/2;
    anzPunkte = 0;
  };
  
  
  void adjust_H_wLen(int windowLength,int _h){
	  wLen = windowLength;
	  h = _h;
	  h_over_2=(h*(h-1))/2;
  }
  
  /** Remove the oldest point added that has not
    * yet been removed.
    */
  void removePunkt();
  
  void addPunkt(double x,double y)
  {
    X.push_back(x);
    Y.push_back(y);
    
    ++anzPunkte;
    if (anzPunkte > wLen)
    {
      removePunkt();
    }
  }
  
  int noOfSols() const
  {
      //XXX! This is used nowhere. Perhaps we should delete
      //everything related to it?
      return _noOfSols;
  }
  
  double fmax_output() const
  {
         //XXX! This is used nowhere. Perhaps we should delete
         //everything related to it?
         return _fmax_output;
  }
  
  void setPoints(list<double> x, list<double> y)
  {
       //XXX! This method is not used anywhere.
       anzPunkte = x.size();
       //XXX! We should check that x.size() == y.size().
    X=x;
    Y=y;
  }
  
  double getSlope()
  {
    return computeSlope();
  }

  RegLine getLQD(double timeZero)
  {
  	    if (h_over_2 - 1 > anzPunkte * (anzPunkte -1) - h_over_2) {
			return RegLine::nullLine;
		}
  
    RegLine line;

    line.steigung = computeSlope();
    list<double>::iterator Xi,Yi;

    Xi=X.begin();
    Yi=Y.begin();

    vector<double> smallMed(anzPunkte);

    for (int i = 0; i != anzPunkte; ++i, ++Xi, ++Yi) {
      smallMed[i] = *Yi - line.steigung * (*Xi-timeZero);
    }
    vector<double>::iterator nth = smallMed.begin() + anzPunkte / 2;

    //compute median
    //XXX! Should we be using the average of the middle two values
    //when anzPunkte is even?
    nth_element(smallMed.begin(), nth, smallMed.end());
    line.y_achse = *nth;
    return line;
  }




};  

//---------------------------------------------------------------------------
#endif
