
// Programmed by Thorsten Bernholt
// University of Dortmund


#ifndef RMquickH
#define RMquickH

#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "RMquick.h"
#include "RegLine.h"
#include <string>
#include <algorithm>
#include <iostream>

void addMemo(char *t);

// Bei mehrdeutigen Medianen nehme den
//kleineren
//#define MEDIAN_LINKS(__MED) (max(0.0,floor((__MED+1.0)/2.0)-1))
//Anzahl der Elemente links vom Median





int MEDIAN_LINKS(int anz);
//int MEDIAN_RECHTS(int anz);

//größeren
//#define MEDIAN_LINKS(__MED) (ceil((__MED+1.0)/2.0)-1)
#ifndef BORLAND_WIN
int random(int a);                // was soll hier passieren ????
#endif

class Median
{

/*  double *medTab;		// Array der Median werte 
  int medTabSize;       // Länge des Arrays

  	void vertausche(int a,int b)
    {
  		double temp=medTab[a];
      	medTab[a]=medTab[b];
      	medTab[b]=temp;
    };*/
    /* Ist eine Implementierung von 
     * Quickselect
     * die den Median zurückliefert ( wenn n_thElement)
     * geeignet gewählt wird)
     * */
    
   /* int random(int a)
    {
      return (int)((((double)rand())/((double)RAND_MAX))*a);
    };*/

    

/*
    double computeMedian(int left, int right, int tempNthElement){
    	int iPivot = random(right - left) + left;
    	int l = left;
    	int r = right-1;
    	while (l!=r){
    		for (;((l!=iPivot)&&(medTab[l]<=medTab[iPivot]));l++);
    		for (;((r!=iPivot)&&(medTab[r]>=medTab[iPivot]));r--);
    		if (l != r){
    			vertausche(l, r);
    			if (l == iPivot) iPivot = r;
    			else
    				if (r == iPivot) iPivot = l;
    			}
    		};
    		int leftSideSize = iPivot - left;
    		if (leftSideSize == tempNthElement){
    			return medTab[iPivot];
    		} else if (tempNthElement < leftSideSize){
    			*//*Search in the leftSide*//*
    			return computeMedian(left, iPivot, tempNthElement);
    		} else if (tempNthElement > leftSideSize){
    			*//*Search in the rightSide*//*
    			return computeMedian(iPivot+1, right, tempNthElement-(leftSideSize+1));
    		}
    		return -1; // kann 
    	}*/
/*    double computeMedian(int L,int R,int n_thElement)
    {
       if (R-L<=0)
         return medTab[L];

       int a=L,b=R;
       int pivo=random(R-L+1)+L;

       while( a<pivo || b>pivo)
       {
         while(medTab[a]<=medTab[pivo] && a<pivo)
           a++;

         while(medTab[b]>medTab[pivo] && b>pivo)
           b--;

         if (a<pivo && b>pivo)
         {
           vertausche(a,b);
           a++;
           b--;
         }
         else
         {
           if (a==pivo)
           {
              vertausche(a,b);
              a++;
              pivo=b;
           }
           else
           if (b==pivo)
           {
              vertausche(a,b);
              b--;
              pivo=a;
           }
         }
       }

       if (pivo<n_thElement) // median_links = berechnet # hälfte des arrays
       {
         return computeMedian(pivo+1,R,n_thElement);
       }
       else
       if (pivo>n_thElement)
       {
         return computeMedian(L,pivo-1,n_thElement);
       }
       else
         return medTab[pivo];

    };*/
 public:
  double getMedian(double *tab,int anz)
  {
    //XXX! This could be implemented in terms of
    //STL's nth_element:
    //http://www.sgi.com/tech/stl/nth_element.html
     if (anz % 2 == 0)   
     { 
    	//double median_links = computeMedian(0,medTabSize-1,MEDIAN_LINKS(medTabSize));
     	//double median_rechts = computeMedian(0,medTabSize-1,MEDIAN_LINKS(medTabSize)-1);
     	std::nth_element(tab,tab + MEDIAN_LINKS(anz), tab+anz);
     	double median_links = tab[MEDIAN_LINKS(anz)];
     	std::nth_element(tab,tab + MEDIAN_LINKS(anz)-1, tab+anz);
     	double median_rechts = tab[MEDIAN_LINKS(anz)-1];
     	return (median_links+median_rechts)/2.0;
     } else {
    	std::nth_element(tab,tab + MEDIAN_LINKS(anz), tab+anz);
     	return tab[MEDIAN_LINKS(anz)];
     }
  }
};

class RMbase
{
  public:
  virtual ~RMbase() {};
  virtual char * getName(void)=0;
  virtual void addPunkt(double x,double y)=0;
  virtual RegLine getRM(double timeZero)=0;
  virtual RegLine getLMS(double timeZero){return(0);};
  virtual void print(void){};
  virtual void init(int windowSize)=0;
  virtual void lout(int nr){};
};

class RMquick:public RMbase
{
public:
  class Punkt
  {
    public:
    double x,y;
  };
  Punkt *tab;
  int windowSize;
  int index;  //zeigt auf das nächste zu ersetzende Element
  int max;
  Median med;
  double *slope;
  double *smallMed;
  double *smallMed2;
  
  double getSlope(int a,int b)
  {
    //double ya=tab[a].y;
    //double yb=tab[b].y;
    //double xa=tab[a].x;
    //double xb=tab[b].x;
    return (tab[a].y-tab[b].y)/(tab[a].x-tab[b].x);
  };

public:
  void init(int s)
  {
    windowSize=s;
    tab=new Punkt[windowSize];
    index=0;
    max=0;

    slope=new double[windowSize];
    smallMed=new double[windowSize];
    smallMed2=new double[windowSize];
  };
  void print(void)
  {
  };
  char * getName(void)
  {
    return ((char*)("RMquick"));
  };
  void addPunkt(double x,double y)
  {
    tab[index].x=x;
    tab[index].y=y;

    if (max<windowSize)
      max++;
    index++;
    if (index==windowSize)
      index=0;
  };
  RegLine getRM(double timeZero)
  {

    //timeZero nicht implementiert
    RegLine l;
    int s;
    for(int k=0;k<max;k++)
    {
      s=0;
      for(int i=0;i<max;i++)
        if (i!=k)
        {
          slope[s]=getSlope(i,k);
          s++;
        }

      smallMed[k]=med.getMedian(slope,max-1);
    }
    const double steigung = med.getMedian(smallMed,max); //Vorzeichen wird in der Schnitt-Funktion verdreht
    for(int k=0;k<max;k++)
      smallMed2[k] = tab[k].y - steigung * tab[k].x;

    l.steigung = steigung;
    l.y_achse=med.getMedian(smallMed2,max);
    l.fitness=0;
    return(l);
  }
};


#endif
