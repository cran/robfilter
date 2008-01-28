
// Programmed by Thorsten Bernholt
// University of Dortmund


#ifndef RMquickH
#define RMquickH

#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "RMquick.h"
#include "RegLine.h"

void addMemo(char *t);

// Bei mehrdeutigen Medianen nehme den
//kleineren
//#define MEDIAN_LINKS(__MED) (max(0.0,floor((__MED+1.0)/2.0)-1))
//Anzahl der Elemente links vom Median





int MEDIAN_LINKS(int anz);
//int MEDIAN_RECHTS(int anz);

//gr��eren
//#define MEDIAN_LINKS(__MED) (ceil((__MED+1.0)/2.0)-1)
#ifndef BORLAND_WIN
int random(int a);                // was soll hier passieren ????
#endif

class Median
{
  double *medTab;		// Array der Median werte 
  int medTabSize;       // L�nge des Arrays

    void vertausche(int a,int b)
    {
     double temp=medTab[a];
      medTab[a]=medTab[b];
      medTab[b]=temp;
    };
    /* Ist eine Implementierung von 
     * Quickselect
     * die den Median zur�ckliefert ( wenn n_thElement)
     * geeignet gew�hlt wird)
     * */
    double computeMedian(int L,int R,int n_thElement)
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

       if (pivo<n_thElement) // median_links = berechnet # h�lfte des arrays
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

    };
 public:
  double getMedian(double *tab,int anz)
  {
    //XXX! This could be implemented in terms of
    //STL's nth_element:
    //http://www.sgi.com/tech/stl/nth_element.html
         
     medTab=tab;
     medTabSize=anz;
     /* falls die Anzahl der X Werte der Schnittpunkte
      * gerade sind wird der Median als Durchschnitt
      * der beiden Werte in der mitte gebildet
      */
     if (medTabSize % 2 == 0)   
     {
     	double median_links = computeMedian(0,medTabSize-1,MEDIAN_LINKS(medTabSize));
     	double median_rechts = computeMedian(0,medTabSize-1,MEDIAN_LINKS(medTabSize)-1);
     	return (median_links+median_rechts)/2.0;
     } else return(computeMedian(0,medTabSize-1,MEDIAN_LINKS(medTabSize)));
  };
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
  int index;  //zeigt auf das n�chste zu ersetzende Element
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
  char *getName(void)
  {
    return((char *)"RMquick");
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
