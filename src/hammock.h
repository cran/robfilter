// Programmed by Thorsten Bernholt
// University of Dortmund
// Nicht für den Einsatz in medizinischen Geräten geeignet

#ifndef HammockH
#define HammockH

#include <cassert>
#include "CircularArray.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "Edge.h"
#include <iostream>
#include <string>
#include "Line.h"
#include "RMquick.h"

using namespace std;

#include <R_ext/Error.h> // for error
#include <R_ext/Print.h>

#include "segment_tree.h"

//INF is defined in Line.h

// This class allows one to walk around in the hammock

class FaceIter
{
  Edge *akt;
  int direction;
public:
  FaceIter(void)
  {
    akt=0;
    direction=0;
  }

  ~ FaceIter() {}

  Edge * start(Edge * anfang,int dir)
  {
     akt = anfang;
     direction = dir;

     return akt;
  }

  // there a two cells adjacent to the current edge
  // jump from one cell to the other cell
  void invertDirection(void)
  {
    direction=1-direction;
  }

  Edge *current(void)
  {
    return akt;
  }

  // walk clockwise in the current cell
  // compute the next edge
  Edge *next(int *dir)
  {
    if (akt!=NULL)
      akt=akt->getNext(direction,&direction);

    *dir=direction;
    return akt;
  }
};

//Reusing old memory
template <class T> class Recycle
{
  T **tab; // Zwischenspeicher für gelöschtes
  T *heap; // großer Speicherblock
  int poolMax;
  int heapMax;
public:
	int poolCount;   // nächste freie Stelle
	  int heapCount;

  Recycle(void)
  {
    heapCount=0;
    heapMax=0;
    poolCount=0;
    poolMax=0;
    tab=0;
    heap=0;
  }

  ~Recycle(void)
  {
    delete[]tab;
    delete[]heap;
  }

  // Es werden enBlock heapSize viele Elemente reserviert
  // Und im Fall einer Löschung poolSize viele aufgehoben
  void setSpace(int poolSize,int heapSize)
  {
     heapMax=heapSize;
     heapCount=0;
     heap=new T[heapMax];
     if (poolSize>poolMax)
     {
       poolMax=poolSize;
       T**tab2=new T*[poolMax];
       for(int i=0;i<poolCount;i++)
          tab2[i]=tab[i];

       delete[]tab;
       tab=tab2;
     }
  }

  void freeUsedMemory(){
	  for (int i = 0; i < heapCount; i++){
		  heap[i].~T();
	  }
	  delete [] tab;
	  tab=new T*[poolMax];

	  //heap = new T[heapMax];
	  //heap=new T[heapMax];
	  /*for (int i = 0; i < poolCount; i++){
		  delete tab[i];
	  }*/

	  heapCount = 0;
	  poolCount = 0;
  }


  void shred(T* o)
  {
    if (poolCount<poolMax)
    {
      tab[poolCount]=o;
      poolCount++;
    }
    else
      if (o->erasewithdelete==1)
        delete o;
  }

  T* neu(void)
  {
    T *p;
    if (poolCount>0)
    {
      poolCount--;
      p=tab[poolCount];
      p->Init();
      p->erasewithdelete=0;
      return(p);
    }
    else
    {
      if (heapCount<heapMax)
      {
        p=&heap[heapCount];
        heapCount++;
        p->Init();
        p->erasewithdelete=0;
        return (p);
      }
      else
      {
         // Rprintf("recycle allocated something using new\n");
        p=new T;
        p->erasewithdelete=1;
        return p;
      }
    }
  }

};

class Hammock:public RMbase
{
  Edge *Lup; //Points to the left top edge (for delete)
  Edge *L; // Points to the left buttom edge
  Edge *R; // Points to the right top edge
  Line *border_L,*border_R; // Linke und rechte Grenzen
  int anzLines;
  int windowSize;	// größes des Fensters das über die statistischen Daten geschoben wird
  double *medTab;

  CircularArray<Line*>* lineTab;
  FaceIter iter;
  Recycle<Edge> bin;
  Median med;
  int h; //Size of the covered area for lms, lts, ...

  bool needReInit;
  bool initDone;

  void reInitHammock();


public:

	bool debugInfo;

  Hammock(void)
  {
    border_L=0;
    border_R=0;
    medTab=0;
    lineTab=0;
    debugInfo = false;
  }

  virtual ~Hammock(void)
  {
    delete border_L;
    delete border_R;
    delete [] medTab;

    if (lineTab!=0)
      for(int i=0; i < anzLines; i++)
        delete lineTab->get(i);
    delete lineTab;
  }

  Edge *getFirstEdge(void)
  {
    int dummy;
    Edge* leftNext = L->getNext(1, &dummy);
    if (leftNext) {
       return leftNext;
    }
    else {
       return L;
    }
  }

  void init(int _windowSize)
  {
    windowSize=_windowSize;
    //XXX! Using integer math, it's always the case
    //that ceil(windowSize/2) == windowSize/2 == floor(windowSize/2).
    //Is this what is intended? Or should we instead divide by 2.0?
    h=ceil(windowSize/2)+1;
    //h=(int)ceil(windowSize/2)+1;
    medTab = new double[windowSize];
    lineTab = new CircularArray<Line*>(windowSize);
    anzLines = 0;
    bin.setSpace(3*(windowSize+2)+10,(windowSize+2)*(windowSize+2+1)+2);

    L=bin.neu(); //new Edge();
    R=bin.neu(); //new Edge();

    const double x = -INF;

    Lup=L;

    border_L=new Line(0,0,this);

    border_L->setTypLeft();
    L->setLine(border_L);
    L->set_X_pre(x,0);
    L->set_X_next(x,0);

    border_R=new Line(0,0,this);

    border_R->setTypRight();
    R->setLine(border_R);
    R->set_X_pre(x,0);
    R->set_X_next(x,0);
    needReInit = false;
  }


void checkLines(void)
{
  int memIntersectionCount = -1;
  for(int i=0;i<anzLines;i++){
  	int aktIntersectionCount = checkLine(i);
  	if (memIntersectionCount !=-1 && aktIntersectionCount != 0 && memIntersectionCount != aktIntersectionCount){
  		std::cout << "Somethings wrong here\n";
  	}
  	memIntersectionCount = aktIntersectionCount;
  }
}

int checkLine(const int nr)
{
  //XXX! This method seems to have no visible output or effects.
  //Maybe we should revert to an earlier version and modify that
  //or remove it. It is not used anywhere in the project.

  Edge *a;
  int d;
  const Line* line = lineTab->get(nr);
  a = line->startE;

  int k=0;


  while(a!=0)
  {
    if (a->get_X_next(0)>-INF && a->get_X_next(0)<+INF)
    {
      k++;
    }
    a=a->getNext(0,&d);
    if (a!=0)
      a=a->getPre(d,&d);
  }
  return k;
}



/*
void printLines(void)
{
  for(int i=0;i<anzLines;i++)
    printLine(i);
}

void printLine(const int nr)
{

  Edge *a;
  int d;
  const Line* line = lineTab->get(nr);
  a= line->startE;
  char *buf=new char[100];

  //sprintf(buf,"%i: (%i) ",nr,lineTab[nr]->links);
  sprintf(buf,"%i: (%.1f,%.1f)(%i) ",nr,line->m, line->b, line->links);
  addMemo(buf);

  while(a!=0)
  {
    if (a->get_X_next(0)>-INF && a->get_X_next(0)<+INF)
    {
      sprintf(buf,"%.2f",a->get_X_next(0));
      addMemo(buf);
      if (a==line->median_)
        addMemo(" #  ");
      else
        addMemo("    ");
    }
      a=a->getNext(0,&d);
      if (a!=0)
        a=a->getPre(d,&d);

  }
  addMemo("\r\n");

};
*/

  void adjust_H_wLen(int newWindowSize, int subsetSize)
  {
	  windowSize = windowSize;
	  h=subsetSize;
  }

private:

  double calcSchnitt(Line *a, Line *b)
  {
    double w;
    if (a && b)
      w=(a->schnittX(b));
    else
      w=INF;

    return w;
  }

  Edge *dissect_L(Line *neueLinie);
  void dissect_R(Edge *lose,Line *neueLinie);
  void delLine(void);
public:
  char * getName(void)
  {
    return ((char*) ("Hammock"));
  }

  RegLine getRM(double timeZero);
  void updateRepeatedMedian();
  void computeLMS_old(void);
  void computeLXX(void);

  int getAnzLines(void) {  return anzLines; }

  // The slope m must be greater than that of all lines before.
  // The slope m must be greater equal Zero
  int addPunkt(double m,double b);

  /** Remove the oldest point added that has not
    * yet been removed. If there is no such point,
    * do nothing.
    */
  void removePunkt();

  int addLine(Line *neuL);
  Edge *dissectEdge(Edge *lose,Edge *e,int e_dir);
  void print(void);

  /**
   * Least Median of Squares (Rousseeuw 1984)
   */
  class _LMS
  {
    public:
      double min,min_y,min_x;
      int min_anz;

    RegLine getRegLine(double timeZero)
    {
      RegLine l;
      l.steigung=min_x; //Vorzeichen wird in der Schnitt-Funktion verdreht

      l.y_achse=min_y;
      //Auf Level umrechnen
      l.y_achse= l.steigung * timeZero + l.y_achse;

      l.fitness=min;

      return l;
    }

    void keep_min(double x,double y1,double y2,Hammock *H)
    {
         //XXX! H is not used!

//      if (outD!=0) fprintf(outD,"    LMS: %f %f %f %f",fabs(y1-y2),x,y1,y2);
       if (fabs(y1-y2)<min)
       {
         min=fabs(y1-y2);
         min_y=y1+(y2-y1)/2.0;
         min_x=x;
       }
    };
    void init(void)
    {
      min=INF;
    };
  }LMS;

  /**
   * Least Trimmed Squares (Rousseeuw 1983)
   */
  class _LTS
  {
    int anz;
      double min,min_m,min_b;
      double sum_x_squared;
      double sum_x;
      double sum_y_squared;
      double sum_y;
      double sum_xy;
    public:

    RegLine getRegLine(double timeZero)
    {
      RegLine l;
      l.steigung=min_m; //Vorzeichen wird in der Schnitt-Funktion verdreht

      l.y_achse=min_b;
      //Auf Level umrechnen
      l.y_achse= l.steigung * timeZero + l.y_achse;

      l.fitness=min;
      return(l);
    }

    void clear(void)
    {
//      if (outD!=0)  fprintf(outD,"clear\n");

      sum_x_squared=0;
      sum_x=0;
      sum_y_squared=0;
      sum_y=0;
      sum_xy=0;
      anz=0;
    };

    void init(void)
    {
      min=INF;
      clear();
    }

    void ins(double x,double y)
    {
//      if (outD!=0)  fprintf(outD,"ins: (%f,%f)\n",x,y);

      sum_x_squared+=x*x;
      sum_x        +=x;
      sum_y_squared+=y*y;
      sum_y        +=y;
      sum_xy       +=x*y;
      anz++;
    };
    void del(double x,double y)
    {
//      if (outD!=0)   fprintf(outD,"del: (%f,%f)\n",x,y);

      sum_x_squared-=x*x;
      sum_x        -=x;
      sum_y_squared-=y*y;
      sum_y        -=y;
      sum_xy       -=x*y;
      anz--;
    };

    void keep_min(double x,float y1,float y2,Hammock *H)
    {

      double m=0, m1=0,m2=0, b=0;

      if (anz==0)
        return;


      {
        // compute classical Regression Line
        m1=(
          sum_xy*((double)anz)
          -sum_x*sum_y
          );

        m2=(
          sum_x_squared*((double)anz)
         -sum_x*sum_x
          );

        if (m2!=0)
          m=m1/m2;
        else
          return;

        b=sum_y/((double)anz) - m * sum_x /((double)anz);

      }

      // Compute Sum of squared residuals for the Line (m,b)

      double w=
              +anz*b*b
              +m*m*sum_x_squared
              -2.0*m*sum_xy
              +2.0*m*sum_x*b
              +sum_y_squared
              -2.0*b*sum_y
              ;

//      if (outD!=0) fprintf(outD,"    LTS: %f",w);

      if (w<min)
      {
        min=w;
        min_m=m;
        min_b=b;
      };
    };
  }LTS;


/*
  class _LTA
  {
      double min,min_y,min_x;
      double sum_abs;
    public:

    RegLine * getRegLine(double timeZero)
    {
      RegLine *l=new RegLine;
      l->steigung=min_x; //Vorzeichen wird in der Schnitt-Funktion verdreht

      l->y_achse=min_y;
      //Auf Level umrechnen
      l->y_achse= l->steigung * timeZero + l->y_achse;

      l->fitness=min;
      return(l);
    };

    void init(void)
    {
      min=INF;
    };
  }LTA;

  */

  class RegDepth
  {
    double max_m;
    double max_b;
    int count; //number of occurence
    int max;
  public:
    RegLine getRegLine(double timeZero)
    {
      RegLine l;
      l.steigung=max_m/(double)count; //Vorzeichen wird in der Schnitt-Funktion verdreht

      l.y_achse=max_b/(double)count;
      //Auf Level umrechnen
      l.y_achse = l.steigung * timeZero + l.y_achse;

      l.fitness=max;
      return(l);
    }

    void keepMin(double m,double b ,int depth)
    {
      if (depth==max)
      {
        max_m+=m;
        max_b+=b;
        count++;
      }
      else
        if (depth>max)
        {
          max_m=m;
          max_b=b;
          count=1;
          max=depth;
        }
    }

    void clear(void)
    {
      count=0;
      max=0;
      max_m=0;
      max_b=0;
    }

    Line *linePassed;

    Edge *nextEdgeOnLine(Edge *e)
    {
      Edge *f;
      int fDir,dummy;

      f=e->getNext(0,&fDir);

      if (f!=0)
      {
        linePassed=f->getLine();
        f=f->getPre(fDir,&dummy);
      }
      return f;
    }
  }regDepth;

  void computeRegDepth(void);
};

void printHammock(Hammock *m);




             /*
void test2(void)
{
  Hammock m;
  Line*line[4];

  m.setWindowSize(4);

  line[0]=m.addLine(1,2);
  line[1]=m.addLine(2,-1);
  line[2]=m.addLine(3,8);
  line[3]=m.addLine(4,7);

  RM rm(4);

  rm.addPoint(1,2);
  rm.addPoint(2,-1);
  rm.addPoint(3,8);
  rm.addPoint(4,7);

  double med=rm.getSlope();

  //m.delLine(line[3]);
  //m.delLine(line[0]);

 // line[0]=m.addLine(1,2);
  //line[3]=m.addLine(4,7);
  //m.addLine(5,20);
  //m.delLine(line[1]);
  //m.delLine(line[0]);
  m.print();
};
  */

#ifdef BORLAND_WIN

//#include <vcl.h>

class Stoppuhr_base
{
//   struct time t1,t2,t3;
   LARGE_INTEGER t1,t2,t3,freq,diff;
public:
   Stoppuhr_base()
   {
     QueryPerformanceFrequency(&freq);
   };
 LARGE_INTEGER getFreq(void)
  {
    return freq;
  };
  void start(void)
  {
    QueryPerformanceCounter(&t1);
    //gettime(&t1);
  };
  void stop1(void)
  {
    QueryPerformanceCounter(&t2);
    //gettime(&t2);
  };
  void stop2(void)
  {
    QueryPerformanceCounter(&t3);
 //   gettime(&t3);
  };

  LONGLONG getMSec1(void)
  {
     diff.QuadPart=(t2.QuadPart-t1.QuadPart)*((LONGLONG)1000000)/freq.QuadPart;
     return diff.QuadPart;
  };
  LONGLONG getMSec2(void)
  {
     diff.QuadPart=(t3.QuadPart-t2.QuadPart)*((LONGLONG)1000000)/freq.QuadPart;
     return diff.QuadPart;
  };

  AnsiString print1(void)
  {
     long l=getMSec1();
     AnsiString s;
     s=s+l;//+" ysec   ";
     return(s.c_str());
  };
  AnsiString print2(void)
  {
     long l=getMSec2();
     AnsiString s;
     s=s+l;//+" ysec   ";
     return(s.c_str());
  };
};


extern Stoppuhr_base Stoppuhr;

void testDaten(int windowsSize,int _anz);

extern Hammock m1;
double * test(RMbase *m);

#else

class Stoppuhr_base
{
public:
  Stoppuhr_base(){ };
  int getFreq(void) { return 0; };
  void start(void) { };
  void stop1(void) { };
  void stop2(void) { };

  long getMSec1(void) {	return 1; };
  long getMSec2(void) { return 1; };

  char * print1(void) { return ((char*)""); };
  char * print2(void) { return ((char*)""); };
};//Stoppuhr; // todo: Fragen was hiermit gemeint ist




#endif

#endif
