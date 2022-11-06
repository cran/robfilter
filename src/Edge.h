#ifndef __EDGE_H__
#define __EDGE_H__

#include "Line.h"

#define DOUBLE_NICHT_SPEICHERN 1

//INF is defined in Line.h

class FaceIter;
//enum Direction {PRE=0,NEXT=1};

class Edge
{
  friend class FaceIter;
  // the two edges
  Edge *kanteA;
  Edge *kanteB;

  Line *line;

  //int nr;
  #if(DOUBLE_NICHT_SPEICHERN==0)
  double x_L,x_R;
  #endif

  // Suppose the edges of the current cell are traversed in clockwise order

  // if the arrowhead of this edge and the arrowhead of the edge kanteA
  // joining the same intersection point, then directionA==1.

  // if the arrowhead of this edge and the stub end of the edge
  // joining the same intersection point, then directionA==0.

  // if the stub end of this edge and the arrowhead of the edge kanteB
  // joining the same intersection point, then directionB==0.

  // if the stub end of this edge and the stub end of the edge kanteB
  // joining the same intersection point, then directionB==1.

  unsigned int directionA:1;
  unsigned int directionB:1;
public:
  unsigned int erasewithdelete:1;

public:

  void Init(void)
  {
    kanteA=0;
    kanteB=0;
    directionA=0;
    directionB=0;
    line=0;
    erasewithdelete=1;
  }

  Edge(void)
  {
    Init();
  }

  ~Edge(void)
  {
  }

  void setLine(Line *l)        {  line=l; }
  Line * getLine(void)         {  return line;  }

  // Suppose the edges of the current cell are traversed in clockwise order
  // if the current edge is directed anticlockwise, then thisDir==1
  // if the next edge is directed anticlockwise, then eDir==1

  Edge * getNext(int thisDir, int *eDir)
  {
       
     if (!thisDir)
     {
       *eDir=directionB;
       return kanteB;
     }
     else
     {
       *eDir=1-directionA;
       return kanteA;
     }
  }

  Edge * getPre(int thisDir,int *eDir)
  {
     if (thisDir==0)
     {
       *eDir=1-directionA;
       return kanteA;
     }
     else
     {
       *eDir=directionB;
       return kanteB;
     }
  }

  // Je nach dem, aus welcher Richtung eine Kante gefunden wird
  // zeigt sie gedanklich in verschieden Richtung
  // thisDir und eDir geben mit =1 an, ob die tatsaechliche Richtung
  // umgekehrt ist.
  // antiDir=0 bedeuted, das ein Pfeilende und ein Pfeilanfang
  // aufeinandertreffen

  void setNext(int thisDir,Edge *e,int eDir,int antiDir)
  {
       
     if (thisDir==0)
     {
       kanteB=e;
       directionB=(thisDir+antiDir+eDir)%2;
     }
     else
     {
       kanteA=e;
       directionA=(thisDir+antiDir+eDir)%2;
     }
  }
  
  void setPre(int thisDir,Edge *e,int eDir,int antiDir)
  {
       
     if (thisDir==0)
     {
       kanteA=e;
       directionA=(thisDir+antiDir+eDir)%2;
     }
     else
     {
       kanteB=e;
       directionB=(thisDir+antiDir+eDir)%2;
     }
  }

  Edge *getNextStrand(int *lowerPoint)
  {
    int d;
    Edge *f;
    f=getNext(0,&d);
    *lowerPoint=0;

    if (f!=0 && d!=0)
    {
      f=f->getPre(d,&d);
      f=f->getPre(d,&d);
      *lowerPoint=1;
    }
    return(f);
  }

  #if(DOUBLE_NICHT_SPEICHERN==1)
  void set_X_pre(double x,int dir) { };
  void set_X_next(double x,int dir) { };
  double get_X_pre(int dir)
  {
    int dummy;
    double w;

    Edge *e = getPre(dir, &dummy);
    if (e)
    {
      Line *quer = e->getLine();
      if (quer->isNormal() && line->isNormal())
        w = line->schnittX(quer);
      else
      {
        if (quer->isRight() || line->isRight())
          w = INF;
        else
          w = -INF;
      }
    }
    else
    {
      if (line->isRight())
        w = INF;
      else
        w = -INF;
    }
    return w;
  }
  
  double get_X_next(int dir)
  {
    int dummy;
    double w;

    Edge *e = getNext(dir, &dummy);
    if (e)
    {
      Line *quer = e->getLine();
      if (quer->isNormal() && line->isNormal())
        w = line->schnittX(quer);
      else
      {
        if (quer->isRight() || line->isRight())
          w = INF;
        else
          w = -INF;
      }
    }
    else
    {
      if (line->isRight())
        w = INF;
      else
        w = -INF;
    }
    return w;
  }

  #else

  void set_X_pre(double x,int dir)
  {
     if (dir==1)
       x_R=x;
     else
       x_L=x;
  }
  
  void set_X_next(double x,int dir)
  {
     if (dir==0)
       x_R=x;
     else
       x_L=x;
  }

  double get_X_pre(int dir)
  {
     if (dir==1)
       return (x_R);
     else
       return (x_L);
  }
  
  double get_X_next(int dir)
  {
     if (dir==0)
       return (x_R);
     else
       return (x_L);
  }
  #endif

};

#endif
