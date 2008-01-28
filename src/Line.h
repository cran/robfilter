#ifndef __LINE_H__
#define __LINE_H__

class Hammock;
class Edge;

#include <limits>

const double INF = std::numeric_limits<double>::infinity();

class Line
{
   //friend class Hammock;

   int links2,rechts2; // werden gesetzt beim Einfügen und Löschen
   int typ; // 0==normal, 1==left, 2==right
   
   int links,rechts; //Number of intersections on both sides of the median
   Edge *median_;    // der Pointer auf die median Kante   
public:
   Hammock *root;
   int nr; //used by regDepth

   double m; //slope of the line
   double b; //y-axis

   Edge *startE;

   int mark;
     //During Delete/Insert:
     //   =0 until now the inserted/deleted Line
     //      has not passed this line here.
     //   =1 the inserted/deleted Line
     //      has passed this line.

   double schnittX(Line *f);
   double schnittY(Line *f);

   void setTypLeft(void)  {  typ=1; }
   void setTypRight(void) {  typ=2; }
   bool isRight(void)      {  return typ==2; }
   bool isNormal(void)     {  return typ==0; }
   double getMedian(int anzLines);
   void attemptToDeleteOrDissect(Edge *e1,Edge *e2)
   {
      if (e1==median_ || e2==median_)
        geheNachLinks();

        //XXX! Why do the same two lines appear twice?

      if (e1==median_ || e2==median_)
        geheNachLinks();
   }
   
   void addSchnitt(Edge *e);
   void delSchnitt(void);
   void geheNachLinks(void);   // verschiebt den Zeiger auf die Kante um eins nach links
   void geheNachRechts(void);

   void updateMedian(void);

   void setMedian(Edge *e)
   {
     median_=e;
   }
   void clearCounts(void)
   {
     links2=0;
     rechts2=0;
   }

   // After calling clearCounts and inserting a line,
   // this function decides if the deleted line was above the median edge
   bool insertedLineWasAbove(void)
   {
     if (rechts2==1)
       return true;

     return false;
   }
   
   bool deletedLineWasAbove(void)
   {
     if (links2==1)
       return true;

     return false;
   }

   ~Line(){}
   Line(double _m, double _b, Hammock *H)
   {
     root=H;
     m=_m;
     b=_b;
     typ=0;
     median_=0;
     links=0;
     rechts=0;
     startE=0;
   }
};

#endif
