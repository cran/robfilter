// Es wird in keinsterweise garantiert, dass dieses Programm fehlerfrei ist.
// Es ist daher nicht für den Einsatz in sicherheitskritischen Anwendungen
// geeignet. Recycle allocated something using new

#include "hammock.h"
#include <iostream>
#include <R_ext/Arith.h> // for ISNA and NA_REALs


  RegLine Hammock::getRM(double timeZero)
  {
    /* The Repeated Median Algorithmn needs at least 4 lines to run correctly
     * if the number of lines fall below 4 the hammock must be reInitiated to
     * work with the getRM function. This is done by calling reInitHammock once
     * the number of lines increases 4 again.
     * The reInitHammock function first memorizes all lines in the Hammock (note not more than 3)
     * than an reInit is Done by deleting the hammock
     * at last the Hammock is reInit and the lines are added to it too.
     * */
	//std::cout << "start getRM anzLines:" << anzLines << "\n";
	if (anzLines < 5){
   		if (!initDone){
   			needReInit = true;
   		}
   		return RegLine::nullLine;
    }else {
    	if (needReInit){
    		reInitHammock();
    	}
    	initDone = false;
    }

    for (int i=0; i<anzLines; i++) {
      Line* l = lineTab->get(i);
      l->updateMedian();
      medTab[i] = l->getMedian(anzLines);
    }

    double *tab2 = new double[anzLines];


    const double steigung = med.getMedian(medTab, anzLines);

    for (int i =0; i < anzLines; i++)
    {
      Line *g = lineTab->get(i);
      tab2[i] = g->b - steigung * g->m;
    }

    //Auf Level umrechnen
    const double y_achse = med.getMedian(tab2, anzLines);

    delete[] tab2;

    return RegLine(y_achse, steigung);

  }
void Hammock::updateRepeatedMedian()
{
	//std::cout << "start updateRepeatedMedian anzLines:" << anzLines << "\n";
	if (anzLines < 5){
   		if (!initDone){
   			needReInit = true;
   		}
   		return;
    }else {
    	if (needReInit){
    		reInitHammock();
    	}
    	initDone = false;
    }



    for (int i=0; i<anzLines; i++) {
      	Line* l = lineTab->get(i);
  		l->updateMedian();
    }
}

void Hammock::reInitHammock(){


	int n = lineTab->size();
	Line ** mem = new Line*[n];
	for (int i = 0; i < n; i++){
		mem[i] = lineTab->get(0);
		lineTab->removeOldest();
	}

	 bin.freeUsedMemory();
	 L=bin.neu(); //new Edge();
	 R=bin.neu(); //new Edge();

	 const double x = -INF;

	 Lup=L;

	 L->setLine(border_L);
	 L->set_X_pre(x,0);
	 L->set_X_next(x,0);

	 R->setLine(border_R);
	 R->set_X_pre(x,0);
	 R->set_X_next(x,0);

	 needReInit = false;
	 initDone = true;
	 anzLines = 0;
	// delete mem;
	 for (int i = 0; i < n; i++){
	 	mem[i]->resetLine();
		addLine(mem[i]);
	 }

}

Edge *Hammock::dissect_L(Line *neueLinie)
{
  Edge* altL=L;
  Edge* neuL=bin.neu();
  L=neuL; // neue untere Kante der linken Begrenzung
  Edge* out=bin.neu();  // neue Kante die auf der linken Seite der linken Begrezung rausguckt
  Edge* lose = bin.neu(); // neue Kante die auf der rechten Seite der linken Begrenzung rausguckt und fuer
  // die in den naechsten Schritten der Test mit den vorhandenen Kanten durchgefuehrt werden muss

  neueLinie->startE = out;

  altL->setNext(0,out,0,1);
  out->setNext(0,neuL,0,0);
  neuL->setPre(0,lose,0,1);
  lose->setPre(0,altL,0,0);

  neuL->setLine(altL->getLine());
  out->setLine(neueLinie);
  lose->setLine(neueLinie);

  const double x = -INF;

  altL->set_X_next(x,0);
  neuL->set_X_pre(x,0);
  neuL->set_X_next(x,0);
  out->set_X_next(x,0);
  out->set_X_pre(x,0);
  lose->set_X_pre(x,0);

//  d1=L;
  return lose;
}

void Hammock::dissect_R(Edge *lose,Line *neueLinie)
{
  Edge *altR, *neuR, *out;
  altR=R;
  neuR=bin.neu();
  R=neuR;
  out=bin.neu();

  altR->setNext(0,out,0,0);
  out->setPre(0,neuR,0,1);
  neuR->setPre(0,lose,0,0);
  lose->setNext(0,altR,0,1);

  neuR->setLine(altR->getLine());
  out->setLine(neueLinie);
  lose->setLine(neueLinie);

  double x=INF;

  altR->set_X_next(x,0);
  neuR->set_X_pre(x,0);
  neuR->set_X_next(x,0);
  out->set_X_next(x,0);
  out->set_X_pre(x,0);
  lose->set_X_next(x,0);

  neueLinie->addSchnitt(lose);//,x);
}

// Gibt neue lose Kante zurück

Edge * Hammock::dissectEdge(Edge *lose, Edge *e, int e_dir)
{
   Edge * e2,  *  neu_lose, *  a1, *   a2,   * a3;
   int    e2_dir,              a1_dir, a2_dir, a3_dir;

   e->getLine()->attemptToDeleteOrDissect(e, NULL);

   e2 = bin.neu();
   neu_lose = bin.neu();

   a1 = e->getNext(e_dir,&a1_dir);
   a2 = a1->getPre(a1_dir,&a2_dir);
   a3 = a2->getPre(a2_dir,&a3_dir);

   const double x = calcSchnitt(lose->getLine(), e->getLine());

   //Falls die Kante e falsch herum durchlaufen wurde
   //muss diese Richtung auf e2 übertragen werden

   //Notation from Bernholt & Fried 2003
   // f_i   = lose
   // f_i+1 = neu_lose
   // e_i   = e
   // e_i'  = e2
   // a_i'  = a1

   e2_dir = e_dir;

   lose->setNext(0, e2, e2_dir, 0);
   e2->setNext(e2_dir, a1, a1_dir, 0);
   a3->setPre(a3_dir, e2, e2_dir, 0);
   e2->setPre(e2_dir, neu_lose, 0, 1);
   neu_lose->setPre(0, e, e_dir, 0);
   e->setNext(e_dir, lose, 0, 1);

   // set the new coordinates
   const double e_x_next = e->get_X_next(e_dir);

   e2->set_X_next(e_x_next, e2_dir);
   e->set_X_next(x, e_dir);
   lose->set_X_next(x, 0);
   neu_lose->set_X_pre(x, 0);
   e2->set_X_pre(x, e2_dir);

   e2->setLine(e->getLine());
   neu_lose->setLine(lose->getLine());

   lose->getLine()->addSchnitt(lose);//,x);
   e2->getLine()->addSchnitt(e2);//,x);

   return neu_lose;
}

void Hammock::computeLXX(void)
{
  int d;
  double y1,y2;
  LMS.init();
  LTS.init();

  Edge *o,*u;
  Edge *o2,*u2;
  int lowP_o=0,lowP_u=0;

  double ow,uw,w;
  w = 0; /*initialize*/
  Line *l1,*l2;

  //XXX! Is it right to do nothing here when k > anz?
  const int anz=getAnzLines();
  const int k=h-1;
  const int numSubsets = anz - k;
  for(int i=0; i<numSubsets; i++)
  {
    LTS.clear();

//    if (outD!=0)  fprintf(outD,"Schleife %i\n",i);

    o = lineTab->get(i)->startE;
    u = lineTab->get(i + k)->startE;
    for (int j = 0; j != h; ++j) {
        const Line* line = lineTab->get(i + j);
        LTS.ins(line->m, line->b);
    }

    o=o->getNext(0,&d);
    o=o->getPre(d,&d);

    u=u->getNext(0,&d);
    u=u->getPre(d,&d);

    o2=o->getNextStrand(&lowP_o);
    u2=u->getNextStrand(&lowP_u);

    ow=o->get_X_next(0);
    uw=u->get_X_next(0);

    bool Schnittpunkt_vorhanden=false;
    bool swap_o=true;

    while(o!=0 && u!=0 && (ow<INF || uw<INF))
    {
      // Segment bearbeiten

      Schnittpunkt_vorhanden=false;

      // Choose that intersecting point with smaller x-coordinate
      if (ow<=uw && o2!=0)
      {
        w=ow;
        Schnittpunkt_vorhanden=true;
        swap_o=true;
      }
      if (uw<ow && u2!=0)
      {
        w=uw;
        Schnittpunkt_vorhanden=true;
        swap_o=false;
      }


      if (Schnittpunkt_vorhanden)
      {
        l1=o->getLine();
        l2=u->getLine();

        y1=l1->m*(-w)+l1->b;
        y2=l2->m*(-w)+l2->b;


        //LTS -----------------------

//        if (outD!=0)  fprintf(outD,"swap_o = %i    lowP_o = %i    lowP_u = %i  ",(int)swap_o,lowP_o,lowP_u);

        if (swap_o)
        {
          l1=o->getLine();
          l2=o2->getLine();
        }
        else
        {
          l1=u->getLine();
          l2=u2->getLine();
        }

//        if (outD!=0) fprintf(outD," w= %f ",w);

        //if (outD)  fprintf(outD," schnitt: % f ",l1->schnittX(l2));

        if ((swap_o && lowP_o==0) || (swap_o==false && lowP_u==1))
        {
          // A new point is discovered, but there are k+1 points in the set
          LTS.del( l1->m , l1->b );
          LTS.ins( l2->m , l2->b );
        }
        else
        {
          //exactly k points inside
          LTS.keep_min( w , y1,y2,this);
          LMS.keep_min( w , y1,y2,this);

          /*
          if (swap_o)
            {if (outD!=0) fprintf(outD,"     Teste(new): (%.0f) (%.0f) (%.0f)  \n",l1->m,l2->m,u->getLine()->m);}
          else
            {if (outD!=0) fprintf(outD,"     Teste(new): (%.0f) (%.0f) (%.0f)  \n",l1->m,l2->m,o->getLine()->m);}
            */
        }

//        if (outD!=0)  fprintf(outD,"\n");
      }

      // Find next intersecting point

      if (ow<=uw)
      {
        o=o2;
        o2=o->getNextStrand(&lowP_o);
        ow=o->get_X_next(0);
      }
      else
      {
        u=u2;
        u2=u->getNextStrand(&lowP_u);
        uw=u->get_X_next(0);
      }

    }
  }
};

/*
void Hammock::testMedian(){
	Hammock m;
	int anz_Linien = 5;
	m.init(anz_Linien);
	m.addPunkt(0,0);
	m.addPunkt(1,-1);
	m.addPunkt(2,-4);
	m.addPunkt(4,-12);
	m.addPunkt(5,-20);

	//m.addPunkt(2,-7);
	m.printLines();
	//m.print();
 	double median_test = m.lineTab[0]->getMedian(anz_Linien);
 	double gesamt_median_test = m.getRM(0)->steigung;
 	std::cout << "Median ist " << median_test << std::endl;
 	std::cout << "Median gesamt ist " << gesamt_median_test << std::endl;
}*/

// Die einzufügende Linie muss eine größere Steigung haben als alle vorhandenen
int Hammock::addPunkt(double m,double b)
{
  // Neue Linie erzeugen
  return addLine(new Line(m,b,this));
}

/** Remove the oldest point added that has not
  * yet been removed. If there is no such point,
  * do nothing.
  */
void Hammock::removePunkt()
{
     if (anzLines < 1) {
        return;
     }

     delLine();
}

int Hammock::addLine(Line *neuL)
{
  neuL->root = this;


  if (anzLines >= windowSize){
	  delLine();
  }
	anzLines++;

  lineTab->append(neuL);

  // In linke Begrenzung einfügen
  Edge *lose = dissect_L(neuL);

  // Flags löschen

  L->getLine()->mark=1;

  for (int i = 0; i<anzLines; i++)
  {
   	lineTab->get(i)->mark = 0;
  }

  iter.start(lose, 1);
  int test_dir;
  Edge *test = iter.next(&test_dir);



  double last_u = 0;
  bool init = true;

  while(test != 0 && test->getLine() != 0 && !test->getLine()->isRight()) //bis bei R angekommen
  {
  	  // calculate intersection
      const double u = neuL->schnittX(test->getLine());
	  // Is the edge "test" intersected by the new line
	  bool leq;
	  bool geq;

      double u1 = test->get_X_pre(0);
      double u2 = test->get_X_next(0);


      leq = (test->get_X_pre(0) < u) || (u1-u==0);
      geq = (u < test->get_X_next(0)) || (u2-u==0);

      //

      /*if (debugInfo){
      	//printf("\npre %1.20f, u %1.20f\n",test->get_X_pre(0),u);
      	//printf("u %1.20f, next %1.20f\n",u,test->get_X_next(0));
      	//printf("leq %i , geq %i\n",leq,geq);

      }*/
      if (leq  && geq)
      {
        if (last_u > u && !(last_u-u==0) && !init) {
        	return -1; // Error in calculation. Please use an other algorithmn
        }
        init = false;
        last_u = u;
        lose = dissectEdge(lose, test, test_dir);
        iter.invertDirection();
      }

      test = iter.next(&test_dir);
  }

  //std::cout << "after:\n";
  //checkLines();

  dissect_R(lose,neuL);
  return 0;
}

void Hammock::delLine(void)
{

  // Flags löschen
    L->getLine()->mark=0;
  //std::cout<< "get Last Line done" << std::endl;

  for (int i=0; i<anzLines; i++) {

    lineTab->get(i)->mark=0;

  }


  int dir;
  Edge *pre = Lup->getNext(0,&dir);
  Line *delL = pre->getLine(); // #### absturz stelle

  anzLines--;
  delL->setMedian(0); //XXX! Is this necessary?

   //Notation from Bernholt & Fried 2003
   // f_i   = lose
   // f_i+1 = neu_lose
   // e_i   = e
   // e_i'  = e2
   // a_i'  = a1

  int   next_dir, e_dir,e2_dir,a1_dir,a2_dir,a3_dir;
  Edge *next,*e,*e2,*a1,*a2,*a3;

  while (pre)
  {
     e2=pre->getNext(0,&e2_dir);
     next=e2->getPre(e2_dir,&next_dir);
     e=next->getPre(next_dir,&e_dir);

     a1=e2->getNext(e2_dir,&a1_dir);

     if (a1)
     {
       a2=a1->getPre(a1_dir,&a2_dir);
       a3=a2->getPre(a2_dir,&a3_dir);

       Line* l=e->getLine();

       l->attemptToDeleteOrDissect(e,e2);

       e->set_X_pre(e2->get_X_next(e2_dir),e_dir);

       e->setPre(e_dir,a1,a1_dir,1);
       a3->setPre(a3_dir,e,e_dir,1);

       l->delSchnitt();
     }
     else
     {
      e->setPre(e_dir,0,0,0);
       bin.shred(next);
       next=0;
     }
     bin.shred(pre);
     bin.shred(e2);

     pre=next;
  }

  lineTab->removeOldest();

  delete delL;


  if (anzLines == 0){
	 reInitHammock();
  }


}



  void Hammock::computeRegDepth(void)
  {
  //XXX! Copying the contents of lineTab might save us some time here by avoiding
  //index checks/dealing with indexes wrapping around the end of the array.

    regDepth.clear();

    border_L->nr=-1;
    border_R->nr=-1;

    for(int i=0;i<anzLines;i++)
      lineTab->get(i)->nr=i;

    Edge *e;
    int value;

    const int twiceAnzLines = 2 * anzLines;

    SegmentTree t(twiceAnzLines);

    for(int i=0; i < twiceAnzLines; i++)
      t.extend(0);
    for(int k=0;k<anzLines;k++)
    {

      for(int i=0; i<twiceAnzLines; i++)
        t.set(i,0);

      for(int i=0;i<anzLines;i++)
      {
        lineTab->get(i)->mark=0;
        t.add(i + 1, i + anzLines, 1);
        //The loop condition guarantees that i is at most anzLines-1.
        //So, the maximum value for i + 1 is anzLines. If anzLines == 0, this
        //loop is not executed. If anzLines >= 1, then twiceAnzLines > anzLines.
        //As well, the maximum value of i + anzLines is then 2 * anzLines - 1.
        //Thus we may safely avoid taking these values modulo twiceAnzLines,
        //as was done before here.
      }

      for(int i=0;i<k;i++)
      {
        lineTab->get(i)->mark=1;

        t.add(i+1, i + anzLines, -1);
        t.add(i + anzLines + 1, i, 1);
        //Again the loop condition guarantees that i + anzLines and i itself
        //will always be less than twiceAnzLines (since k is also always <= anzLines - 1).
        //Further since i < k < anzLines, that is,
        //the maximum value for i is anzLines - 2, i + 1 is also always < twiceAnzLines.
        //This also means that the maximum of i + anzLines + 1 is 2 * anzLines - 1.
      }

      // We are on the line k
      t.add((anzLines+k+1)%(2*anzLines), k, 1);
      //k is guaranteed to be less than anzLines by the loop condition,
      //which in turn means that k is always < twiceAnzLines, thus we
      //can safely avoid taking k % twiceAnzLines.


      e=lineTab->get(k)->startE;

      //Jump over L;
      e=regDepth.nextEdgeOnLine(e);

      //Jump over next line
      e=regDepth.nextEdgeOnLine(e);

      while(e!=0 && regDepth.linePassed->nr!=-1)
      {
        int nr=regDepth.linePassed->nr;
        // enter the line
        if (regDepth.linePassed->mark==0)
          t.add((anzLines+nr+1)%twiceAnzLines,(nr)%twiceAnzLines,1);
        else
          t.add((nr+1)%twiceAnzLines,(nr+anzLines)%twiceAnzLines,1);

        t.getMinPos(0, twiceAnzLines - 1);
        value=t.getLastValue();

        Line* lineK = lineTab->get(k);
        double m = lineK->schnittX(regDepth.linePassed);
        double b = lineK->schnittY(regDepth.linePassed);

        regDepth.keepMin(m,b,value);

        //exit the line
        if (regDepth.linePassed->mark==0)
          t.add((nr+1)%twiceAnzLines,(nr+anzLines)%twiceAnzLines,-1);
        else
          t.add((anzLines+nr+1)%twiceAnzLines,(nr)%twiceAnzLines,-1);
        regDepth.linePassed->mark=1-regDepth.linePassed->mark;

        e=regDepth.nextEdgeOnLine(e);
      }
    }
  }


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

#include <vcl.h>

double *X,*Y;
int len,anz;
int ausg=1;

void testDaten(int windowsSize,int _anz)
{
  len=windowsSize;
  anz=_anz;
  X=new double[anz+len+2];
  Y=new double[anz+len+2];
  int st=1;
  for(int i=0;i<anz+len+1;i++)
  {
    X[i]=st++;
    Y[i]=random(10000);
  }
  if (ausg==1)
  {
    ausg=0;
    char *buf=new char[100];
    sprintf(buf,"Fenster=%i   Updates=%i\r\n",len,anz);
    addMemo(buf);
  }
}
extern Hammock m1;
double * test(RMbase *m)
{
  m->init(len);

  // Anlaufphase

  int st=0;

  double *med=new double[anz+len];
  RegLine l;

  for(int i=0;i<len;i++)
  {
     m->addPunkt(X[st],Y[st]);
     //l=m->getRM();
     //med[st]=l->steigung;
     st++;
  }

  Stoppuhr.start();

  for(int i=0;i<anz;i++)
  {
    m->addPunkt(X[st],Y[st]);
    //l=m->getRM();
    l=m->getRM(1);
    med[st]=l->steigung;
//    med[st]=m->getSlope();
//    med[st]=m->getLMS();

  /*
  {
    for (int i=0;i<len;i++)
      m1.lout(i);
    addMemo("\r\n");
  }
  */

    st++;
  }

//  double med2=m->getSlope();
  Stoppuhr.stop1();
  char *buf=new char[100];
//  sprintf(buf,"%s:  %s   ",m->getName(),Stoppuhr.print1().c_str());
  sprintf(buf," %s",Stoppuhr.print1().c_str());
  addMemo(buf);
  return med;
};
#else

#endif

Stoppuhr_base Stoppuhr;

class debug_test_class
{
public:
  double *med;
  double schnitt;
  int anz;
  int count;

  void putRM(double *medTab,int anzahl)
  {
    anz=anzahl;
    delete[]med;
    med=new double[anz];

    for(int i=0;i<anz;i++)
        med[i]=medTab[i];

    for(int j=0;j<anz;j++)
      for(int i=1;i<anz;i++)
        if (med[i-1]<med[i])
          {double t=med[i-1];med[i-1]=med[i];med[i]=t;};
  };
  void putLMS(double s)
  {
    int nichts=1;
    count++;
    schnitt=s;
    for(int i=0;i<anz;i++)
      if (med[i]==s)
        {
          nichts=0;
        }
  };
  debug_test_class(void)
  {
    count=0;
    med=0;
  };
}debug_test;


void Hammock::print(void)
{
/*
  Hammock *m=this;
  Edge *e;
  char *buf=new char[100];
  */
  /*
 if (DEBUG)
 {
  for(int i=0;i<Edge::count;i++)
  {
    if (i!=0)
      e=Edge::tab[i];
    else
      e=m->L;
    if (e!=0)
    {
      e->setAlias(0);

      int typ=-1;
      double _m=0;
      if (e->getLine()!=0)
        _m=e->getLine()->m;



      if (e->getLine()!=0)
        typ=e->getLine()->typ;

      sprintf(buf,"Kante %3i: Pre=%3i   AntiPre=%3i   Next=%3i   AntiNext=%3i  typ=%3i      line=%3.0f\r\n",
        e->getNr(),
        e->getA(PRE)->getNr(),
        e->getA(ANTIPRE)->getNr(),
        e->getA(NEXT)->getNr(),
        e->getA(ANTINEXT)->getNr(),
        typ,
        //e->getPreVertex()->getNr(),
        //e->getNextVertex()->getNr(),
        _m
      );
      addMemo(buf);
    }
  }
  */
  /*
  Vertex *v;
  for(int i=1;i<Vertex::count;i++)
  {
    if (i!=0)
      v=Vertex::tab[i];
    if (v!=0)
    {
      v->setAlias(0);

      Line *l1,*l2;
      float m1=0,m2=0;
      l1=v->getLine();
      l2=v->getAntiLine();
      if (l1!=0)
        m1=l1->m;
      if (l2!=0)
        m2=l2->m;


      sprintf(buf,"Knoten %3i: Pre=%3i   Next=%3i   Left=%3i  Right=%3i  line=(%3.0f,%3.0f)  x=%3f\r\n",
        v->getNr(),
        v->getPre()->getNr(),
        v->getNext()->getNr(),
        v->getLeft()->getNr(),
        v->getRight()->getNr(),
        m1,m2,
        (float)v->x
      );
      addMemo(buf);
    }
  }
   */
  /*
  for(int i=1;i<Edge::count;i++)
  {
    e=Edge::tab[i];
    if (e!=0)
    {
      e->setAlias(0);

      int typ=-1;
      double _m=0;
      if (e->getLine()!=0)
        _m=e->getLine()->m;

      if (e->getLine()!=0)
        typ=e->getLine()->typ;

      sprintf(buf,"Kante %3i:    xL=%3.1f      xR=%3.1f \r\n",
        e->getNr(),
//        e->getPreVertex()->getNr(),
//        e->getNextVertex()->getNr(),
        e->x_L,
        e->x_R
      );
      addMemo(buf);
    }
  }
  */
          /*
//  Vertex *v;
  for(int i=1;i<Vertex::count;i++)
  {
    if (i!=0)
      v=Vertex::tab[i];
    if (v!=0)
    {
      v->setAlias(0);

      Line *l1,*l2;
      float m1=0,m2=0;
      l1=v->getLine();
      l2=v->getAntiLine();
      if (l1!=0)
        m1=l1->m;
      if (l2!=0)
        m2=l2->m;


      sprintf(buf,"Knoten %3i: x=%3.1f\r\n",
        v->getNr(),
        (float)v->x
      );
      addMemo(buf);
    }
  }

  */

//}

/*
  sprintf(buf,"%i\r\n",getAnzLines());
  addMemo(buf);

  Line *l=lineList;


  while(l!=0)
  {
    l->updateMedian();
    sprintf(buf,"(%3.0f) Median %3f vertex= %i  [%i %i] \r\n",l->m-1,l->getMedian(),l->median_->getNr(),l->links,l->rechts );
    addMemo(buf);
    l=l->next;
  };
  */

}
