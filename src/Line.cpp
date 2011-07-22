#include "hammock.h"
#include "Line.h"
#include <R_ext/Arith.h> // for ISNA and NA_REALs

double Line::schnittX(Line *f)
   {
     return(f->b - b)/(f->m - m );
   }
   
double Line::schnittY(Line *f)
   {
     return ( m * f->b - f->m * b )/( m - f->m );
   }

double Line::getMedian(int anzLines)
{
  if (median_)  // Test ob der median schon berechnet wurde
  {
    
	const double right = median_->get_X_next(0);
  
    /* falls die Anzahl der Linien ungerade ist gibt
     * es eine gerade Anzahl von Schnittpunkten, sodass
     * der Median als Durchschnitt der beiden Werte
     * in der Mitte berechnet wird 
     * */ 
    if (anzLines & 1) {
		const double left = median_->get_X_pre(0);
		return R_FINITE(left) ? (left + right) / 2.0 : right; 
	}
    else
      return right; // ansonsten wird einfach der Median genommen
  }
  else {
    return NA_REAL;
  }
}

   void Line::addSchnitt(Edge *e)
   {
     if (this)
     {
       mark=1;
       if (!median_)
       {
         int med = MEDIAN_LINKS(root->getAnzLines() - 1);
         if (links == med)
           median_ = e;
         else
           links++;
       }
       else
       {
         int d;
         Line *lq = median_->getNext(0, &d)->getLine();
         if ((m > lq->m && lq->mark == 0) || (m < lq->m && lq->mark == 1))
            links++;
          else
            rechts++;
       }
     }
   }

   void Line::delSchnitt(void)
   {
     if (this)
     {
       int d;
       mark = 1;

       if (median_)
       {
         Line *lq = median_->getNext(0,&d)->getLine();
         if ((m > lq->m && lq->mark == 1) || (m < lq->m && lq->mark == 0))
           links--;
         else
           rechts--;
       }
     }
   }

   void Line::geheNachLinks(void)
   {
     int e_dir,m_dir;

     Edge *e = median_->getPre(0, &e_dir);
     median_ = e->getPre(e_dir, &m_dir);
     links--;
     rechts++;
   }
   
   void Line::geheNachRechts(void)
   {
     int e_dir,m_dir;

     Edge *e = median_->getNext(0,&e_dir);
     median_ = e->getPre(e_dir,&m_dir);
     links++;
     rechts--;
   }

   void Line::updateMedian(void)
   {
      if (median_)
      {
        const int l = MEDIAN_LINKS(root->getAnzLines()-1);
	//	std::cout << "Anzahl Lines:" << root->getAnzLines() << "\n";
    //    std::cout << "l:" << l << " links:" << links << " rechts:" << rechts << "\n";
        
        while (l < links)
          geheNachLinks();

        while (l > links)
          geheNachRechts();

	//	std::cout <<"----------------------\n";

    //    std::cout << "l:" << l << " links:" << links << " rechts:" << rechts << "\n";
        
      }
   }

