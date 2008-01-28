
// Programmed by Thorsten Bernholt
// University of Dortmund

//---------------------------------------------------------------------------

#ifndef medianFilterH
#define medianFilterH

#include "CircularArray.h"
#include "RMquick.h"

class MedianFilter
{
  Median median;
  CircularArray<double> fenster;

  double *temp;

public:

  MedianFilter(int windowWidth);
       
  ~MedianFilter();
  
  void add(double wert);

  void remove();

  double getMedian();
};


//---------------------------------------------------------------------------
#endif
