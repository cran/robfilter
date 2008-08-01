//---------------------------------------------------------------------------
// Programmed by Thorsten Bernholt
// University of Dortmund

#include "MedianFilter.h"


MedianFilter::MedianFilter(int windowWidth) : fenster(windowWidth), temp(new double[windowWidth])
{
}
       
MedianFilter::~MedianFilter()
{
  delete []temp;
}
  
void MedianFilter::add(double wert)
{
 fenster.append(wert);
}

void MedianFilter::remove()
{
 fenster.removeOldest();
 }

double MedianFilter::getMedian()
{
 fenster.copyInto(temp);
 return median.getMedian(temp, fenster.size());
}
