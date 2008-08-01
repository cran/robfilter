// Programmed by Thorsten Bernholt
// University of Dortmund

#include "RMquick.h"
//#include <iostream>

int MEDIAN_LINKS(int anz)
{
  //XXX! For anz >= 0, this is always equal to anz/2
  //(using integer math). And anz is always >= 0.
    
  int m=(int)ceil((anz+1.0)/2.0)-1;
  if (m<0)
    m=0;
  
  return (m);
};

/*int MEDIAN_RECHTS(int anz)
{
	int m = (int)floor(((double)anz)/2.0)+1;
	
	if (m > anz -1) 
	  m=anz-1;
	  
	return m;
};*/


#ifndef BORLAND_WIN
int random(int a)
{
  return (int)((((double)rand())/((double)RAND_MAX))*a);
};
#endif
