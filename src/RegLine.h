#ifndef __REGLINE_H__
#define __REGLINE_H__

class RegLine
{
  public:
    /**
     * Value returned in case there is no sensible numeric
     * value. When linking with R, this is equivalent to
     * setting every field to NA_REAL.
     */
    static RegLine nullLine;
    
    /**
     * Intercept
     */
    double y_achse;
    
    /**
     * Slope
     */
    double steigung;
    
    /**
     * For LMS
     */
    double fitness;//nur LMS
    
    RegLine(double intercept = 0.0, double slope = 0.0, double fitness = 0.0);
};

#endif
