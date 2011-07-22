#include "RegLine.h"
#include <R_ext/Arith.h> // for ISNA and NA_REAL

RegLine::RegLine(double intercept, double slope, double fit) :
 y_achse(intercept), steigung(slope), fitness(fit) {}
 
 RegLine RegLine::nullLine = RegLine(NA_REAL, NA_REAL, NA_REAL);
//Note: If we want to avoid requiring the R libraries here,
//we could use numeric_limits<double>::quiet_NaN() from the limits
//header in case R is not present.
