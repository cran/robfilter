#ifndef BACKPORTS_H
#define BACKPORTS_H

#include <Rinternals.h>
#include <Rversion.h>

#if R_VERSION < R_Version(4, 5, 0)
#define Rf_isDataFrame(x) Rf_isFrame(x)
#endif

#endif
