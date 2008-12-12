/*
 * Programmed by Robin Nunkesser
 * University of Dortmund
 */

#ifndef COMPILE_AS_EXE
#define COMPILE_AS_DLL 1
#endif

#ifdef COMPILE_AS_DLL

#include "CircularArray.h"
#include <cmath>
#include "hammock.h"
#include "LQDAdvanced.h"
#include "MedianFilter.h"

#include "RegLine.h"
#include <Rinternals.h>    // R einbinden
#include <R_ext/Arith.h> // for ISNA and NA_REAL
#include <R_ext/Error.h> // for error
#include <R_ext/Print.h> // for Rprintf
#include <set>
#include <string>
#include <vector>
using namespace std;

const int numberOfRegressionMethods = 6;
const char* regressionMethodNames[numberOfRegressionMethods] =  {"LQD", "RM", "LMS", "LTS", "DR", "MED"};
// const char* regressionMethodNames[numberOfRegressionMethods] =  {"DR", "LMS", "LQD", "LTS", "MED", "RM"};

double randomGaussianValues(double mean, double Stdev){
//Could also use rnorm from <RMath.h>.
	double v1 = (double)(rand()) / (double)(RAND_MAX);
	double v2 = (double)(rand()) / (double)(RAND_MAX);
	double p;
	do{
		p = (2 * v1 - 1) * (2* v1 - 1) + (2* v2 - 1) * (2* v2 -1);
	} while (p >1);
	double x = (2 * v1 - 1) * sqrt(- log(p) / p );
	return sqrt(Stdev) * x + mean;
}

class RobustReg {

/**
 * Index i in nonNaTracker is true if the ith oldest element in the current
 * is a non-NA element, ie, x for which ISNA(x) is false.
 */
CircularArray<bool>* nonNaTracker;

/**
 * Correctness amounts to ensuring that the data structures
 * contain the same non-NA values as naTracker.
 *
 * Initialisation: Let n = fensterbreite.
 * Say that nonNaTracker has been filled with the
 * bools corresponding to the first n values out of responseVector and
 * that each non-NA value out of the responseVector
 * has been inserted into the various data structures.
 * (Let us take the Hammock H to stand for all of them.)
 *
 * Now let x be the next value to be inserted, ie, that
 * at position n in the responseVector.
 * If the oldest value in nonNaTracker is true/non-NA,
 * then the oldest value in H must be deleted from H.
 * -- If H contains n values and x is not an NA, then this will be taken care of
 *    by the code already in place when x is inserted.
 * -- Otherwise, if the count of non-NA values < n or x is an NA,
 *    we must take care of deleting it ourselves.
 * If the oldest value is an NA,
 * then there is no work to do because there is no corresponding value
 * in H.
 */

LQDAdvanced computeLQDAdvanced;
Hammock H;
const double *responseVector;

long fensterbreite;

long timeC;

//false
const bool addNoise;

//1.0
const double noise;

/**
 * Indicates if there are methods in use which use the Hammock.
 */
bool usingHammock;

/**
 * Indicates if we are doing least quartile difference regression.
 */
bool usingLQD;

//0.0
const double searchExponent;

//2
const char lqdDetermineStart;

//4
const char lqdDescendMethod;

//true
const bool lqdApprox;

//The number of values != NA in the current window.
int numNonNAs;

MedianFilter* medianFilter;

public:

RobustReg() : nonNaTracker(0), responseVector(0),
fensterbreite(0),
timeC(0),
addNoise(false),
noise(1.0),
usingHammock(false),
usingLQD(false),
searchExponent(0.0),
lqdDetermineStart(2),
lqdDescendMethod(4),
lqdApprox(true),
numNonNAs(0),
medianFilter(0)
{
}

~RobustReg()
{
  delete medianFilter;
}

private:
/**
 * Insert the next point in responseVector into the various
 * data structures.
 * Return true if the next value was a good value, ie, not an NA.
 */
bool insertNext(int &errorState)
{
  //std::cout << "entering insertNext() \n";
   //std::cout << "accesing responseVector" << std::endl;
   double w = responseVector[timeC];
   ++timeC;

   const bool isNotNA = !ISNA(w);
   if (addNoise && isNotNA)
   {
     w += randomGaussianValues(0,noise);
   }

  nonNaTracker->append(isNotNA);

   if (isNotNA) {
      ++numNonNAs;
   }
  else {
  	//std::cout << "leaving insertNext() --> NA inserted \n";
    return false;
  }

  //Note that timeC has already been incremented so that for the element
  //at i in the C array, we pass i+1 as the corresponding time, thus
  //implicitly matching the sequence numbers/indexes used in R.

  if (usingHammock) {
	  //std::cout << "adding point to hammock" << std::endl;
	  errorState = H.addPunkt(timeC, w);
  }

  if (usingLQD) {
	  //std::cout << "adding point to lqd" << std::endl;
	  computeLQDAdvanced.addPunkt(timeC, w);
  }

  if (medianFilter) {
	  //std::cout << "adding point to median filter" << std::endl;
	  medianFilter->add(w);
  }

  //std::cout << "leaving insertNext() --> value inserted \n";

  return true;
}

/**
 * Remove the oldest point currently in the various data
 * structures.
 */
void removeOldest()
{
	//std::cout << "entering removeOldest() \n";
  if (usingHammock) {
       H.removePunkt();
  }

  if (usingLQD) {
       computeLQDAdvanced.removePunkt();
  }

  if (medianFilter) {
    medianFilter->remove();
  }
  //std::cout << "leaving removeOldest() \n";
}

void recordEstimate(vector<RegLine*>& results, int row) {
    if (results[lmsIndex] || results[ltsIndex]) {
       H.computeLXX();
    }

    if (results[drIndex]) {
       H.computeRegDepth();
    }

    //We take care of shifting from estimating the intercept
    //to estimating the level in the robustRegression method.
    const double timeZero = 0;

    if (results[rmIndex]) {
       results[rmIndex][row] = H.getRM(timeZero);
    }

    if (results[drIndex]) {
       results[drIndex][row] = H.regDepth.getRegLine(timeZero);
    }

    if (results[lmsIndex]) {
       results[lmsIndex][row] = H.LMS.getRegLine(timeZero);
    }

    if (results[ltsIndex]) {
       results[ltsIndex][row] = H.LTS.getRegLine(timeZero);
    }

    if (results[lqdIndex]) {
       results[lqdIndex][row] = computeLQDAdvanced.getLQD(timeZero);
    }

    if (medianFilter) {
       results[medIndex][row] = RegLine(medianFilter->getMedian());
    }
}


void recordNAs(vector<RegLine*>& results, int row) {


    //We take care of shifting from estimating the intercept
    //to estimating the level in the robustRegression method.


	if (results[rmIndex]) {
		H.updateRepeatedMedian();
    }

    if (results[rmIndex]) {
       results[rmIndex][row] = NA_REAL;
    }

    if (results[drIndex]) {
       results[drIndex][row] = NA_REAL;
    }

    if (results[lmsIndex]) {
       results[lmsIndex][row] = NA_REAL;
    }

    if (results[ltsIndex]) {
       results[ltsIndex][row] = NA_REAL;
    }

    if (results[lqdIndex]) {
       results[lqdIndex][row] = NA_REAL;
    }

    if (medianFilter) {
       results[medIndex][row] = NA_REAL;
    }
}

    public:


/**
 * Indexes of the names of the regression names and corresponding
 * storage for results, if it exists.
 */
enum { lqdIndex, rmIndex, lmsIndex, ltsIndex, drIndex, medIndex};

       //When centre is true,
		//we estimate the last value in each window.
		//When centre is false, we estimate the middle
		//value.
    /**
     * This method returns a vector of pointers to arrays, each of
     * which corresponds to the estimate of the corresponding
     * regression method. If the corresponding method is not included
     * in regressionMethods, that element of the vector is null.
     *
     * Note that this method returns estimate of intercept not level.
     *
     * regressionMethods -- the regression methods to apply.
     *                       Available keys: "RM" -- repeated median, "DR" -- deepest regression,
     *                                    "LMS" -- least median of squares, "LTS" -- least trimmed
     *                                     squares regression, "LQD" -- least quartile difference
     *                                     regression, "MED" -- estimate level with the median of each window
     */
    vector<RegLine*> robustRegression(const double* response,
                              int responseSize,
                              int windowWidth,
                              const set<string>& regressionMethods,
                              bool centre,
                              int subsetsize,
                              bool extrapolation,
                              int minNumNonNAs,
                              int & errorState)
    {
     errorState = 0;
     responseVector = response;
     fensterbreite = windowWidth;

     if (windowWidth < 3) {
        error("width must be >= 3");
     }


     if (fensterbreite <= 0) {
            //XXX! When we call error, does the destructor get called?
            //No! Thus to avoid memory leaks, all these errors should
            //be converted to throws.
            error("width <= 0");
     }
     if (fensterbreite > responseSize) {
            error("width > responseSize");
     }

     vector<RegLine*> results(numberOfRegressionMethods);
     for (int i = 0; i != numberOfRegressionMethods; ++i) {
         if (regressionMethods.count(regressionMethodNames[i])) {
           results[i] = new RegLine[responseSize];
		   for (int j = 0; j != responseSize; ++j) {
			results[i][j] = RegLine::nullLine;
		   }
         }
         //Otherwise, result[i] already has the default value
         //for pointers, which is 0.
     }

     	// Initialisiere die verschiedenen Algorithmen

     	usingLQD = results[lqdIndex];
	    if (usingLQD)
	    {
	    //	std::cout << "initialisiere lqd" << std::endl;

            const double epsilon = 0.01;
	    	computeLQDAdvanced.init(fensterbreite, subsetsize, epsilon, searchExponent,
                           lqdDetermineStart, lqdApprox, lqdDescendMethod, DBL_EPSILON * 3);
	    }
	    usingHammock = results[rmIndex] || results[lmsIndex] || results[ltsIndex] || results[drIndex];
	    if (usingHammock)
	    {
	    	H.init(fensterbreite);
	    	H.adjust_H_wLen(fensterbreite,subsetsize);
	    }

	    if (results[medIndex]) {
	    	medianFilter = new MedianFilter(fensterbreite);
        }

	    nonNaTracker = new CircularArray<bool>(fensterbreite);

	    //Insert fensterbreite values.

	    for (int i = 0; i != fensterbreite && errorState == 0; ++i) {
            insertNext(errorState);
        }

        int row;
        if (centre) {
           row = fensterbreite / 2;
        }
        else {
            row = fensterbreite - 1;
        }

//        std::cout << "erster record" << std::endl;

        if (numNonNAs >= minNumNonNAs && errorState == 0) {
           recordEstimate(results, row);
        } else {
        	recordNAs(results,row);
        }
        ++row;

        //Process the remaining responseSize - fensterbreite values,
        //recording estimates for each.

        const int remainder = responseSize - fensterbreite;
        //std::cout << "begin der for schleife" << remainder << std::endl;

        for (int j = 0; j != remainder; ++j, ++row) {
        	//std::cout<< "schleife: " << j << "\n";
        	//H.debugInfo =  j > 2149 && j < 2154;
        	//H.debugInfo = j < 120;
        	if (errorState == 0){
				bool oldestWasNotNA = false;
				if (nonNaTracker->size() != 0) {
					oldestWasNotNA = nonNaTracker->oldest();
				}
				const bool wasNotFull = numNonNAs != fensterbreite;


				const bool insertedNA = !insertNext(errorState);
				if (errorState == 0) {
					//Whenever we insert an non-NA into the window, we increment
					//this count. Whenever we remove such a value, we decrement
					//the count.
					if (oldestWasNotNA) {
					   --numNonNAs;
					}

					//We "represent" NAs simply by not inserting a value.
					//This means that we must manually delete a non-NA value
					//when it is pushed off the end of a window containing NAs.
					if (oldestWasNotNA && (insertedNA || wasNotFull)) {
						removeOldest();
					}

					int newH = (subsetsize*(numNonNAs))/(fensterbreite);
					if (usingHammock)
					{
						H.adjust_H_wLen(numNonNAs,newH); // Aenderung der internen Variable Fensterbreite von Hammock duerfte keine auswirkungen auf den algorithmus haben
					}
					if (usingLQD)
					{
						computeLQDAdvanced.adjust_H_wLen(numNonNAs,newH);
					}
					// if the number of elements unequal to NA fall below minNumNonNAs
					// an NA will be recorded instead of the current estimates.
					if (numNonNAs >= minNumNonNAs) {
					   recordEstimate(results, row);
					} else {
						recordNAs(results,row);
					}
				} else {
					recordNAs(results,row);
				}
        	} else {
        		recordNAs(results,row);
        	}
        }

        if (extrapolation) {
            //Fill in elements for which an estimate has not been provided
            //using the nearest element that has.
            //
            //The first and last indexes with an estimate. Initialize them
            //to the values they would have for centre == false.
            int first = fensterbreite - 1;
            int last = responseSize - 1;
            if (centre) {
               first = fensterbreite / 2;
               last = first + responseSize - fensterbreite;
            }

            for (int i = 0; i != numberOfRegressionMethods; ++i) {
                if (results[i]) {
                   RegLine firstValue = results[i][first];
                   for (int j = first - 1; j >= 0; --j) {
                       results[i][j] = firstValue;
                   }

                   RegLine lastValue = results[i][last];
                   for (int j = last + 1; j < responseSize; ++j) {
                       results[i][j] = lastValue;
                   }
                }
            }
        }

        delete nonNaTracker;
        return results;
    }

    /**
     * Return an object containing estimates of slope and level for
     * each method named in the vector regressionMethods.
     * centreExp ^= !online
     */
	SEXP robustRegression(SEXP response, SEXP windowWidth, SEXP regressionMethods,
                               SEXP centreExp, SEXP h, SEXP extrapolation, SEXP minNumNonNAsExp)
	{

     if (!isVector(response) && !isVectorizable(response) && !isFrame(response)) {
        error("response is neither a data.frame nor vectorizable!");
     }
     //XXX! Why not just always require a double vector, forcing the caller
     //to take care of the conversions in R? Otherwise, should we warn the
     //caller when the matrix or data frame has more than one column?

	// Daten aus response einlesen
	    SEXP column = 0;
		if (isFrame(response)) {
			PROTECT(column = coerceVector(VECTOR_ELT(response, 0), REALSXP));
		}
		else {
             PROTECT(column = coerceVector(response, REALSXP));
        }

	    int responseSize;
        if (isMatrix(response)) {
			responseSize = INTEGER(getAttrib(response,R_DimSymbol))[0];
		}
		else {
             responseSize = length(column);
        }

		double* responseVector = new double[responseSize];
        for (int i = 0; i != responseSize; ++i) {
            responseVector[i] = REAL(column)[i];
        }
        UNPROTECT(1);

        //The R expression c() evaluates to something with class NULL,
        //not a vector.
        if (!isVector(regressionMethods) && !isNull(regressionMethods)) {
          error("Did not get vector for regressionMethods");
        }
        const int numMethods = length(regressionMethods);
        set<string> methods;
        for (int i = 0; i < numMethods; ++i) {
            methods.insert(CHAR(STRING_ELT(regressionMethods, i)));
        }
        int errorState = 0;
        vector<RegLine*> results = robustRegression(responseVector,
                                                responseSize,
                                                INTEGER(windowWidth)[0],
                                                methods,
                                                LOGICAL(centreExp)[0],
                                                INTEGER(h)[0],
                                                LOGICAL(extrapolation)[0],
                                                INTEGER(minNumNonNAsExp)[0],
                                                errorState);
		int numEstimates = 0;
		for (int i = 0; i != numberOfRegressionMethods; ++i) {
			if (results[i]) {
				++numEstimates;
			}
		}

		if (numEstimates != numMethods) {
			warning("Got back fewer estimates than expected.");
		}




		SEXP answer = 0;
		PROTECT(answer = allocVector(VECSXP, numEstimates));

		SEXP names = 0;
		PROTECT(names = allocVector(STRSXP, numEstimates));
		int answerIndex = 0;
		for (int i = 0; i != numberOfRegressionMethods; ++i) {
			if (results[i]) {
				SET_STRING_ELT(names, answerIndex, mkChar(regressionMethodNames[i]));

				SEXP estimate = 0;
				PROTECT(estimate = allocVector(VECSXP, 2));
				SEXP estimateNames = 0;
				PROTECT(estimateNames = allocVector(STRSXP, 2));
				SET_STRING_ELT(estimateNames, 0, mkChar("level"));
				SET_STRING_ELT(estimateNames, 1, mkChar("slope"));
				namesgets(estimate, estimateNames);
				UNPROTECT(1);

				SEXP level = 0;
				PROTECT(level = allocVector(REALSXP, responseSize));
				SEXP slope = 0;
				PROTECT(slope = allocVector(REALSXP, responseSize));
				for (int j = 0; j != responseSize; ++j) {
                    const double beta = results[i][j].steigung;
					REAL(slope)[j] = beta;
					//Convert from an estimate of intercept to an estimate of level.
					REAL(level)[j] = ISNA(beta) ? NA_REAL :
                                                  (results[i][j].y_achse + beta * (j + 1));
				}
				SET_VECTOR_ELT(estimate, 0, level);
				SET_VECTOR_ELT(estimate, 1, slope);
				UNPROTECT(2);

				SET_VECTOR_ELT(answer, answerIndex, estimate);
				UNPROTECT(1);

				++answerIndex;
			}
		}
		namesgets(answer, names);

		SEXP mainAnswer = 0;
		SEXP mainNames  = 0;
		PROTECT(mainAnswer = allocVector(VECSXP, 2));
		PROTECT(mainNames  = allocVector(STRSXP, 2));
		SET_STRING_ELT(mainNames, 0, mkChar("Error_State"));
		SET_STRING_ELT(mainNames, 1, mkChar("Estimation_Result"));
		SEXP eState = 0;

		PROTECT(eState = allocVector(INTSXP, 1));
		SET_VECTOR_ELT(mainAnswer, 0, eState);
		SET_VECTOR_ELT(mainAnswer, 1, answer);

		INTEGER(eState)[0] = errorState;
		namesgets(mainAnswer, mainNames);
		UNPROTECT(3);
		UNPROTECT(2);
 

		delete[] responseVector;
		for (int i = 0; i < numberOfRegressionMethods; ++i) {
            delete[] results[i];
        }

		return mainAnswer;
}

};

extern "C"
//By default, DevC++ defines the BUILDING_DLL macro
//when building a DLL
#ifdef BUILDING_DLL
//Used to mark a symbol for export in DevC++.
__declspec(dllexport)
#endif
	SEXP robustRegression(SEXP response, SEXP windowWidth, SEXP regressionMethods,
                               SEXP centreExp, SEXP h, SEXP epsilonI, SEXP minNumNonNAs)
{
   SEXP result = 0;
   try {
    RobustReg r;
    result = r.robustRegression(response, windowWidth, regressionMethods,
                                       centreExp,
                                       h,
                                       epsilonI,
                                       minNumNonNAs);
  }
  catch (std::exception e) {
    error(e.what());
  }
  return result;
}
#endif


