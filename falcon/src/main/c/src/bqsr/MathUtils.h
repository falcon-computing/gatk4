#ifndef BQSR_MATH_UTILS_H
#define BQSR_MATH_UTILS_H

#include <limits>

#include "SimpleTimer.h"

// Implementation of selected routines in
// GATK MathUtils and QualityUtils used for BQSR and PR
class MathUtils {
  public:
    // initialize caches
    MathUtils();

    ~MathUtils();

    double bayesianEstimateOfEmpiricalQuality(
          uint64_t nObservations,
          uint64_t nErrors,
          double QReported);

    double log10QempPrior(double Qempirical, double Qreported);
    double log10QempLikelihood(double Qempirical,
              uint64_t nObservations, uint64_t nErrors);

  protected:
    inline double log10BinomialCoefficient(int n, int k);
    inline double log10BinomialProbability(int n, int k, double log10p);
    inline double log10Factorial(int x);

    const int MAX_PHRED_SCORE = 93;
    const double RESOLUTION_BINS_PER_QUAL = 1.0;
    const int8_t MAX_REASONABLE_Q_SCORE = 60;
    const int8_t MAX_GATK_USABLE_Q_SCORE = 40;

    double* log10QempPriorCache;
    double* Log10FactorialCache;
    const int Log10FactorialCacheSize = 1024*1024;

    SimpleTimer timer_;
};
#endif
