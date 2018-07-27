#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <stdio>

#include "Common.h"
#include "MathUtils.h"

MathUtils::MathUtils() {
  log10QempPriorCache = new double[MAX_GATK_USABLE_Q_SCORE + 1];

  const double GF_a = 0.0;
  const double GF_b = 0.9;
  const double GF_c = 0.0;
  const double GF_d = 0.5;

  for (int i = 0; i <= MAX_GATK_USABLE_Q_SCORE; i++) {
    double iMc = (double)(i) - GF_c;
    double value = GF_a + GF_b*exp(-1*iMc * iMc / (2.0 * (GF_d * GF_d)));
    double log10Prior = log10(value);
    if (std::isinf(log10Prior)) {
      log10Prior = -(std::numeric_limits<double>::max());
    }
    log10QempPriorCache[i] = log10Prior;
  }

  Log10FactorialCache = new double[Log10FactorialCacheSize];
  Log10FactorialCache[0] = 0.0;
  for (int i = 1; i < Log10FactorialCacheSize; i++) {
    Log10FactorialCache[i] = Log10FactorialCache[i-1] + std::log10(i);
  }
}

MathUtils::~MathUtils() {
  delete [] log10QempPriorCache;
  delete [] Log10FactorialCache;
}

double MathUtils::bayesianEstimateOfEmpiricalQuality(
    uint64_t nObservations,
    uint64_t nErrors,
    double QReported)
{
  int numBins = (MAX_REASONABLE_Q_SCORE + 1) * (int)RESOLUTION_BINS_PER_QUAL;

  double* log10Posteriors = (double*)malloc(numBins*sizeof(double));

  // NOTE: based on timer setup, this for loop takes most of the time.
  double max_val = -(std::numeric_limits<double>::max());
  for (int bin = 0; bin < numBins; bin++) {
    double QEmpOfBin = (double)bin / RESOLUTION_BINS_PER_QUAL;
    log10Posteriors[bin] = log10QempPrior(QEmpOfBin, QReported) +
                           log10QempLikelihood(QEmpOfBin, nObservations, nErrors);
    if (max_val < log10Posteriors[bin]) {
      max_val = log10Posteriors[bin];
    }
  }

  // double[] normalizedPosteriors = MathUtils.normalizeFromLog10(log10Posteriors);
  // int MLEbin = MathUtils.maxElementIndex(normalizedPosteriors);
  double sum = 0.0;
  for (int bin = 0; bin < numBins; bin++) {
    log10Posteriors[bin] = std::pow(10, log10Posteriors[bin] - max_val);
    sum += log10Posteriors[bin];
  }
  int MLEbin = 0;
  for (int bin = 0; bin < numBins; bin++) {
    log10Posteriors[bin] = log10Posteriors[bin] / sum;
    if (log10Posteriors[bin] > log10Posteriors[MLEbin]) {
      MLEbin = bin;
    }
  }
  free(log10Posteriors);

  return MLEbin / RESOLUTION_BINS_PER_QUAL;
}


double MathUtils::log10QempPrior(double Qempirical, double Qreported) {
  int difference = std::min(std::abs((int)(Qempirical - Qreported)),
                            (int)MAX_GATK_USABLE_Q_SCORE);
  fprintf(stderr, difference);
  return log10QempPriorCache[difference];
}

static inline double qualToErrorProbLog10(double qual) {
  return qual * -0.1;
}

static inline int HI(double x) {
  uint64_t ret;
  // copy the bits of x to result
  memcpy(&ret, &x, sizeof(x));
  return (int)(ret >> 32);
}

static inline int LO(double x) {
  uint64_t ret;
  // copy the bits of x to result
  memcpy(&ret, &x, sizeof(x));
  return (int)(ret);
}

static double zero = 0.0, one = 1.0, half = .5, a0 = 7.72156649015328655494e-02, a1 = 3.22467033424113591611e-01, a2 = 6.73523010531292681824e-02, a3 = 2.05808084325167332806e-02, a4 = 7.38555086081402883957e-03, a5 = 2.89051383673415629091e-03, a6 = 1.19270763183362067845e-03, a7 = 5.10069792153511336608e-04, a8 = 2.20862790713908385557e-04, a9 = 1.08011567247583939954e-04, a10 = 2.52144565451257326939e-05, a11 = 4.48640949618915160150e-05, tc = 1.46163214496836224576e+00, tf = -1.21486290535849611461e-01, tt = -3.63867699703950536541e-18, t0 = 4.83836122723810047042e-01, t1 = -1.47587722994593911752e-01, t2 = 6.46249402391333854778e-02, t3 = -3.27885410759859649565e-02, t4 = 1.79706750811820387126e-02, t5 = -1.03142241298341437450e-02, t6 = 6.10053870246291332635e-03, t7 = -3.68452016781138256760e-03, t8 = 2.25964780900612472250e-03, t9 = -1.40346469989232843813e-03, t10 = 8.81081882437654011382e-04, t11 = -5.38595305356740546715e-04, t12 = 3.15632070903625950361e-04, t13 = -3.12754168375120860518e-04, t14 = 3.35529192635519073543e-04, u0 = -7.72156649015328655494e-02, u1 = 6.32827064025093366517e-01, u2 = 1.45492250137234768737e+00, u3 = 9.77717527963372745603e-01, u4 = 2.28963728064692451092e-01, u5 = 1.33810918536787660377e-02, v1 = 2.45597793713041134822e+00, v2 = 2.12848976379893395361e+00, v3 = 7.69285150456672783825e-01, v4 = 1.04222645593369134254e-01, v5 = 3.21709242282423911810e-03, s0 = -7.72156649015328655494e-02, s1 = 2.14982415960608852501e-01, s2 = 3.25778796408930981787e-01, s3 = 1.46350472652464452805e-01, s4 = 2.66422703033638609560e-02, s5 = 1.84028451407337715652e-03, s6 = 3.19475326584100867617e-05, r1 = 1.39200533467621045958e+00, r2 = 7.21935547567138069525e-01, r3 = 1.71933865632803078993e-01, r4 = 1.86459191715652901344e-02, r5 = 7.77942496381893596434e-04, r6 = 7.32668430744625636189e-06, w0 = 4.18938533204672725052e-01, w1 = 8.33333333333329678849e-02, w2 = -2.77777777728775536470e-03, w3 = 7.93650558643019558500e-04, w4 = -5.95187557450339963135e-04, w5 = 8.36339918996282139126e-04, w6 = -1.63092934096575273989e-03;

static inline double lnGamma(double x) {
  double t, y, z, p, p1, p2, p3, q, r, w;
  int i;

  int hx = HI(x);
  int lx = LO(x);

  /* purge off +-inf, NaN, +-0, and negative arguments */
  int ix = hx & 0x7fffffff;
  if (ix >= 0x7ff00000)
    return std::numeric_limits<double>::infinity();

  if ((ix | lx) == 0 || hx < 0)
    return std::numeric_limits<double>::signaling_NaN();

  if (ix < 0x3b900000) {    /* |x|<2**-70, return -log(|x|) */
    return -std::log(x);
  }

  /* purge off 1 and 2 */
  if ((((ix - 0x3ff00000) | lx) == 0) || (((ix - 0x40000000) | lx) == 0))
  r = 0;
  /* for x < 2.0 */
  else if (ix < 0x40000000) {
    if (ix <= 0x3feccccc) {     /* lgamma(x) = lgamma(x+1)-log(x) */
      r = -std::log(x);
      if (ix >= 0x3FE76944) {
        y = one - x;
        i = 0;
      }
      else if (ix >= 0x3FCDA661) {
        y = x - (tc - one);
        i = 1;
      }
      else {
        y = x;
        i = 2;
      }
    }
    else {
      r = zero;
      if (ix >= 0x3FFBB4C3) {
        y = 2.0 - x;
        i = 0;
      } /* [1.7316,2] */
      else if (ix >= 0x3FF3B4C4) {
        y = x - tc;
        i = 1;
      } /* [1.23,1.73] */
      else {
        y = x - one;
        i = 2;
      }
    }

    switch (i) {
      case 0:
        z = y * y;
        p1 = a0 + z * (a2 + z * (a4 + z * (a6 + z * (a8 + z * a10))));
        p2 = z * (a1 + z * (a3 + z * (a5 + z * (a7 + z * (a9 + z * a11)))));
        p = y * p1 + p2;
        r += (p - 0.5 * y);
        break;
      case 1:
        z = y * y;
        w = z * y;
        p1 = t0 + w * (t3 + w * (t6 + w * (t9 + w * t12)));    /* parallel comp */
        p2 = t1 + w * (t4 + w * (t7 + w * (t10 + w * t13)));
        p3 = t2 + w * (t5 + w * (t8 + w * (t11 + w * t14)));
        p = z * p1 - (tt - w * (p2 + y * p3));
        r += (tf + p);
        break;
      case 2:
        p1 = y * (u0 + y * (u1 + y * (u2 + y * (u3 + y * (u4 + y * u5)))));
        p2 = one + y * (v1 + y * (v2 + y * (v3 + y * (v4 + y * v5))));
        r += (-0.5 * y + p1 / p2);
        break;
      default: ;
    }
  }
  else if (ix < 0x40200000) {             /* x < 8.0 */
    i = (int) x;
    t = zero;
    y = x - (double) i;
    p = y * (s0 + y * (s1 + y * (s2 + y * (s3 + y * (s4 + y * (s5 + y * s6))))));
    q = one + y * (r1 + y * (r2 + y * (r3 + y * (r4 + y * (r5 + y * r6)))));
    r = half * y + p / q;
    z = one;    /* lgamma(1+s) = log(s) + lgamma(s) */
    switch (i) {
      case 7:
        z *= (y + 6.0);    /* FALLTHRU */
      case 6:
        z *= (y + 5.0);    /* FALLTHRU */
      case 5:
        z *= (y + 4.0);    /* FALLTHRU */
      case 4:
        z *= (y + 3.0);    /* FALLTHRU */
      case 3:
        z *= (y + 2.0);    /* FALLTHRU */
        r += std::log(z);
      break;
    }
    /* 8.0 <= x < 2**58 */
  }
  else if (ix < 0x43900000) {
    t = std::log(x);
    z = one / x;
    y = z * z;
    w = w0 + z * (w1 + y * (w2 + y * (w3 + y * (w4 + y * (w5 + y * w6)))));
    r = (x - half) * (t - one) + w;
  }
  else {
    /* 2**58 <= x <= inf */
    r = x * (std::log(x) - one);
  }
  return r;
}

static inline double log10Gamma(double x) {
  return lnGamma(x) * std::log10(M_E);
}

inline double MathUtils::log10Factorial(int x) {
  if (x >= Log10FactorialCacheSize || x < 0) {
    return log10Gamma(x+1);
  }
  else {
    return Log10FactorialCache[x];
  }
}

inline double MathUtils::log10BinomialCoefficient(int n, int k) {
  if (n < 0) {
    throw new std::runtime_error("n: Must have non-negative number of trials");
  }
  if (k > n || k < 0) {
    throw new std::runtime_error("k: Must have non-negative number of successes, and no more successes than number of trials");
  }
  return log10Factorial(n) - log10Factorial(k) - log10Factorial(n - k);
}

double MathUtils::log10BinomialProbability(int n, int k, double log10p) {
  if (log10p > 1e-18) {
    throw new std::runtime_error("log10p: Log-probability must be 0 or less");
  }
  double log10OneMinusP = std::log10(1 - std::pow(10, log10p));
  return log10BinomialCoefficient(n, k) + log10p * k + log10OneMinusP * (n - k);
}

double MathUtils::log10QempLikelihood(double Qempirical,
                            uint64_t nObservations,
                            uint64_t nErrors)
{
  if (nObservations == 0) return 0.0;

  // the binomial code requires ints as input (because it does caching).  This should theoretically be fine because
  // there is plenty of precision in 2^31 observations, but we need to make sure that we don't have overflow
  // before casting down to an int.
  uint64_t MAX_NUMBER_OF_OBSERVATIONS = std::numeric_limits<int>::max() - 1;
  if (nObservations > MAX_NUMBER_OF_OBSERVATIONS) {
    // we need to decrease nErrors by the same fraction that we are decreasing nObservations
    double fraction = (double)MAX_NUMBER_OF_OBSERVATIONS / (double)nObservations;
    nErrors = std::round((double)nErrors * fraction);
    nObservations = MAX_NUMBER_OF_OBSERVATIONS;
  }

  // this is just a straight binomial PDF
  double log10Prob = log10BinomialProbability((int)nObservations, (int)nErrors,
                        qualToErrorProbLog10(Qempirical));
  if (std::isinf(log10Prob) || std::isnan(log10Prob)) {
    log10Prob = -(std::numeric_limits<double>::max());
  }
  return log10Prob;
}
