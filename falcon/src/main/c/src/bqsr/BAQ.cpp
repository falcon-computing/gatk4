#include <iostream>
#include <iterator>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>
#include <stdexcept>

#include "Common.h"
#include "BAQ.h"

BAQ::BAQ(): cd(-1), ce(0.1), cb(7),
    MAX_PHRED_SCORE(93),
    EM(0.33333333333),
    EI(0.25),
    includeClippedBases(false),
    minBaseQual(4),
    hmm_time_ns(0),
    baq_time_ns(0),
    fracerror_time_ns(0)
{
  // NOTE: Here initialize BAQ parameters using the default value
  // TODO: need to support configurable parameters in the future
  const int DEFAULT_GOP = 40;
  cd = std::pow(10, (-DEFAULT_GOP)/10.);

  double qual2prob[256];
  for (int i = 0; i < 256; ++i) {
    qual2prob[i] = std::pow(10, (double)-i/10.);
  }

  for (int i = 0; i < 256; i++) {
    for (int j = 0; j < 256; j++) {
      for (int q = 0; q <= MAX_PHRED_SCORE; q++ ) {
        EPSILONS[i][j][q] = 1.0;
      }
    }
  }

  const char bases[] = {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'};
  for (int i = 0; i < 8; i++) {
    char b1 = bases[i];
    for (int j = 0; j < 8; j++) {
      char b2 = bases[j];
      for (int q = 0; q <= MAX_PHRED_SCORE; q++) {
        double qual = qual2prob[q < minBaseQual ? minBaseQual : q];
        double e = std::tolower(b1) == std::tolower(b2) ? 1 - qual : qual * EM;
        EPSILONS[(int)b1][(int)b2][q] = e;
      }
    }
  }
}

BAQ::~BAQ() {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Total time to calculate hmm: " << (double)hmm_time_ns / 1e6 << " ms";
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Total time to calculate baq: " << (double)baq_time_ns / 1e6 << " ms";
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Total time to calculate errors: " << (double)fracerror_time_ns / 1e6 << " ms";
}

// helper routines for BAQ calculation
static inline bool stateIsIndel(int state) {
  return (state & 3) != 0;
}

/** decode the bit encoded state array values */
static inline int stateAlignedPosition(int state) {
  return state >> 2;
}

inline int8_t BAQ::capBaseByBAQ(int8_t oq,
                                int8_t bq,
                                int state,
                                int expectedPos)
{
  int8_t b;
  bool isIndel = stateIsIndel(state);
  int pos = stateAlignedPosition(state);
  if (isIndel || pos != expectedPos) // we are an indel or we don't align to our best current position
    b = minBaseQual; // just take b = minBaseQuality
  else
    b = bq < oq ? bq : oq;

  return b;
}

int BAQ::calculateBAQ(int8_t* bqTag,
        int8_t* bases,
        int8_t* basesQuals,
        int8_t* refBases,
        int8_t* cigarOps,
        int* cigarLens,
        int refLength,
        int readLength,
        int numCigar,
        int refOffset)
{
  // calculate query range (BAQ.calculateQueryRange())
  int queryStart = -1;
  int queryStop = -1;
  int readI = 0;
  // iterate over the cigar elements to determine the start and stop of the read bases for the BAQ calculation
  //for (CigarElement elt : cigarElements) {
  for (int i = 0; i < numCigar; i++) {
    int prev;
    switch (cigarOps[i]) {
      case 'N':
        return 1; // cannot handle these
      case 'H':
      case 'P':
      case 'D':
        break; // ignore pads, hard clips, and deletions
      case 'I':
      case 'S':
      case 'M':
      case '=': //case EQ:
      case 'X':
        prev = readI;
        readI += cigarLens[i];
        if ( includeClippedBases || cigarOps[i] != 'S') {
          if ( queryStart == -1 )
            queryStart = prev;
          queryStop = readI;
        }
        // in the else case we aren't including soft clipped bases, so we don't update
        // queryStart or queryStop
        break;
      default:
        throw std::runtime_error("BUG: Unexpected CIGAR element in read");
    }
  }

  if (queryStop == queryStart) {
    // this read is completely clipped away, and yet is present in the file for some reason
    // usually they are flagged as non-PF, but it's possible to push them through the BAM
    //System.err.printf("WARNING -- read is completely clipped away: " + read.format());
    return 1;
  }

  //BAQ.BAQCalculationResult baqResult = new BAQ.BAQCalculationResult(query, quals, ref);
  int queryLen = queryStop - queryStart;
  int* baqState = (int*)calloc(readLength, sizeof(int));
  int8_t* baqBq = (int8_t*)calloc(readLength, sizeof(int));
  int error = 0;

  //DLOG(INFO) << "refLength = " << refLength;
  //DLOG(INFO) << "queryStart = " << queryStart;
  //DLOG(INFO) << "queryLen = " << queryLen;
  //DLOG(INFO) << "basesQuals = ";

  //std::cerr << "refBases = ";
  //std::copy(refBases, refBases + refLength,
  //    std::ostream_iterator<char>(std::cerr, " "));
  //std::cerr << std::endl;
//
  //std::cerr << "bases = ";
  //std::copy(bases, bases + readLength,;
  //    std::ostream_iterator<char>(std::cerr, " "));
  //std::cerr << std::endl;
//
  //std::cerr << "basesQuals = ";
  //std::copy(basesQuals, basesQuals + readLength,
  //    std::ostream_iterator<int>(std::cerr, " "));
  //std::cerr << std::endl;

  uint64_t start_ns = getNs();
  hmm_glocal(refLength,
      refBases, bases,
      queryStart, queryLen,
      basesQuals, baqState, baqBq);
  hmm_time_ns += getNs() - start_ns;

  //std::cerr << "baqBq = ";
  //std::copy(baqBq + queryStart, baqBq + queryStop,
  //    std::ostream_iterator<int>(std::cerr, " "));
  //std::cerr << std::endl;

  readI = 0;
  int refI = 0;
  //for (CigarElement elt : cigarElements) {
  for (int k = 0; k < numCigar; k++) {
    int l = cigarLens[k];
    int expectedPos;
    switch (cigarOps[k]) {
      case 'N': // cannot handle these
        error = 1;
        goto func_return;
      case 'H':
      case 'P': // ignore pads and hard clips
        break;
      case 'S':
        refI += l; // move the reference too, in addition to I
      case 'I':
        // todo -- is it really the case that we want to treat I and S the same?
        for (int i = readI; i < readI + l; i++) baqBq[i] = basesQuals[i];
        readI += l;
        break;
      case 'D':
        refI += l;
        break;
      case 'M':
      case '=':
      case 'X':
        for (int i = readI; i < readI + l; i++) {
          expectedPos = refI - refOffset + (i - readI);
          baqBq[i] = capBaseByBAQ(basesQuals[i], baqBq[i], baqState[i], expectedPos);
        }
        readI += l; refI += l;
        break;
      default:
        throw std::runtime_error("BUG: Unexpected CIGAR element in read");
    }
  }
  if (readI != readLength) {// odd cigar string
    DLOG(INFO) << "Odd cigar string";
    memcpy(baqBq, basesQuals, readLength);
  }

  //BAQ.encodeBQTag(read, hmmResult.bq);
  for (int i = 0; i < readLength; i++) {
    int bq = (int)basesQuals[i] + 64;
    int baq_i = (int)baqBq[i];
    int tag = bq - baq_i;
    // problem with the calculation of the correction factor; this is our problem
    if (tag < 0) {
      throw std::runtime_error(
        "BAQ tag calculation error. BAQ value above base quality");
    }
    // the original quality is too high, almost certainly due to using the wrong encoding in the BAM file
    if (tag > SCHAR_MAX) {
      throw std::runtime_error(
        "we encountered an extremely high quality score with "
        "BAQ correction factor");
    }
    bqTag[i] = (int8_t)tag;
  }

func_return:
  free(baqState);
  free(baqBq);
  return error;
}

/**
 * helper routine for hmm_glocal
 *
 * @param b
 * @param i
 * @param k
 * @return
 */
static inline int set_u(int b, int i, int k) {
  int x = i - b;
  x = x > 0 ? x : 0;
  return (k + 1 - x) * 3;
}

int BAQ::hmm_glocal(int ref_len,
          int8_t* ref,
          int8_t* query,
          int qstart,
          int l_query,
          int8_t* _iqual,
          int* state,
          int8_t* q)
{
  int i, k;

  /*** initialization ***/
  // change coordinates
  int l_ref = ref_len;

  // set band width
  int bw2, bw = l_ref > l_query? l_ref : l_query;
  if (cb < std::abs(l_ref - l_query)) {
    bw = std::abs(l_ref - l_query) + 3;
    //System.out.printf("SC  cb=%d, bw=%d%n", cb, bw);
  }
  if (bw > cb) bw = cb;
  if (bw < std::abs(l_ref - l_query)) {
    //int bwOld = bw;
    bw = std::abs(l_ref - l_query);
    //System.out.printf("old bw is %d, new is %d%n", bwOld, bw);
  }
  //System.out.printf("c->bw = %d, bw = %d, l_ref = %d, l_query = %d\n", cb, bw, l_ref, l_query);
  bw2 = bw * 2 + 1;

  // allocate the forward and backward matrices f[][] and b[][] and the scaling array s[]
  double** f = (double**)malloc((l_query+1)*sizeof(double*));
  double** b = (double**)malloc((l_query+1)*sizeof(double*));
  for (int i = 0; i <= l_query; i++) {
    f[i] = (double*)calloc(bw2*3 + 6, sizeof(double));
    b[i] = (double*)calloc(bw2*3 + 6, sizeof(double));
  }
  double*  s = (double*)calloc(l_query+2, sizeof(double));

  // initialize transition probabilities
  double sM, sI, bM, bI;
  sM = sI = 1. / (2 * l_query + 2);
  bM = (1 - cd) / l_ref; bI = cd / l_ref; // (bM+bI)*l_ref==1

  double m[9];
  m[0*3+0] = (1 - cd - cd) * (1 - sM); m[0*3+1] = m[0*3+2] = cd * (1 - sM);
  m[1*3+0] = (1 - ce) * (1 - sI); m[1*3+1] = ce * (1 - sI); m[1*3+2] = 0.;
  m[2*3+0] = 1 - ce; m[2*3+1] = 0.; m[2*3+2] = ce;

  /*** forward ***/
  // f[0]
  f[0][set_u(bw, 0, 0)] = s[0] = 1.;
  { // f[1]
    double* fi = f[1];
    double sum;
    int beg = 1, end = l_ref < bw + 1? l_ref : bw + 1, _beg, _end;
    for (k = beg, sum = 0.; k <= end; ++k) {
      int u;
      double e = calcEpsilon(ref[k-1], query[qstart], _iqual[qstart]);
      u = set_u(bw, 1, k);
      fi[u+0] = e * bM; fi[u+1] = EI * bI;
      sum += fi[u] + fi[u+1];
    }
    // rescale
    s[1] = sum;
    _beg = set_u(bw, 1, beg); _end = set_u(bw, 1, end); _end += 2;
    for (k = _beg; k <= _end; ++k) fi[k] /= sum;
  }

  // f[2..l_query]
  for (i = 2; i <= l_query; ++i) {
    double* fi = f[i];
    double* fi1 = f[i-1];
    double sum;
    int beg = 1, end = l_ref, x, _beg, _end;
    int8_t qyi = query[qstart+i-1];
    x = i - bw; beg = beg > x? beg : x; // band start
    x = i + bw; end = end < x? end : x; // band end
    for (k = beg, sum = 0.; k <= end; ++k) {
      int u, v11, v01, v10;
      double e = calcEpsilon(ref[k-1], qyi, _iqual[qstart+i-1]);
      u = set_u(bw, i, k); v11 = set_u(bw, i-1, k-1); v10 = set_u(bw, i-1, k); v01 = set_u(bw, i, k-1);
      fi[u+0] = e * (m[0] * fi1[v11+0] + m[3] * fi1[v11+1] + m[6] * fi1[v11+2]);
      fi[u+1] = EI * (m[1] * fi1[v10+0] + m[4] * fi1[v10+1]);
      fi[u+2] = m[2] * fi[v01+0] + m[8] * fi[v01+2];
      sum += fi[u] + fi[u+1] + fi[u+2];
      //System.out.println("("+i+","+k+";"+u+"): "+fi[u]+","+fi[u+1]+","+fi[u+2]);
    }
    // rescale
    s[i] = sum;
    _beg = set_u(bw, i, beg); _end = set_u(bw, i, end); _end += 2;
    for (k = _beg, sum = 1./sum; k <= _end; ++k) fi[k] *= sum;
  }
  { // f[l_query+1]
    double sum;
    for (k = 1, sum = 0.; k <= l_ref; ++k) {
      int u = set_u(bw, l_query, k);
      if (u < 3 || u >= bw2*3+3) continue;
      sum += f[l_query][u+0] * sM + f[l_query][u+1] * sI;
    }
    s[l_query+1] = sum; // the last scaling factor
  }

  /*** backward ***/
  // b[l_query] (b[l_query+1][0]=1 and thus \tilde{b}[][]=1/s[l_query+1]; this is where s[l_query+1] comes from)
  for (k = 1; k <= l_ref; ++k) {
    int u = set_u(bw, l_query, k);
    double* bi = b[l_query];
    if (u < 3 || u >= bw2*3+3) continue;
    bi[u+0] = sM / s[l_query] / s[l_query+1]; bi[u+1] = sI / s[l_query] / s[l_query+1];
  }
  // b[l_query-1..1]
  for (i = l_query - 1; i >= 1; --i) {
    int beg = 1, end = l_ref, x, _beg, _end;
    double* bi = b[i];
    double* bi1 = b[i+1];
    double y = (i > 1)? 1. : 0.;
    int8_t qyi1 = query[qstart+i];
    x = i - bw; beg = beg > x? beg : x;
    x = i + bw; end = end < x? end : x;
    for (k = end; k >= beg; --k) {
      int u, v11, v01, v10;
      u = set_u(bw, i, k); v11 = set_u(bw, i+1, k+1); v10 = set_u(bw, i+1, k); v01 = set_u(bw, i, k+1);
      double e = (k >= l_ref? 0 : calcEpsilon(ref[k], qyi1, _iqual[qstart+i])) * bi1[v11];
      bi[u+0] = e * m[0] + EI * m[1] * bi1[v10+1] + m[2] * bi[v01+2]; // bi1[v11] has been folded into e.
      bi[u+1] = e * m[3] + EI * m[4] * bi1[v10+1];
      bi[u+2] = (e * m[6] + m[8] * bi[v01+2]) * y;
    }
    // rescale
    _beg = set_u(bw, i, beg); _end = set_u(bw, i, end); _end += 2;
    for (k = _beg, y = 1./s[i]; k <= _end; ++k) bi[k] *= y;
  }

  //double pb;
  { // b[0]
    int beg = 1, end = l_ref < bw + 1? l_ref : bw + 1;
    double sum = 0.;
    for (k = end; k >= beg; --k) {
      int u = set_u(bw, 1, k);
      double e = calcEpsilon(ref[k-1], query[qstart], _iqual[qstart]);
      if (u < 3 || u >= bw2*3+3) continue;
      sum += e * b[1][u+0] * bM + EI * b[1][u+1] * bI;
    }
    // pb is used for debugging only
    //pb = b[0][set_u(bw, 0, 0)] = sum / s[0]; // if everything works as is expected, pb == 1.0
  }


  /*** MAP ***/
  for (i = 1; i <= l_query; ++i) {
    double sum = 0., max = 0.;
    double* fi = f[i];
    double* bi = b[i];
    int beg = 1, end = l_ref, x, max_k = -1;
    x = i - bw; beg = beg > x? beg : x;
    x = i + bw; end = end < x? end : x;
    for (k = beg; k <= end; ++k) {
      int u = set_u(bw, i, k);
      double z;
      sum += (z = fi[u+0] * bi[u+0]); if (z > max) { max = z; max_k = (k-1)<<2 | 0; }
      sum += (z = fi[u+1] * bi[u+1]); if (z > max) { max = z; max_k = (k-1)<<2 | 1; }
    }
    max /= sum; sum *= s[i]; // if everything works as is expected, sum == 1.0
    if (state != NULL) state[qstart+i-1] = max_k;
    if (q != NULL) {
      k = (int)(-4.343 * std::log(1. - max) + .499); // = 10*log10(1-max)
      q[qstart+i-1] = (int8_t)(k > 100? 99 : (k < minBaseQual ? minBaseQual : k));
    }
    //System.out.println("("+pb+","+sum+")"+" ("+(i-1)+","+(max_k>>2)+","+(max_k&3)+","+max+")");
  }
  for (int i = 0; i <= l_query; i++) {
    free(f[i]);
    free(b[i]);
  }
  free(f);
  free(b);
  free(s);
  return 0;
}


static inline void calcErrorBlockSkipIndel(
    int bound, int blockStartIndex,
    int8_t* snpArray,
    //int8_t* insertArray,
    //int8_t* deleteArray,
    double* snpErrors)
    //double* snpErrors,
    //double* insertErrors,
    //double* deleteErrors)
{
  int totalSnpErrors = 0;
  int totalInsertErrors = 0;
  int totalDeleteErrors = 0;
  for (int i = std::max(0, blockStartIndex - 1); i <= bound; i++) {
    totalSnpErrors    += (int)snpArray[i];
    //totalInsertErrors += (int)insertArray[i];
    //totalDeleteErrors += (int)deleteArray[i];
  }
  for (int i = std::max(0, blockStartIndex - 1); i <= bound; i++) {
    snpErrors[i] = ((double)totalSnpErrors) / ((double)(bound - std::max(0, blockStartIndex-1) + 1));
    //insertErrors[i] = ((double)totalInsertErrors) / ((double)(bound - std::max(0, blockStartIndex-1) + 1));
    //deleteErrors[i] = ((double)totalDeleteErrors) / ((double)(bound - std::max(0, blockStartIndex-1) + 1));
  }
}

static inline void calcErrorBlock(
    int bound, int blockStartIndex,
    int8_t* snpArray,
    int8_t* insertArray,
    int8_t* deleteArray,
    double* snpErrors,
    double* insertErrors,
    double* deleteErrors)
{
  int totalSnpErrors = 0;
  int totalInsertErrors = 0;
  int totalDeleteErrors = 0;
  for (int i = std::max(0, blockStartIndex - 1); i <= bound; i++) {
    totalSnpErrors    += (int)snpArray[i];
    totalInsertErrors += (int)insertArray[i];
    totalDeleteErrors += (int)deleteArray[i];
  }
  for (int i = std::max(0, blockStartIndex - 1); i <= bound; i++) {
    snpErrors[i] = ((double)totalSnpErrors) / ((double)(bound - std::max(0, blockStartIndex-1) + 1));
    insertErrors[i] = ((double)totalInsertErrors) / ((double)(bound - std::max(0, blockStartIndex-1) + 1));
    deleteErrors[i] = ((double)totalDeleteErrors) / ((double)(bound - std::max(0, blockStartIndex-1) + 1));
  }
}

void BAQ::calculateFractionalErrorArraySkipIndel(
        double* snpErrors,
        //double* insertErrors,
        //double* deleteErrors,
        int8_t* snpArray,
        //int8_t* insertArray,
        //int8_t* deleteArray,
        int8_t* baqArray,
        int readLength)
{
  const int BLOCK_START_UNSET = -1;
  const int8_t NO_BAQ_UNCERTAINTY = (int8_t)'@';

  // for some reason memset does not work
  for (int i = 0; i < readLength; i++) {
    snpErrors[i] = 0.;
    //insertErrors[i] = 0.;
    //deleteErrors[i] = 0.;
  }

  bool inBlock = false;
  int blockStartIndex = BLOCK_START_UNSET;
  int iii = 0;
  for( iii = 0; iii < readLength; iii++ ) {
    if( baqArray[iii] == NO_BAQ_UNCERTAINTY ) {
      if (!inBlock) {
        snpErrors[iii]    = (double)snpArray[iii];
        //insertErrors[iii] = (double)insertArray[iii];
        //deleteErrors[iii] = (double)deleteArray[iii];
      }
      else {
        //calcErrorBlock(iii, blockStartIndex,
        //    snpArray, insertArray, deleteArray,
        //    snpErrors, insertErrors, deleteErrors);
        calcErrorBlockSkipIndel(iii, blockStartIndex,
                    snpArray,
                    snpErrors);

        inBlock = false; // reset state variables
        blockStartIndex = BLOCK_START_UNSET; // reset state variables
      }
    }
    else {
      inBlock = true;
      if (blockStartIndex == BLOCK_START_UNSET) blockStartIndex = iii;
    }
  }
  if (inBlock) {
    //calcErrorBlock(iii-1, blockStartIndex,
    //    snpArray, insertArray, deleteArray,
    //    snpErrors, insertErrors, deleteErrors);
    calcErrorBlockSkipIndel(iii-1, blockStartIndex,
            snpArray,
            snpErrors);
  }
}

void BAQ::calculateFractionalErrorArray(
        double* snpErrors,
        double* insertErrors,
        double* deleteErrors,
        int8_t* snpArray,
        int8_t* insertArray,
        int8_t* deleteArray,
        int8_t* baqArray,
        int readLength)
{
  const int BLOCK_START_UNSET = -1;
  const int8_t NO_BAQ_UNCERTAINTY = (int8_t)'@';

  // for some reason memset does not work
  for (int i = 0; i < readLength; i++) {
    snpErrors[i] = 0.;
    insertErrors[i] = 0.;
    deleteErrors[i] = 0.;
  }

  bool inBlock = false;
  int blockStartIndex = BLOCK_START_UNSET;
  int iii = 0;
  for( iii = 0; iii < readLength; iii++ ) {
    if( baqArray[iii] == NO_BAQ_UNCERTAINTY ) {
      if (!inBlock) {
        snpErrors[iii]    = (double)snpArray[iii];
        insertErrors[iii] = (double)insertArray[iii];
        deleteErrors[iii] = (double)deleteArray[iii];
      }
      else {
        calcErrorBlock(iii, blockStartIndex,
            snpArray, insertArray, deleteArray,
            snpErrors, insertErrors, deleteErrors);

        inBlock = false; // reset state variables
        blockStartIndex = BLOCK_START_UNSET; // reset state variables
      }
    }
    else {
      inBlock = true;
      if (blockStartIndex == BLOCK_START_UNSET) blockStartIndex = iii;
    }
  }
  if (inBlock) {
    calcErrorBlock(iii-1, blockStartIndex,
        snpArray, insertArray, deleteArray,
        snpErrors, insertErrors, deleteErrors);
  }
}

static inline int getBaseIdx(int8_t base) {
  switch (base) {
    case 'A':
    case 'a':
    case '*':
      return 0;
    case 'C':
    case 'c':
      return 1;
    case 'G':
    case 'g':
      return 2;
    case 'T':
    case 't':
      return 3;
    default:
      return -1;
  }
}

int BAQ::calculateErrorsSkipIndel(
        int     readLength,
        int     refLength,
        int     refOffset,
        int8_t* bases,
        int8_t* quals,
        int8_t* refForBAQ,
        int8_t* refBases,
        int     numCigarElements,
        int8_t* cigarOps,
        int*    cigarLens,
        bool    isNegativeStrand,
        bool    isExcludeFromBAQ,
        int8_t* readBAQArray,
        double* snpErrors,
        bool enableBAQ)
        //double* snpErrors,
        //double* insertErrors,
        //double* deleteErrors)
{
  // first calculate error event arrays
  int nErrors = 0;
  int readPos = 0;
  int refPos = 0;

  int8_t* isSnp = (int8_t*)calloc(readLength, sizeof(int8_t));
  int8_t* isInd = (int8_t*)calloc(readLength, sizeof(int8_t));
  int8_t* isDel = (int8_t*)calloc(readLength, sizeof(int8_t));

  for (int i = 0; i < numCigarElements; i++) {
    int index = 0;
    switch (cigarOps[i]) {
      case 'M':
      case '=':
      case 'X':
        for (int j = 0; j < cigarLens[i]; j++) {
          if (getBaseIdx(refBases[refPos]) != getBaseIdx(bases[readPos])) {
            isSnp[readPos] = 1;
            nErrors ++;
          }
          readPos++;
          refPos++;
        }
        break;
      case 'D':
        index = readPos;
        if (!isNegativeStrand) {
          index--;
        }
        if (index >= 0 && index < readLength){
          isDel[index] = 1;
          nErrors ++;
        }
        refPos += cigarLens[i];
        break;
      case 'N':
        refPos += cigarLens[i];
        break;
      case 'I':
        if (!isNegativeStrand) {
          if (readPos >= 1 && readPos-1 < readLength){
            isInd[readPos-1] = 1;
            nErrors ++;
          }
        }
        readPos += cigarLens[i];
        if (isNegativeStrand) {
          if (readPos>=0 && readPos < readLength){
            isInd[readPos] = 1;
            nErrors ++;
          }
        }
        break;
      case 'S':
        readPos += cigarLens[i];
        break;
      case 'H':
      case 'P':
        break;
      default:
        DLOG(ERROR) << "Unsupported cigar operator!";
    }
  }



  int8_t* baqArray = (int8_t*)malloc(readLength);
  // Peipei Debug:
  // TODO: please pass flag in for recalArgs.enableBAQ
  //bool enableBAQ = false;

  bool isBAQAvailable = true;
  if (!enableBAQ || nErrors == 0) { // use flatBAQArray
    const int8_t NO_BAQ_UNCERTAINTY = (int8_t)'@';
    for (int i = 0; i < readLength; i++) {
      baqArray[i] = NO_BAQ_UNCERTAINTY;
    }
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "read has no errors, use flat BAQ array";
  }
  else {
    // if read is excluded from BAQ, and read has a BAQ arrays
    // use the read baq array instead
    if (isExcludeFromBAQ) {
      if (readBAQArray) {
        memcpy(baqArray, readBAQArray, sizeof(readLength));
        DLOG_IF(INFO, VLOG_IS_ON(1)) << "read is excluded for BAQ, "
                                     << "use readBAQArray";
      }
      else {
        DLOG_IF(INFO, VLOG_IS_ON(1)) << "read is excluded for BAQ, "
                                     << "and does not have BAQ tag";
        isBAQAvailable = false;
      }
    }
    else {
      // need to check BAQ is calculatable
      if (refForBAQ == NULL) {
        DLOG_IF(INFO, VLOG_IS_ON(1)) << "This baqRead returns null";
        isBAQAvailable = false;
      }
      else {
        uint64_t start_ns = getNs();
        int baq_ret = this->calculateBAQ(baqArray,
                        bases, quals, refForBAQ,
                        cigarOps, cigarLens,
                        refLength, readLength, numCigarElements,
                        refOffset);
        baq_time_ns += getNs() - start_ns;

        if (baq_ret != 0) {
          DLOG_IF(INFO, VLOG_IS_ON(1)) << "This read is not BAQ-able";
          isBAQAvailable = false;
        }
      }
    }
  }
  if (!isBAQAvailable) {
    free(baqArray);
    free(isSnp);
    free(isInd);
    free(isDel);

    return 1;
  }
  //std::cerr << "baqArray = ";
  //std::copy(baqArray, baqArray + readLength,
  //    std::ostream_iterator<int>(std::cerr, " "));
  //    std::cerr << std::endl;

  uint64_t start_ns = getNs();
  // calculate fractional errors
  //calculateFractionalErrorArray(
  //    snpErrors, insertErrors, deleteErrors,
  //    isSnp, isInd, isDel,
  //   baqArray,
  //    readLength);

  calculateFractionalErrorArraySkipIndel(
        snpErrors,
        isSnp,
        baqArray,
        readLength);
  fracerror_time_ns += getNs() - start_ns;

  //std::cerr << "snpErrors = ";
  //std::copy(snpErrors, snpErrors + readLength,
  //    std::ostream_iterator<double>(std::cerr, " "));
  //    std::cerr << std::endl;
  //std::cerr << "insertErrors = ";
  //std::copy(insertErrors, insertErrors + readLength,
  //    std::ostream_iterator<double>(std::cerr, " "));
  //    std::cerr << std::endl;
  //std::cerr << "deleteErrors = ";
  //std::copy(deleteErrors, deleteErrors + readLength,
  //    std::ostream_iterator<double>(std::cerr, " "));
  //    std::cerr << std::endl;

  free(baqArray);
  free(isSnp);
  free(isInd);
  free(isDel);

  return 0;
}




int BAQ::calculateErrors(
        int     readLength,
        int     refLength,
        int     refOffset,
        int8_t* bases,
        int8_t* quals,
        int8_t* refForBAQ,
        int8_t* refBases,
        int     numCigarElements,
        int8_t* cigarOps,
        int*    cigarLens,
        bool    isNegativeStrand,
        bool    isExcludeFromBAQ,
        int8_t* readBAQArray,
        double* snpErrors,
        double* insertErrors,
        double* deleteErrors,
        bool enableBAQ)
{
  // first calculate error event arrays
  int nErrors = 0;
  int readPos = 0;
  int refPos = 0;

  int8_t* isSnp = (int8_t*)calloc(readLength, sizeof(int8_t));
  int8_t* isInd = (int8_t*)calloc(readLength, sizeof(int8_t));
  int8_t* isDel = (int8_t*)calloc(readLength, sizeof(int8_t));

  for (int i = 0; i < numCigarElements; i++) {
    int index = 0;
    switch (cigarOps[i]) {
      case 'M':
      case '=':
      case 'X':
        for (int j = 0; j < cigarLens[i]; j++) {
          if (getBaseIdx(refBases[refPos]) != getBaseIdx(bases[readPos])) {
            isSnp[readPos] = 1;
            nErrors ++;
          }
          readPos++;
          refPos++;
        }
        break;
      case 'D':
        index = readPos;
        if (!isNegativeStrand) {
          index--;
        }
        if (index >= 0 && index < readLength){
          isDel[index] = 1;
          nErrors ++;
        }
        refPos += cigarLens[i];
        break;
      case 'N':
        refPos += cigarLens[i];
        break;
      case 'I':
        if (!isNegativeStrand) {
          if (readPos >= 1 && readPos-1 < readLength){
            isInd[readPos-1] = 1;
            nErrors ++;
          }
        }
        readPos += cigarLens[i];
        if (isNegativeStrand) {
          if (readPos>=0 && readPos < readLength){
            isInd[readPos] = 1;
            nErrors ++;
          }
        }
        break;
      case 'S':
        readPos += cigarLens[i];
        break;
      case 'H':
      case 'P':
        break;
      default:
        DLOG(ERROR) << "Unsupported cigar operator!";
    }
  }

  //DLOG(INFO) << "nErrors = " << nErrors;
  //std::cerr << "isSNP = ";
  //std::copy(isSnp, isSnp + readLength,
  //    std::ostream_iterator<int>(std::cerr, " "));
  //    std::cerr << std::endl;

  //std::cerr << "isInd = ";
  //std::copy(isInd, isInd + readLength,
  //    std::ostream_iterator<int>(std::cerr, " "));
  //    std::cerr << std::endl;

  //std::cerr << "isDel = ";
  //std::copy(isDel, isDel + readLength,
  //    std::ostream_iterator<int>(std::cerr, " "));
  //    std::cerr << std::endl;

  int8_t* baqArray = (int8_t*)malloc(readLength);
  // Peipei Debug:
  // TODO: please pass flag in for recalArgs.enableBAQ
  //bool enableBAQ = false;

  bool isBAQAvailable = true;
  if (!enableBAQ || nErrors == 0) { // use flatBAQArray
    const int8_t NO_BAQ_UNCERTAINTY = (int8_t)'@';
    for (int i = 0; i < readLength; i++) {
      baqArray[i] = NO_BAQ_UNCERTAINTY;
    }
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "read has no errors, use flat BAQ array";
  }
  else {
    // if read is excluded from BAQ, and read has a BAQ arrays
    // use the read baq array instead
    if (isExcludeFromBAQ) {
      if (readBAQArray) {
        memcpy(baqArray, readBAQArray, sizeof(readLength));
        DLOG_IF(INFO, VLOG_IS_ON(1)) << "read is excluded for BAQ, "
                                     << "use readBAQArray";
      }
      else {
        DLOG_IF(INFO, VLOG_IS_ON(1)) << "read is excluded for BAQ, "
                                     << "and does not have BAQ tag";
        isBAQAvailable = false;
      }
    }
    else {
      // need to check BAQ is calculatable
      if (refForBAQ == NULL) {
        DLOG_IF(INFO, VLOG_IS_ON(1)) << "This baqRead returns null";
        isBAQAvailable = false;
      }
      else {
        uint64_t start_ns = getNs();
        int baq_ret = this->calculateBAQ(baqArray,
                        bases, quals, refForBAQ,
                        cigarOps, cigarLens,
                        refLength, readLength, numCigarElements,
                        refOffset);
        baq_time_ns += getNs() - start_ns;

        if (baq_ret != 0) {
          DLOG_IF(INFO, VLOG_IS_ON(1)) << "This read is not BAQ-able";
          isBAQAvailable = false;
        }
      }
    }
  }
  if (!isBAQAvailable) {
    free(baqArray);
    free(isSnp);
    free(isInd);
    free(isDel);

    return 1;
  }
  //std::cerr << "baqArray = ";
  //std::copy(baqArray, baqArray + readLength,
  //    std::ostream_iterator<int>(std::cerr, " "));
  //    std::cerr << std::endl;

  uint64_t start_ns = getNs();
  // calculate fractional errors
  calculateFractionalErrorArray(
      snpErrors, insertErrors, deleteErrors,
      isSnp, isInd, isDel,
      baqArray,
      readLength);
  fracerror_time_ns += getNs() - start_ns;

  //std::cerr << "snpErrors = ";
  //std::copy(snpErrors, snpErrors + readLength,
  //    std::ostream_iterator<double>(std::cerr, " "));
  //    std::cerr << std::endl;
  //std::cerr << "insertErrors = ";
  //std::copy(insertErrors, insertErrors + readLength,
  //    std::ostream_iterator<double>(std::cerr, " "));
  //    std::cerr << std::endl;
  //std::cerr << "deleteErrors = ";
  //std::copy(deleteErrors, deleteErrors + readLength,
  //    std::ostream_iterator<double>(std::cerr, " "));
  //    std::cerr << std::endl;

  free(baqArray);
  free(isSnp);
  free(isInd);
  free(isDel);

  return 0;
}
