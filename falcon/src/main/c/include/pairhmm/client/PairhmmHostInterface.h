#ifndef PAIRHMM_HOST_INTERFACE_H
#define PAIRHMM_HOST_INTERFACE_H
#include <cstdint>
#include <string>
#include <vector>
#include <time.h>
#include <string>

#include "pairhmm/common/headers.h"

#define MAX_READ_LEN 192
#define MAX_HAP_LEN 1024
#define MAX_RSDATA_NUM 2048   //default 4096, minimum 32
#define MAX_HAPDATA_NUM 128   
#define DEP_DIST 42
#define READ_BLOCK_SIZE 2
#define HAP_BLOCK_SIZE 4

#define FPGA_perf 18
#define AVX_perf 0.3


#define TRANS_PROB_ARRAY_LENGTH 6

#define TRANSITION_matchToMatch 0
#define TRANSITION_indelToMatch 1
#define TRANSITION_matchToInsertion 2
#define TRANSITION_insertionToInsertion 3
#define TRANSITION_matchToDeletion 4
#define TRANSITION_deletionToDeletion 5

#define MM 0
#define GapM 1
#define MX 2
#define XX 3
#define MY 4
#define YY 5


#define MAX_QUAL 254
#define MAX_JACOBIAN_TOLERANCE 8.0
#define JACOBIAN_LOG_TABLE_STEP 0.0001
#define JACOBIAN_LOG_TABLE_INV_STEP (1.0 / JACOBIAN_LOG_TABLE_STEP)
#define MAXN 70000
#define LOG10_CACHE_SIZE  (4*MAXN)  // we need to be able to go up to 2*(2N) when calculating some of the coefficients
#define JACOBIAN_LOG_TABLE_SIZE ((int) (MAX_JACOBIAN_TOLERANCE / JACOBIAN_LOG_TABLE_STEP) + 1)

#define SET_MATCH_TO_MATCH_PROB(output, insQual, delQual)                       \
{                                                                               \
      int minQual = delQual;                                                        \
      int maxQual = insQual;                                                        \
      if (insQual <= delQual)                                                       \
      {                                                                             \
              minQual = insQual;                                                          \
              maxQual = delQual;                                                          \
            }                                                                             \
      (output) = (MAX_QUAL < maxQual) ?                                             \
      ((NUMBER)1.0) - ctx.POW(((NUMBER)10), ctx.approximateLog10SumLog10(((NUMBER)-0.1)*minQual, ((NUMBER)-0.1)*maxQual))       \
      : ctx.matchToMatchProb[((maxQual * (maxQual + 1)) >> 1) + minQual];           \
}


typedef struct {
  int   len;
  const char* _b;
  const char* _q;
  const char* _i;
  const char* _d;
  const char* _c;
} read_t;

typedef struct {
  int len;
  const char* _b;
} hap_t;

uint64_t serialize(void* buf, const read_t* reads, int num);
uint64_t serialize(void* buf, const hap_t* haps, int num);

#endif
