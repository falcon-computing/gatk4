#ifndef BQSR_BAQ_H
#define BQSR_BAQ_H
#include <cstdint>
#include "gtest/gtest_prod.h"
#include "Common.h"

class BAQ {
  public:
    BAQ();
    /*
     * Same function as BaseRecalibrator.calculateBAQArray()
     * - Input:
     *     refBases, bases, baseQuals
     *     cigarOps, cigarLens
     *     refOffset
     * - Output:
     *     bqTag
     * - Return:
     *     true -- normal
     *     false -- null
     */
    int calculateBAQ(
        int8_t* bqTag,
        int8_t* bases,
        int8_t* baseQuals,
        int8_t* refBases,
        int8_t* cigarOps,
        int* cigarLens,
        int refLength,
        int readLength,
        int numCigar,
        int refOffset);

    /*
         * Same function as BaseRecalibrator.calculateFractionalErrorArray()
         * and also contains this->calculate(), which calculates BAQ Before
         * calculating fractional errors.
         */
        int calculateErrorsSkipIndel(
            int readLength,
            int refLength,
            int refOffset,
            int8_t* bases,
            int8_t* quals,
            int8_t* refForBAQ,
            int8_t* refBases,
            int     numCigarElements,
            int8_t* cigarOps,
            int* cigarLens,
            bool isNegativeStrand,
            bool isExcludeFromBAQ,
            int8_t* readBAQArray,
            double* snpErrors,
            bool enableBAQ);
            //double* snpErrors,
            //double* insertErrors,
            //double* deleteErrors);


    /*
     * Same function as BaseRecalibrator.calculateFractionalErrorArray()
     * and also contains this->calculate(), which calculates BAQ Before
     * calculating fractional errors.
     */
    int calculateErrors(
        int readLength,
        int refLength,
        int refOffset,
        int8_t* bases,
        int8_t* quals,
        int8_t* refForBAQ,
        int8_t* refBases,
        int     numCigarElements,
        int8_t* cigarOps,
        int* cigarLens,
        bool isNegativeStrand,
        bool isExcludeFromBAQ,
        int8_t* readBAQArray,
        double* snpErrors,
        double* insertErrors,
        double* deleteErrors,
        bool enableBAQ);

    ~BAQ();

  protected:
    FRIEND_TEST(BAQTest, Test_hmmglocal);
    FRIEND_TEST(BAQTest, Test_CalcErrorArray);

    inline int8_t capBaseByBAQ(int8_t oq, int8_t bq, int state, int expectedPos);

    /*
     * THIS CODE IS SYNCHRONIZED WITH CODE IN THE SAMTOOLS REPOSITORY.
     *
     * The following arrays all have length base_len:
     *   query: read bases
     *   iqual: base qualities
     *   state: return state
     *   q: return qual
     *
     * The following arrays have length ref_len
     *   ref: reference bases
     */
    int hmm_glocal(
          int ref_len,
          int8_t* ref,
          int8_t* query,
          int qstart,
          int l_query,
          int8_t* _iqual,
          int* state,
          int8_t* q);

    inline double calcEpsilon(int8_t ref, int8_t read, int8_t qualB) {
      return EPSILONS[ref][read][qualB];
    }

    void calculateFractionalErrorArraySkipIndel(
              double* snpErrors,
              //double* insertErrors,
              //double* deleteErrors,
              int8_t* snpArray,
              //int8_t* insertArray,
              //int8_t* deleteArray,
              int8_t* baqArray,
              int readLength);

    void calculateFractionalErrorArray(
          double* snpErrors,
          double* insertErrors,
          double* deleteErrors,
          int8_t* snpArray,
          int8_t* insertArray,
          int8_t* deleteArray,
          int8_t* baqArray,
          int readLength);

    double cd;     // gap open probability [1e-3]
    double ce;     // gap extension probability [0.1]
    int cb;        // band width [7]

    int MAX_PHRED_SCORE = 93;
    double EM = 0.33333333333;
    double EI = 0.25;

    // 94 = MAX_PHRED_SCORE + 1
    double EPSILONS[256][256][94];

    bool   includeClippedBases;
    int8_t minBaseQual;

    // performance meters
    uint64_t hmm_time_ns;
    uint64_t baq_time_ns;
    uint64_t fracerror_time_ns;
};

#endif
