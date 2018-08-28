#ifndef SMITHWATERMAN_COMMON_H
#define SMITHWATERMAN_COMMON_H

#if defined(_MSC_VER)
  #include <intrin.h> // SIMD intrinsics for Windows
#else
  #include <x86intrin.h> // SIMD intrinsics for GCC
#endif

#define __STDC_LIMIT_MACROS

#include <stdint.h>
#include <string.h>
#include <immintrin.h>


#define CAT(X,Y) X##Y
#define CONCAT(X,Y) CAT(X,Y)

#define MATCH 0
#define INSERT 1
#define DELETE 2
#define INSERT_EXT 4
#define DELETE_EXT 8
#define SOFTCLIP 9
#define INDEL 10
#define LEADING_INDEL 11
#define IGNORE 12

#define OVERHANG_STRATEGY_SOFTCLIP 0
#define OVERHANG_STRATEGY_INDEL 1
#define OVERHANG_STRATEGY_LEADING_INDEL 2
#define OVERHANG_STRATEGY_IGNORE 3

#define STATE_MATCH 0       //equal to CigarOperator.M
#define STATE_INSERTION 1   //equal to CigarOperator.I
#define STATE_DELETION 2    //equal to CigarOperator.D
#define STATE_CLIP 4        //equal to CigarOperator.S

#define W_MATCH 200
#define W_MISMATCH -150
#define W_OPEN -260
#define W_EXTEND -11

#define MAX_SEQ_LENGTH 1536
struct CigarElement{
    int length;
    int state;
};

struct Cigar{
    struct CigarElement cigarElements[MAX_SEQ_LENGTH];
    int CigarElementNum;
};

struct sw_ptrs{
    int32_t* E_;
    int16_t* backTrack_;
    int16_t* cigarBuf_;
};



typedef struct dnaSeqPair
{
        int32_t id;
        uint8_t *seq1;
        uint8_t *seq2;
        int16_t len1, len2;
        int8_t overhangStrategy;
        int32_t score;
        int16_t max_i;
        int16_t max_j;
        bool update_max_j;
        int16_t *btrack;
        struct Cigar* cigar;
        int16_t alignmentOffset;
}SeqPair;

//the maximum DNA sequence length
#define MAX_SEQ_LEN MAX_SEQ_LENGTH
#define MAX_NUM_PAIRS 800000
#define MATRIX_MIN_CUTOFF -100000000
#define LOW_INIT_VALUE (INT32_MIN/2)
#define max(x, y) ((x)>(y)?(x):(y))
#define min(x, y) ((x)<(y)?(x):(y))
#define DUMMY1 'B'
#define DUMMY2 'D'


#endif
