#ifndef BQSR_BAQTEST_H
#define BQSR_BAQTEST_H

#include <iostream>
#include <gtest/gtest.h>
#include <glog/logging.h>
#include <string>

#include "BAQ.h"

class BAQTest : public ::testing::Test {
  protected:
    ;
};

int string_split(std::string str, int8_t* dst) {
  int length = 0;
  std::stringstream ss(str);
  while (!ss.eof()) {
    ss >> dst[length++];
  }
  return length;
}

int string_split_num(std::string str, int8_t* dst) {
  int length = 0;
  std::stringstream ss(str);
  while (!ss.eof()) {
    std::string num;
    ss >> num;
    dst[length++] = (char)std::stoi(num);
  }
  return length;
}

int string_split_num(std::string str, double* dst) {
  int length = 0;
  std::stringstream ss(str);
  while (!ss.eof()) {
    std::string num;
    ss >> num;
    dst[length++] = (double)std::stod(num);
  }
  return length;
}

TEST_F(BAQTest, Test_hmmglocal) {
  BAQ* baq = new BAQ();

  // prepare input
  int readLength = 101;
  int refLength = 104;
  int queryStart = 0;
  int queryLen = 97;
  std::string str_refBases = "C C C T A A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A A C C";
  std::string str_bases = "T A A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A C C C T A A C C C T A A C C C T A A C C C T A A C C C T A C C C C T A A C C C T A A C C C T A";
  std::string str_quals = "34 34 34 37 37 37 37 37 38 39 39 37 39 40 41 40 41 40 41 41 41 40 40 40 41 38 40 40 41 38 38 40 40 41 40 39 33 37 35 39 38 40 34 34 38 38 38 38 34 38 40 40 39 28 37 31 21 31 34 24 30 32 28 7 21 7 22 22 30 32 30 29 29 26 28 24 26 26 30 7 11 7 24 20 20 27 30 32 32 7 18 7 23 27 27 27 30";
  std::string str_baq = "46 62 77 81 79 81 94 84 84 81 79 81 94 84 84 81 79 81 94 84 84 81 79 81 94 84 84 81 79 81 94 84 84 81 79 81 94 84 84 81 79 81 94 83 73 70 68 67 42 4 43 44 47 65 65 65 65 64 65 65 65 65 65 64 65 65 65 65 65 64 65 65 65 65 65 64 65 65 63 63 63 63 64 65 65 65 65 64 65 65 65 65 65 64 64 60 39";

  int8_t* refBases = (int8_t*)malloc(refLength);
  int8_t* bases = (int8_t*)malloc(readLength);
  int8_t* quals = (int8_t*)malloc(readLength);
  int*  state = (int*)malloc(readLength*sizeof(int));
  int8_t* q = (int8_t*)calloc(readLength, sizeof(int8_t));
  int8_t* q_base = (int8_t*)calloc(readLength, sizeof(int8_t));

  string_split(str_refBases, refBases);
  string_split(str_bases, bases);
  string_split_num(str_quals, quals);
  string_split_num(str_baq, q_base);

  for (int i = 0; i < 10; i++) {
    baq->hmm_glocal(refLength,
          refBases, bases,
          queryStart, queryLen,
          quals, state, q);
  }

  for (int i = 0; i < queryLen; i++) {
    ASSERT_EQ(q[i], q_base[i]);
  }

  free(refBases);
  free(bases);
  free(quals);
  free(state);
  free(q);
  free(q_base);

  delete baq;
}

TEST_F(BAQTest, Test_CalcErrorArray) {
  BAQ* baq = new BAQ();

  int maxReadLength = 101;
  std::string str_isSNP = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
  std::string str_isInd = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
  std::string str_isDel = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
  std::string str_baqArray = "64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 72 64 95 93 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64";

  std::string str_snpErrors = "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0";
  std::string str_insertionErrors = "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0";
  std::string str_deletionErrors = "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3333333333333333 0.3333333333333333 0.25 0.25 0.25 0.25 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0";

  int8_t* isSnp = (int8_t*)malloc(maxReadLength);
  int8_t* isInd = (int8_t*)malloc(maxReadLength);
  int8_t* isDel = (int8_t*)malloc(maxReadLength);
  int8_t* baqArray = (int8_t*)malloc(maxReadLength);

  double* snpErrors       = (double*)calloc(maxReadLength, sizeof(double));
  double* insertionErrors = (double*)calloc(maxReadLength, sizeof(double));
  double* deletionErrors  = (double*)calloc(maxReadLength, sizeof(double));

  int readLength = 0;
  DLOG(INFO) << "isSnp" << "length = " << (readLength = string_split_num(str_isSNP, isSnp));
  ASSERT_EQ(readLength, string_split_num(str_isInd, isInd));
  ASSERT_EQ(readLength, string_split_num(str_isDel, isDel));
  ASSERT_EQ(readLength, string_split_num(str_baqArray, baqArray));

  baq->calculateFractionalErrorArray(snpErrors, insertionErrors, deletionErrors,
          isSnp, isInd, isDel, baqArray, readLength);

  double* base1 = (double*)malloc(readLength*sizeof(double));
  double* base2 = (double*)malloc(readLength*sizeof(double));
  double* base3 = (double*)malloc(readLength*sizeof(double));

  DLOG(INFO) << "snpErrors length = " << string_split_num(str_snpErrors, base1);
  DLOG(INFO) << "insertErrors length = " << string_split_num(str_insertionErrors, base2);
  DLOG(INFO) << "deletionErrors length = " << string_split_num(str_deletionErrors, base3);

  for (int i = 0; i < readLength; i++) {
    ASSERT_EQ(base1[i], snpErrors[i]) << "i = " << i;
    ASSERT_EQ(base2[i], insertionErrors[i]) << "i = " << i;
    ASSERT_EQ(base3[i], deletionErrors[i]) << "i = " << i;
  }

  delete baq;

  free(isSnp);
  free(isInd);
  free(isDel);
  free(baqArray);
  free(snpErrors);
  free(insertionErrors);
  free(deletionErrors);
}
#endif
