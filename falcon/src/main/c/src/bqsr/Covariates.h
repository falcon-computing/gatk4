#ifndef BQSR_COVARIATES_H
#define BQSR_COVARIATES_H
#include <cstdint>
#include <unordered_map>
#include <string>

class Covariates {

  public:
    Covariates(int _numEvents,
      int _numCovariates,
      int _mismatchesContextSize,
      int _indelsContextSize,
      uint8_t _lowQualTail,
      int _maxCycleValue,
      int _cushionForIndels);

    void initReadGroup(const std::string readGroup);

    void compute(int* keys,
      const int readLength,
      const std::string readGroup,
      const int8_t* bases,
      const int8_t* basesQuals,
      const int8_t* insertionQuals,
      const int8_t* deleteionQuals,
      const int platformType,
      const bool isNegativeStrand,
      const bool isReadPaired,
      const bool isSecondOfPair);

    void computeReadGroupCovariates(int* keys,
      const int readLength,
      const std::string readGroup);

    void computeQualCovariates(int* keys,
      const int readLength,
      const int8_t* basesQuals,
      const int8_t* insertionQuals,
      const int8_t* deleteionQuals);

    void computeContextCovariates(int* keys,
      const int readLength,
      const bool isNegativeStrand,
      const int8_t* bases,
      const int8_t* quals);

    void computeCycleCovariates(int* keys,
      const int readLength,
      const int platformType,
      const bool isNegativeStrand,
      const bool isReadPaired,
      const bool isSecondOfPair);

    const std::vector<std::string>& getRGKeys() const;

  private:
    inline void setCovariate(
      int* keys,
      const int cov_idx,
      const int base_idx,
      const int snp,
      const int insertion,
      const int deletion);

    inline int keyFromCycle(int cycle);

    // TODO: this part needs to be locked
    std::unordered_map<std::string, int> readGroupTable;
    std::vector<std::string> readGroupIdxTable;
    int readGroupIdx;

    int numEvents;
    int numCovariates;
    int mismatchesContextSize;
    int indelsContextSize;
    int mismatchesKeyMask;
    int indelsKeyMask;

    uint8_t LOW_QUAL_TAIL;
    int MAX_CYCLE_VALUE;
    int CUSHION_FOR_INDELS;
};
#endif
