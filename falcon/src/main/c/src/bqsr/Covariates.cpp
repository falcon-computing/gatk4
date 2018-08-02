#include <cstdlib>
#include <stdexcept>

#include "Common.h"
#include "Covariates.h"
#include <iostream>
#include <string>

const int LENGTH_BITS = 4;
const int LENGTH_MASK = 15;

static inline int8_t simpleComplement(int8_t base) {
  switch (base) {
    case 'A':
    case 'a':
      return 'T';
    case 'C':
    case 'c':
      return 'G';
    case 'G':
    case 'g':
      return 'C';
    case 'T':
    case 't':
      return 'A';
    default:
      return base;
  }
}

static inline int simpleBaseToBaseIndex(int8_t base) {
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
    default: ;
  }
  return -1;
}

static inline int keyFromContext(int8_t* dna, int start, int end) {
  int key = end - start;
  int bitOffset = LENGTH_BITS;
  for (int i = start; i < end; i++) {
    int baseIndex = simpleBaseToBaseIndex(dna[i]);
    if (baseIndex == -1) { // ignore non-ACGT bases
      return -1;
    }
    key |= (baseIndex << bitOffset);
    bitOffset += 2;
  }
  return key;
}

static inline void contextWith(
    int* keys,
    int readLength,
    int8_t* bases,
    int contextSize,
    int mask)
{
  memset(keys, 0, readLength*sizeof(int));

  int tempIndex = 0;
  for (int i = 1; (i < contextSize && i <= readLength); i++) {
    keys[tempIndex++] = -1;
  }

  if (readLength < contextSize) {
    return;
  }

  int newBaseOffset = 2 * (contextSize - 1) + LENGTH_BITS;

  // get (and add) the key for the context starting at the first base
  int currentKey = keyFromContext(bases, 0, contextSize);
  keys[tempIndex++] = currentKey;

  // if the first key was -1 then there was an N in the context; figure out how many more consecutive contexts it affects
  int currentNPenalty = 0;

  if (currentKey == -1) {
    currentKey = 0;
    currentNPenalty = contextSize - 1;
    int offset = newBaseOffset;

    while (bases[currentNPenalty] != 'N') {
      int baseIndex = simpleBaseToBaseIndex(bases[currentNPenalty]);
      currentKey |= (baseIndex << offset);
      offset -= 2;
      currentNPenalty--;
    }
  }
  for (int currentIndex = contextSize; currentIndex < readLength; currentIndex++) {
    int baseIndex = simpleBaseToBaseIndex(bases[currentIndex]);
    if (baseIndex == -1) { // ignore non-ACGT bases
      currentNPenalty = contextSize;
      currentKey = 0; // reset the key
    }
    else {
      // push this base's contribution onto the key: shift everything 2 bits, mask out the non-context bits, and add the new base and the length in
      currentKey = (currentKey >> 2) & mask;
      currentKey |= (baseIndex << newBaseOffset);
      currentKey |= contextSize;
    }

    if (currentNPenalty == 0) {
      keys[tempIndex++] = currentKey;
    }
    else {
      currentNPenalty --;
      keys[tempIndex++] = -1;
    }
  }
}

Covariates::Covariates(int _numEvents,
      int _numCovariates,
      int _mismatchesContextSize,
      int _indelsContextSize,
      uint8_t _lowQualTail,
      int _maxCycleValue,
      int _cushionForIndels): readGroupIdx(0),
                              numEvents(_numEvents),
                              numCovariates(_numCovariates),
                              mismatchesContextSize(_mismatchesContextSize),
                              indelsContextSize(_indelsContextSize),
                              mismatchesKeyMask(0),
                              indelsKeyMask(0),
                              LOW_QUAL_TAIL(_lowQualTail),
                              MAX_CYCLE_VALUE(_maxCycleValue),
                              CUSHION_FOR_INDELS(_cushionForIndels)
{
  // compute key masks
  for (int i = 0; i < mismatchesContextSize; i++) {
    mismatchesKeyMask = (mismatchesKeyMask << 2) | 3;
  }
  mismatchesKeyMask = mismatchesKeyMask << LENGTH_BITS;
  for (int i = 0; i < indelsContextSize; i++) {
      indelsKeyMask = (indelsKeyMask << 2) | 3;
  }
  indelsKeyMask = indelsKeyMask << LENGTH_BITS;
}

// helper function to set output covariates
inline void Covariates::setCovariate(
    int* keys,
    const int cov_idx,
    const int base_idx,
    const int snp,
    const int insertion,
    const int deletion) {

  // in this version the layout of keys is the same as table update
  // keys: readLength x [qual context, cycle] x numEvents
  keys[base_idx*numCovariates*numEvents + cov_idx*numEvents + 0] = snp;
  keys[base_idx*numCovariates*numEvents + cov_idx*numEvents + 1] = insertion;
  keys[base_idx*numCovariates*numEvents + cov_idx*numEvents + 2] = deletion;
}


void Covariates::compute(int* keys,
      const int readLength,
      const std::string readGroup,
      const int8_t* bases,
      const int8_t* baseQuals,
      const int8_t* insertionQuals,
      const int8_t* deletionQuals,
      const int platformType,
      const bool isNegativeStrand,
      const bool isReadPaired,
      const bool isSecondOfPair) {

  computeReadGroupCovariates(keys, readLength, readGroup);
  computeQualCovariates(keys, readLength,
     baseQuals, insertionQuals, deletionQuals);
  computeContextCovariates(keys, readLength,
     isNegativeStrand, bases, baseQuals);
  computeCycleCovariates(keys, readLength, platformType,
     isNegativeStrand, isReadPaired, isSecondOfPair);
}

void Covariates::initReadGroup(const std::string readGroup) {

  if (readGroupTable.count(readGroup)) {
    //throw std::runtime_error("duplicated Read Group name");
    LOG(WARNING) << "Duplicated read group name";
    return;
  }
  readGroupTable[readGroup] = readGroupIdx;
  readGroupIdx++;

  // record the idx for the string
  readGroupIdxTable.push_back(readGroup);
}

//added by Peipei in bqsr init table
void Covariates::addReadGroup(const std::string readGroup){
    int readGroupId;
    if (readGroupTable.count(readGroup)){
        readGroupId = readGroupTable[readGroup];
    }
    else{
        std::cout<<"First time See: "<<readGroup<<" map to Id: "readGroupIdx<<std::endl;
        readGroupTable[readGroup] = readGroupIdx;
        readGroupId = readGroupIdx;
        readGroupIdx++;
        readGroupIdxTable.push_back(readGroup);

    }
}

void Covariates::computeReadGroupCovariates(int* keys,
    const int readLength,
    const std::string readGroup) {

  int readGroupId;
  if (readGroupTable.count(readGroup)) {
    readGroupId = readGroupTable[readGroup];
  }
  else {
    readGroupTable[readGroup] = readGroupIdx;
    readGroupId = readGroupIdx;
    readGroupIdx++;

    // record the idx for the string
    readGroupIdxTable.push_back(readGroup);
  }
  for (int i = 0; i < readLength; i++) {
    setCovariate(keys, 0, i,
        readGroupId, readGroupId, readGroupId);
  }
}

void Covariates::computeQualCovariates(int* keys,
    const int readLength,
    const int8_t* baseQuals,
    const int8_t* insertionQuals,
    const int8_t* deleteionQuals) {

  for (int i = 0; i < readLength; i++) {
    setCovariate(keys, 1, i,
        baseQuals[i],
        insertionQuals[i],
        deleteionQuals[i]);
  }
}

void Covariates::computeContextCovariates(int* keys,
    const int readLength,
    const bool isNegativeStrand,
    const int8_t* bases,
    const int8_t* quals)
{
  // clip reads
  int8_t* clipped_bases = (int8_t*)malloc(readLength);

  int leftClipIndex = 0;
  int rightClipIndex = readLength - 1;

  // check how far we can clip both sides
  while (rightClipIndex >= 0 && quals[rightClipIndex] <= LOW_QUAL_TAIL)
    rightClipIndex --;

  while (rightClipIndex < readLength && quals[leftClipIndex] <= LOW_QUAL_TAIL)
    leftClipIndex ++;

  // the read is empty, set all covariates to zero
  if (leftClipIndex > rightClipIndex) {
    for (int i = 0; i < readLength; i++) {
      setCovariate(keys, 2, i, 0, 0, 0);
    }
    free(clipped_bases);
    return;
  }

  // clip the low quality reads on both ends
  for (int i = 0; i < readLength; i++) {
      if (i < leftClipIndex || i > rightClipIndex)
        clipped_bases[i] = 'N';
      else
        clipped_bases[i] = bases[i];
  }

  // reverse read if necessary
  if (isNegativeStrand) {
    int8_t* temp_bases = (int8_t*)malloc(readLength);
    memcpy(temp_bases, clipped_bases, readLength);

    for (int i = 0; i < readLength; i++) {
      clipped_bases[i] = simpleComplement(temp_bases[readLength - i - 1]);
    }

    free(temp_bases);
  }

  int* mismatchKeys = (int*)malloc(readLength*sizeof(int));
  int* indelKeys = (int*)malloc(readLength*sizeof(int));
  contextWith(mismatchKeys, readLength, clipped_bases, mismatchesContextSize, mismatchesKeyMask);
  contextWith(indelKeys, readLength, clipped_bases, indelsContextSize, indelsKeyMask);

  //Note: duplicated the loop to avoid checking recordIndelValues on each iteration
  for (int i = 0; i < readLength; i++) {
    int readOffset = i;
    if (isNegativeStrand){
      readOffset = readLength - i - 1;
    }
    setCovariate(keys, 2, readOffset, mismatchKeys[i], indelKeys[i], indelKeys[i]);
  }
  free(clipped_bases);
  free(mismatchKeys);
  free(indelKeys);
}

inline int Covariates::keyFromCycle(int cycle) {
  // no negative values because values must fit
  // into the first few bits of the long
  int result = cycle;
  if (result < 0){
    result = -result;
  }
  if (result > MAX_CYCLE_VALUE) {
    return -1;
  }

  result <<= 1; // shift so we can add the "sign" bit
  if (cycle < 0) {
    result++; // negative cycles get the lower-most bit set
  }
  return result;
}

void Covariates::computeCycleCovariates(int* keys,
    const int readLength,
    const int platformType,
    const bool isNegativeStrand,
    const bool isReadPaired,
    const bool isSecondOfPair)
{

  // TODO: need to do flow type platforms as well
  if (platformType != 0) { // Currently only supports DISCRETE platform
    LOG(WARNING) << "Only DISCRETE platforms are supported in this version.";
    throw std::runtime_error("unsupported platform");
  }

  int readOrderFactor = isReadPaired && isSecondOfPair ? -1 : 1;
  int increment;
  int cycle;

  if (isNegativeStrand) {
    cycle = readLength * readOrderFactor;
    increment = -1 * readOrderFactor;
  } else {
    cycle = readOrderFactor;
    increment = readOrderFactor;
  }
  int MAX_CYCLE_FOR_INDELS = readLength - CUSHION_FOR_INDELS - 1;

  for (int i = 0; i < readLength; i++) {
    int substitutionKey = keyFromCycle(cycle);
    int indelKey = (i < CUSHION_FOR_INDELS || i > MAX_CYCLE_FOR_INDELS) ? -1 : substitutionKey;
    setCovariate(keys, 3, i, substitutionKey, indelKey, indelKey);
    cycle += increment;
  }
}

const std::vector<std::string>& Covariates::getRGKeys() const {
  return readGroupIdxTable;
}
