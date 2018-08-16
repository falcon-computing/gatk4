#include <cstdlib>
#include <stdexcept>

#include "Common.h"
#include "RecalibrationTable.h"

RecalibrationTable::RecalibrationTable(
    int numReadGroups,
    int numEvents,
    int numCovariates,
    int* dims): math_utils_(),
                DatumTables_(numCovariates, NULL),
                RecalDatumTables_(numCovariates, NULL),
                DatumTableDimensions_(numCovariates),
                tableSizes_(numCovariates),
                numReadGroups_(numReadGroups),
                numEvents_(numEvents),
                numCovariates_(numCovariates) {

  //
  numReadsProcessed = 0;

  int qualDimension = dims[1];
  // skip the first RG table
  for (int i = 0; i < numCovariates; i++) {

    // all tables starts with numReadGroups x numEvents
    DatumTableDimensions_[i].push_back(numEvents);
    DatumTableDimensions_[i].push_back(numReadGroups);
    if (i > 0) DatumTableDimensions_[i].push_back(qualDimension);
    if (i > 1) DatumTableDimensions_[i].push_back(dims[i]);

    // allocate RecalDatum tables
    uint64_t tableSize = 1;
    for (auto dim : DatumTableDimensions_[i]) {
      tableSize *= dim;
    }
    tableSizes_[i] = tableSize;

    // skip the RG table here since we are doing recalibration
    if (i == 0) continue;
    DatumTables_[i] = (Datum*)calloc(tableSize, sizeof(Datum));
  }
  DLOG(INFO) << "Initialized RecalibrationTable";
}

RecalibrationTable::RecalibrationTable(
    int  numReadGroups,
    int  numEvents,
    int  numCovariates,
    int* dims,
    int  quantTableSize,
    int8_t* quantizationTable,
    int staticQuantizedMappingSize,
    int8_t* staticQuantizedMapping,
    bool    disableIndelQuals,
    int     preserveQLessThan,
    double  globalQScorePrior,
    bool    emitOriginalQuals):
        RecalibrationTable(numReadGroups, numEvents, numCovariates, dims)
{
  numReadsProcessed = 0;
  disableIndelQuals_ = disableIndelQuals;
  preserveQLessThan_ = preserveQLessThan;
  globalQScorePrior_ = globalQScorePrior;
  emitOriginalQuals_ = emitOriginalQuals;

  // all other fields should be initialized by overloaded constructor
  if (quantTableSize) {
    quantizationTable_ = (int8_t*)malloc(quantTableSize);
    memcpy(quantizationTable_, quantizationTable, quantTableSize);
  }
  if (staticQuantizedMapping) {
    staticQuantizedMapping_ = (int8_t*)malloc(staticQuantizedMappingSize);
    memcpy(staticQuantizedMapping_, staticQuantizedMapping, staticQuantizedMappingSize);
    DLOG(INFO) << "within native staticQuantizedMapping is not NULL, size is " << staticQuantizedMappingSize;
    for(int i = 0; i < staticQuantizedMappingSize; i++){
        DLOG(INFO) << i <<" : "<< unsigned(staticQuantizedMapping[i]);
    }
  }

  for (int i = 0; i < numCovariates; i++) {
    // free DatumTable to save space
    if (i > 0) {
      free(DatumTables_[i]);
      DatumTables_[i] = NULL;
    }

    RecalDatumTables_[i] = (RecalDatum*)calloc(tableSizes_[i], sizeof(RecalDatum));
    for (int k = 0; k < tableSizes_[i]; k++) {
      // mark empiricalQuality as uninitialized
      RecalDatumTables_[i][k].numOccurance = 0;
      RecalDatumTables_[i][k].numMismatches = .0;
      RecalDatumTables_[i][k].estimatedQReported = .0;
      RecalDatumTables_[i][k].empiricalQuality = -1.0;
    }
  }

  DLOG(INFO) << "Initialized RecalibrationTable";
}

RecalibrationTable::~RecalibrationTable() {
  // release DatumTables_
  for (auto table = DatumTables_.begin(); table != DatumTables_.end(); ++table) {
    if (*table) {
      free(*table);
      *table = NULL;
    }
  }

  // release RecalDatumTables_
  for (auto table = RecalDatumTables_.begin(); 
       table != RecalDatumTables_.end(); ++table) 
  {
    if (*table) {
      free(*table);
      *table = NULL;
    }
  }

  if (quantizationTable_) free(quantizationTable_);
  if (staticQuantizedMapping_) free(staticQuantizedMapping_);

  DLOG(INFO) << "Free RecalibrationTable";
}

void RecalibrationTable::put(RecalDatum &e, int* keys, int cov_idx) {
  if (cov_idx >= numCovariates_) {
    throw std::runtime_error("invalid table index");
  }
  // calculate index from keys
  // we need to change the layout a little, putting numEvents from
  // the last element to the first
  int num_dims = 2;

  // for covariate 1 (qual), the dimensions are 3
  if (cov_idx > 0) num_dims ++;

  // for covariates after 1, the dimensions are 4
  if (cov_idx > 1) num_dims ++;

  // shifting numEvents from last element to first
  // then start the idx calculation
  int idx = keys[num_dims - 1];
  int pitch = numEvents_;
  for (int i = 0; i < num_dims - 1; i++) {
    idx += pitch*keys[i];
    pitch *= DatumTableDimensions_[cov_idx][i+1];
  }
  if (idx >= tableSizes_[cov_idx]) {
    DLOG(ERROR) << "invalid idx = " << idx << " for cov#" << cov_idx;
  }

  RecalDatumTables_[cov_idx][idx].numOccurance       = e.numOccurance;
  RecalDatumTables_[cov_idx][idx].numMismatches      = e.numMismatches;
  RecalDatumTables_[cov_idx][idx].estimatedQReported = e.estimatedQReported;
  //RecalDatumTables_[cov_idx][idx].empiricalQuality   = e.empiricalQuality;
}

// Convert keys from computeCovariates to index to recaltables
// this function is used by both update() and recalibrate()
inline int RecalibrationTable::keysToIndex(int* keys,
      int cov_idx, int rd_idx, int event_idx)
{
  int dims[3] = {0};
  dims[0] = keys[rd_idx*numCovariates_*numEvents_ + 0*numEvents_ + event_idx];
  dims[1] = keys[rd_idx*numCovariates_*numEvents_ + 1*numEvents_ + event_idx];
  if (cov_idx > 1) {
    dims[2] = keys[rd_idx*numCovariates_*numEvents_ + cov_idx*numEvents_ + event_idx];
  }
  int idx = event_idx;
  int pitch = 1;
  int num_dims = DatumTableDimensions_[cov_idx].size()-1;

  for (int d = 0; d < num_dims; d++) {
    if (dims[d] < 0) return -1; // negative index means negative keys
    pitch *= DatumTableDimensions_[cov_idx][d];
    idx += dims[d]*pitch;
  }
  return idx;
}

void RecalibrationTable::update(int readLength,
    int*    keys,
    uint8_t*  skips,
    //uint8_t** quals,
    double**  isErrors) {
  /**
   * Covariates:
   *  - RGCov (omitted here)
   *  - QualCov: 3-dimensions, Qual x RG x Events
   *  - OptCov: 4-dimensions, Cov x Qual x RG x Events
   * keys: readLength x [RG, Qual, Cov1, Cov2] x numEvents
   *  - select two or three from keys to update different cov tables
   * skips: numEvents x readLength
   * DatumTables: numCovariates x [[Cov] x Qual x RG x Events]
   * DatumTablesDim: numCovariates x [Events, RG, Qual, [Cov]]
   */
  int saveIdx[3] = {0};
    for (int j = 0; j < readLength; j++) {
      if (skips[j]) continue;
      int newEvents=1;
       for (int k = 0; k < newEvents; k++) {
         for (int i = 1; i < numCovariates_; i++) {

      //for (int k = 0; k < numEvents_; k++) {
        //int* key = &keys[readLength*numCovariates_*k + numCovariates_*j];
        int idx = keysToIndex(keys, i, j, k);
        saveIdx[i-1]=idx;
        if (idx < 0) continue;

        DatumTables_[i][idx].numOccurance += 1;
        DatumTables_[i][idx].numMismatches += isErrors[k][j];

      }
     DLOG(INFO) << "read "<< numReadsProcessed << ", offset: "<<j<<", keys: "<< savedIdx[0] << " " << savedIdx[1] << " "<< savedIdx[2]<<" "<< " isError "<<isErrors[k][j];
    }

  }

  numReadsProcessed++;


}

// const number in bayesian calculations below
const int MAX_PHRED_SCORE = 93;
const int8_t MAX_RECALIBRATED_Q_SCORE = MAX_PHRED_SCORE;

inline double RecalibrationTable::getEmpiricalQuality(
      int cov_idx, int idx, double prior)
{
  const int SMOOTHING_CONSTANT = 1;

  RecalDatum* q = &RecalDatumTables_[cov_idx][idx];
  if (q->empiricalQuality == -1.0) {
    // calculate empiricalQuality
    // smoothing is one error and one non-error observation
    uint64_t mismatches = (uint64_t)(q->numMismatches + 0.5) + SMOOTHING_CONSTANT;
    uint64_t observations = q->numOccurance + SMOOTHING_CONSTANT + SMOOTHING_CONSTANT;

    //timer_.start(1);
    double empiricalQual = math_utils_.bayesianEstimateOfEmpiricalQuality(
                             observations, mismatches, prior);
    //timer_.stop(1);

    // This is the old and busted point estimate approach:
    //final double empiricalQual = -10 * Math.log10(getEmpiricalErrorRate());
     q->empiricalQuality = std::min(empiricalQual, (double)MAX_RECALIBRATED_Q_SCORE);
   }
   return q->empiricalQuality;
}

inline double RecalibrationTable::hierarchicalBayesianQualityEstimate(
          double epsilon, int* datum_indexes) {

  double ret = epsilon;
  double delta_prior = epsilon;

  for (int i = 0; i < numCovariates_; i++) {
    if (datum_indexes[i] < tableSizes_[i]) {
      int idx = datum_indexes[i];
      double q = idx < 0 ? 0.0 : getEmpiricalQuality(i, idx, delta_prior)
                 - delta_prior;
      ret += q;
      if (i < 2) delta_prior += q; // only need to add delta Qs for RG and Qual
    }
  }
  return ret;
}

static inline int fastRound(double d) {
  return (d > 0.0) ? (int)(d + 0.5) : (int)(d - 0.5);
}

static inline int8_t boundQual(int qual, int8_t maxQual) {
  return std::max(std::min(qual, maxQual & 0xFF), 1) & 0xFF;
}

// recalibrate read's qualities using RecalibrationTables
// quals: numEvents x readLength
//        input qualities for different events
// recal_quals: output results, allocated here
int RecalibrationTable::recalibrate(int readLength,
    int* keys,
    int8_t** quals,
    int8_t** recal_quals) {

    //DLOG(INFO)<< "within falc src, original quals: [";
    //    for (int j = 0; j < readLength; j++){
    //        DLOG(INFO)<<unsigned(quals[0][j])<<" ";
    //    }
    //DLOG(INFO)<<"]\n";
    //DLOG(INFO)<<"disableIndelQuals is "<< disableIndelQuals_<<"\n";

    for (int k = 0; k < numEvents_; k++) {
        if (disableIndelQuals_ && k > 0) {
      // skip events 1 and 2, which means:
      // BASE_INSERTION and BASE_DELETION
      continue;
    }
    // the rg key is constant over the whole read, the global deltaQ is too
    int rgIdx = keysToIndex(keys, 0, 0, k); // keys[0][0][event]
    if (rgIdx >= tableSizes_[0] || rgIdx < 0) {
      DLOG(INFO) << "invalid read group";
      // equivalent to empiricalQualRG == null
      memcpy(recal_quals[k], quals[k], readLength);
      continue;
    }
    RecalDatum empiricalQualRG = RecalDatumTables_[0][rgIdx];
    double epsilon = (globalQScorePrior_ > 0.0 &&
                      k == 0 // event = BASE_SUBSTITUTION
                     ) ? globalQScorePrior_ :
                         empiricalQualRG.estimatedQReported;

    //DLOG(INFO) << "epsilon = " << epsilon;

    int* empirical_quals_idx = (int*)malloc(numCovariates_*sizeof(int));
    for (int j = 0; j < readLength; j++) {
      if (quals[k][j] < preserveQLessThan_) {
        recal_quals[k][j] = quals[k][j];
        continue;
      }

      // here in hierarchicalBayesianQualityEstimate() function,
      // the values of RecalDatum will be modified, and reused in future
      // references, therefore we pass the idx rather than actual values
      // to the function.
      for (int i = 0; i < numCovariates_; i++) {
        empirical_quals_idx[i] = keysToIndex(keys, i, j, k);
      }
      //timer_.start(0);
      double recalibratedQualDouble = hierarchicalBayesianQualityEstimate(
            epsilon, empirical_quals_idx);
      //timer_.stop(0);
      //DLOG(INFO) << "recalibratedQualDouble = " << recalibratedQualDouble;


      // recalibrated quality is bound between 1 and MAX_QUAL
      int8_t recalibratedQual = boundQual(fastRound(recalibratedQualDouble),
                                          MAX_RECALIBRATED_Q_SCORE);

      // return the quantized version of the recalibrated quality
      int8_t recalibratedQualityScore = quantizationTable_[recalibratedQual];

      // Bin to static quals
      if(staticQuantizedMapping_ != NULL) {
        recal_quals[k][j] = staticQuantizedMapping_[recalibratedQualityScore];
      }
      else {
        recal_quals[k][j] = recalibratedQualityScore;
      }
      //offset: 62, recalibratedQualDouble: 33.000000, recalibratedQualityScore: 33
      //DLOG(INFO) << "offset: "<< j << ", recalibratedQualDouble: "<<recalibratedQualDouble<<", recalibratedQual: "<<unsigned(recalibratedQual)<<", recalibratedQualityScore: "<<unsigned(recalibratedQualityScore)<<"\n";
    }

    //DLOG(INFO) << (staticQuantizedMapping_ == NULL) ;
    free(empirical_quals_idx);
  }
  //DLOG_EVERY_N(INFO, TIMER_SAMPLE_RATE)
  //                    << "Running " << timer_.get_laps(0) << "times of "
  //                    << "hierarchicalBayesianQualityEstimate() takes "
  //                    << timer_.get_ms(0) << " ms";
  //DLOG_EVERY_N(INFO, TIMER_SAMPLE_RATE)
  //                    << "Running " << timer_.get_laps(1) << "times of "
  //                    << "bayesianEstimateOfEmpiricalQuality() takes "
  //                    << timer_.get_ms(1) << " ms";

  return 0;
}

Datum* RecalibrationTable::getTable(int idx) {
  if (idx >= (int)DatumTables_.size()) {
    throw std::runtime_error("invalid table size");
  }
  return DatumTables_[idx];
}

RecalDatum* RecalibrationTable::getFullTable(int idx) {
  if (idx >= (int)RecalDatumTables_.size()) {
    throw std::runtime_error("invalid table size");
  }
  return RecalDatumTables_[idx];
}

int RecalibrationTable::getTableSize(int idx) {
  return tableSizes_[idx];
}

std::vector<int> RecalibrationTable::getTableDimensions(int idx) {
  return DatumTableDimensions_[idx];
}
