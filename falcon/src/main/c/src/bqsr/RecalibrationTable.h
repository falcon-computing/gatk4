#include <cstdint>
#include <vector>

#include "MathUtils.h"
#include "SimpleTimer.h"

// This struct is used table construct to save space
typedef struct {
  uint64_t numOccurance;
  double   numMismatches;
} Datum;

// this struct is used in table lookup after table is constructed
typedef struct {
  uint64_t numOccurance;
  double   numMismatches;
  double   estimatedQReported;
  double   empiricalQuality;
} RecalDatum;

class RecalibrationTable {
  public:
    // construct an empty table for the recalibrator
    RecalibrationTable(int  numReadGroups,
                       int  numEvents,
                       int  numCovariates,
                       int* dims);

    RecalibrationTable(int  numReadGroups,
                       int  numEvents,
                       int  numCovariates,
                       int* dims,
                       int quantTableSize,
                       int8_t* quantizationTable,
                       int staticQuantizedMappingSize,
                       int8_t* staticQuantizedMapping,
                       bool disableIndelQuals,
                       int preserveQLessThan,
                       double globalQScorePrior,
                       bool emitOriginalQuals);

    ~RecalibrationTable();

    // insert a recal datum into the
    // - keys: the covariates organizing using GATK format,
    //         rg, qual, cov1/cov2, numEvents_
    // - cov_idx: index to the covariates
    void put(RecalDatum& e, int* keys, int cov_idx);

    // histogram from input read
    // used in GATK BaseRecalibrator
    void update(int      readLength,
                int*     keys,
                uint8_t* skips,
                //uint8_t** quals,
                double** isErrors);

    // recalibrate read by contents in the table
    // used in GATK PrintReads
    int recalibrate(int readLength,
                    int* keys,       // computed by Covariates.compute()
                    int8_t** quals,
                    int8_t** recal_quals); // output quals (snp, indel, delete)

    Datum* getTable(int idx);
    RecalDatum* getFullTable(int idx);

    int getTableSize(int idx);
    std::vector<int> getTableDimensions(int idx);

  private:

    // helper function for update() and recalibrate()
    inline int keysToIndex(int* keys,
        int cov_idx, int rd_idx, int event_idx);

    // helper functions for recalibrate()
    inline double getEmpiricalQuality(
        int cov_idx, int idx, double prior);

    inline double hierarchicalBayesianQualityEstimate(
        double epsilon, int* datum_indexes);

    // used for math routines in this class
    MathUtils math_utils_;

    std::vector<Datum*> DatumTables_;
    std::vector<RecalDatum*> RecalDatumTables_;
    std::vector<std::vector<int> > DatumTableDimensions_;
    std::vector<int> tableSizes_;

    int numReadGroups_;
    int numEvents_;
    int numCovariates_;

    int8_t* quantizationTable_      = NULL;
    int8_t* staticQuantizedMapping_ = NULL;

    // parameters for recalibrate
    bool   disableIndelQuals_;
    int    preserveQLessThan_;
    double globalQScorePrior_;
    bool   emitOriginalQuals_;

    SimpleTimer timer_;
};
