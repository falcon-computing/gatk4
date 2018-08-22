#include <string>

#include "com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine.h"
#include "Common.h"
#include "BAQ.h"
#include "Covariates.h"
#include "MathUtils.h"
#include "RecalibrationTable.h"
#include "SimpleTimer.h"
#include <iostream>

int g_numReadGroups = 0;
int g_numEvents = 0;
int g_numCovariates = 0;

RecalibrationTable* table = NULL;
Covariates*         cov   = NULL;
BAQ*                baq   = NULL;

// for performance measuring
uint64_t total_baq_time = 0;
uint64_t total_baq_compute_time = 0;

uint64_t total_update_baq_time = 0;
uint64_t total_update_covariate_time = 0;
uint64_t total_update_compute_time = 0;
uint64_t total_update_total_time = 0;
int      total_update_num_calls = 0;

int lap_recalibrate_num_calls = 0;
int total_recalibrate_num_calls = 0;

SimpleTimer timer;

// init function for BaseRecalibrator
JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_initNative__III_3IBIIII(
    JNIEnv *env, jclass cls,
    jint numReadGroups,
    jint numEvents,
    jint numCovariates,
    jintArray covariatesDimensions,
    jbyte LOW_QUAL_TAIL,
    jint MISMATCHES_CONTEXT_SIZE,
    jint INDELS_CONTEXT_SIZE,
    jint MAXIMUM_CYCLE_VALUE,
    jint CUSHION_FOR_INDEL) {

  if (license_verify() != 0) {
    throwAccError(env, "license check failed");
    return;
  }

  // init
  int* dims = (int*) env->GetIntArrayElements(covariatesDimensions, NULL);

  // setup global variables
  g_numReadGroups = numReadGroups;
  g_numEvents = numEvents;
  g_numCovariates = numCovariates;

  // Here need to make sure the init function is only called once.
  if (table || cov || baq) {
    DLOG(ERROR) << "static objects already initialized before init()";
    throwAccError(env, "Falcon internal error in initNative()");
    return;
    //throw std::runtime_error("Currently the init function cannot be called twice");
  }

  cov = new Covariates(numEvents, numCovariates,
                MISMATCHES_CONTEXT_SIZE,
                INDELS_CONTEXT_SIZE,
                LOW_QUAL_TAIL,
                MAXIMUM_CYCLE_VALUE,
                CUSHION_FOR_INDEL);

  DLOG(INFO) << "Initialize for BQSR";
  table = new RecalibrationTable(numReadGroups,
                numEvents, numCovariates,
                dims);

  baq = new BAQ();

  env->ReleaseIntArrayElements(covariatesDimensions, dims, 0);
}

// init function for PrintReads
//JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_initNative__I_3Lorg_broadinstitute_gatk_engine_recalibration_covariates_Covariate_2_3B_3BZIDZBIIII(
JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_initNative__I_3Lorg_broadinstitute_hellbender_utils_recalibration_covariates_Covariate_2_3B_3BZIDZBIIII(

    JNIEnv *env, jclass cls,
    jint numEvents,
    jobjectArray jcovariates,
    jbyteArray jquantizationTable,
    jbyteArray jstaticQuantizedMapping,
    jboolean disableIndelQuals,
    jint     preserveQLessThan,
    jdouble  globalQScorePrior,
    jboolean emitOriginalQuals,
    jbyte LOW_QUAL_TAIL,
    jint MISMATCHES_CONTEXT_SIZE,
    jint INDELS_CONTEXT_SIZE,
    jint MAXIMUM_CYCLE_VALUE,
    jint CUSHION_FOR_INDEL) {

  if (license_verify() != 0) {
    throwAccError(env, "license check failed");
    return;
  }

  if (!jcovariates) {
    throwAccError(env, "Falcon internal error in initNative()");
    return;
    //throw std::runtime_error("Covariates input is empty");
  }

  if (table || cov) {
    DLOG(ERROR) << "static objects already initialized before init()";
    throwAccError(env, "Falcon internal error in initNative()");
    return;
    //throw std::runtime_error("Currently the init function cannot be called twice");
  }

  int numCovariates = env->GetArrayLength(jcovariates);
  int numReadGroups = 0;

  cov = new Covariates(numEvents, numCovariates,
                MISMATCHES_CONTEXT_SIZE,
                INDELS_CONTEXT_SIZE,
                LOW_QUAL_TAIL,
                MAXIMUM_CYCLE_VALUE,
                CUSHION_FOR_INDEL);

  int* dims = (int*)malloc(numCovariates*sizeof(int));
  // initialize covariates
  for (int i = 0; i < numCovariates; i++) {
    jobject jcov = (jobject)env->GetObjectArrayElement(jcovariates, i);
    jclass cls = env->GetObjectClass(jcov);
    jmethodID maxValue_method = env->GetMethodID(cls, "maximumKeyValue", "()I");

    if (!maxValue_method) {
      DLOG(ERROR) << "Problem getting maximumKeyValue() method for covariate #" << i;
      throwAccError(env, "Falcon internal error in initNative()");
      return;
    }

    dims[i] = env->CallIntMethod(jcov, maxValue_method) + 1;
    //DLOG(INFO) << "Covariate #" << i << " has " << dims[i] << " values";
    if (i == 0) {
      // this is the RG covariate
      numReadGroups = dims[i];

      // initialize the existing readgroups
      for (int j = 0; j < numReadGroups; j++) {
        jmethodID method = env->GetMethodID(cls, "formatKey", "(I)Ljava/lang/String;");
        if (!method) {
          DLOG(ERROR) << "Problem getting formatKey() method for covariate #" << i;
          throwAccError(env, "Falcon internal error in initNative()");
          return;
        }
        jstring rg_name = (jstring)env->CallObjectMethod(jcov, method, j);
        const char* readGroup = env->GetStringUTFChars(rg_name, NULL);

        std::string rg(readGroup);
        cov->initReadGroup(rg);

        //DLOG(INFO) << "read group #" << j << " is " << readGroup;

        env->ReleaseStringUTFChars(rg_name, readGroup);
      }
    }
  }
  DLOG(INFO) << "numReadGroups = " << numReadGroups;

  // setup global variables
  g_numReadGroups = numReadGroups;
  g_numEvents = numEvents;
  g_numCovariates = numCovariates;

  DLOG(INFO) << "Initialize for PR";

  int qtable_size = 0;
  int sqmap_size  = 0;
  int8_t* quantizationTable = NULL;
  int8_t* staticQuantizedMapping = NULL;

  if (jquantizationTable) {
    qtable_size = env->GetArrayLength(jquantizationTable);
    quantizationTable = (int8_t*)env->GetByteArrayElements(jquantizationTable, NULL);
  }
  if (jstaticQuantizedMapping) {
    sqmap_size = env->GetArrayLength(jstaticQuantizedMapping);
    staticQuantizedMapping = (int8_t*)env->GetByteArrayElements(jstaticQuantizedMapping, NULL);
  }
  DLOG(INFO) << "quantizedTable = " << qtable_size;
  DLOG(INFO) << "staticQuantizedMapping = " << sqmap_size;

  table = new RecalibrationTable(numReadGroups,
                numEvents, numCovariates,
                dims,
                qtable_size,
                quantizationTable,
                sqmap_size,
                staticQuantizedMapping,
                disableIndelQuals,
                preserveQLessThan,
                globalQScorePrior,
                emitOriginalQuals);

  free(dims);
  //env->ReleaseIntArrayElements(covariatesDimensions, dims, 0);
  if (qtable_size) env->ReleaseByteArrayElements(jquantizationTable, quantizationTable, 0);
  if (sqmap_size) env->ReleaseByteArrayElements(jstaticQuantizedMapping, staticQuantizedMapping, 0);
}

JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_bqsrInitReadGroupNative(
    JNIEnv *env, jobject,
    jstring readGroupId){
    const char* cname;
    cname = (env)-> GetStringUTFChars(readGroupId, 0);
    std::string str(cname);
    //std::cout << "in Native C"<< str<<std::endl;
    cov->addReadGroup(str);
    env->ReleaseStringUTFChars(readGroupId, cname);

}

JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_setTableNative(
    JNIEnv *env, jobject obj,
    jlong     numOccurance,
    jdouble   numMismatches,
    jdouble   estimatedQReported,
    jintArray jkeys,
    jint      cov_idx) {

  int num_keys = env->GetArrayLength(jkeys);
  // verify if the dimensions of keys are correct
  if (cov_idx == 0 && num_keys != 2) {
    DLOG(ERROR) << "input keys are not compatible with RG";
    return;
  }
  else if (cov_idx == 1 && num_keys != 3) {
    DLOG(ERROR) << "input keys are not compatible with QUAL";
    return;
  }
  else if (cov_idx >= 2 && num_keys != 4) {
    DLOG(ERROR) << "input keys are not compatible with other covariates";
    return;
  }

  RecalDatum e;
  e.numOccurance = numOccurance;
  e.numMismatches = numMismatches;
  e.estimatedQReported = estimatedQReported;
  e.empiricalQuality = -1.0;

  if (!table) {
    DLOG(ERROR) << "RecalibrationTable is uninitialized";
    throwAccError(env, "Falcon internal error in setTableNative()");
    return;
    //throw std::runtime_error("unexpected error");
  }

  int* keys = env->GetIntArrayElements(jkeys, 0);

  table->put(e, keys, cov_idx);

  env->ReleaseIntArrayElements(jkeys, keys, 0);
}

JNIEXPORT jintArray JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_computeContextCovariatesNative(
    JNIEnv *env, jobject obj,
    jbyteArray jbases,
    jbyteArray jquals,
    jboolean isNegativeStrand) {

  int readLength = env->GetArrayLength(jbases);
  int8_t* bases = (int8_t*)env->GetByteArrayElements(jbases, 0);
  int8_t* quals = (int8_t*)env->GetByteArrayElements(jquals, 0);

  jintArray ret = (jintArray)env->NewIntArray(g_numEvents*readLength*g_numCovariates);
  int* keys = env->GetIntArrayElements(ret, 0);

  cov->computeContextCovariates(keys, readLength, isNegativeStrand, bases, quals);

  env->ReleaseByteArrayElements(jbases, bases, 0);
  env->ReleaseByteArrayElements(jquals, quals, 0);
  env->ReleaseIntArrayElements(ret, keys, 0);

  return ret;
}

JNIEXPORT jintArray JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_computeCycleCovariatesNative(
    JNIEnv *env, jobject obj,
    jint readLength,
    jint platformType,
    jboolean isNegativeStrand,
    jboolean isReadPaired,
    jboolean isSecondOfPair) {

  jintArray ret = (jintArray)env->NewIntArray(g_numEvents*readLength*g_numCovariates);
  int* keys = env->GetIntArrayElements(ret, 0);

  try {
    cov->computeCycleCovariates(keys, readLength, platformType,
        isNegativeStrand, isReadPaired, isSecondOfPair);
  }
  catch (std::runtime_error &e) {
    throwAccError(env, e.what());
  }

  env->ReleaseIntArrayElements(ret, keys, 0);
  return ret;
}

JNIEXPORT jintArray JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_computeCovariatesNative(
    JNIEnv *env, jobject obj,
    jbyteArray jbases,
    jbyteArray jbaseQuals,
    jbyteArray jinsertionQuals,
    jbyteArray jdeletionQuals,
    jstring  jreadGroup,
    jboolean isNegativeStrand,
    jboolean isReadPaired,
    jboolean isSecondOfPair,
    jint platformType)
{

  int readLength = env->GetArrayLength(jbases);
  int8_t* bases = (int8_t*)env->GetByteArrayElements(jbases, 0);
  int8_t* baseQuals = (int8_t*)env->GetByteArrayElements(jbaseQuals, 0);
  int8_t* insertionQuals = (int8_t*)env->GetByteArrayElements(jinsertionQuals, 0);
  int8_t* deletionQuals = (int8_t*)env->GetByteArrayElements(jdeletionQuals, 0);
  const char* readGroup = env->GetStringUTFChars(jreadGroup, NULL);

  jintArray ret = (jintArray)env->NewIntArray(g_numEvents*readLength*g_numCovariates);
  int* keys = env->GetIntArrayElements(ret, 0);

  try {
    cov->compute(keys, readLength, std::string(readGroup),
          bases, baseQuals, insertionQuals, deletionQuals,
          platformType,
          isNegativeStrand, isReadPaired, isSecondOfPair);
  }
  catch (std::runtime_error &e) {
    throwAccError(env, e.what());
  }

  env->ReleaseByteArrayElements(jbases, bases, 0);
  env->ReleaseByteArrayElements(jbaseQuals, baseQuals, 0);
  env->ReleaseByteArrayElements(jinsertionQuals, insertionQuals, 0);
  env->ReleaseByteArrayElements(jdeletionQuals, deletionQuals, 0);
  env->ReleaseStringUTFChars(jreadGroup, readGroup);
  env->ReleaseIntArrayElements(ret, keys, 0);
  return ret;
}

JNIEXPORT jbyteArray JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_calculateBAQArrayNative(
    JNIEnv *env, jobject obj,
    jbyteArray jrefBases,
    jbyteArray jbases,
    jbyteArray jbaseQuals,
    jbyteArray jcigarOps,
    jintArray  jcigarLens,
    jint refOffset) {

  uint64_t start_ns = getNs();
  uint64_t elapsed_ns = 0;

  int refLength = env->GetArrayLength(jrefBases);
  int readLength = env->GetArrayLength(jbases);
  int numCigar = env->GetArrayLength(jcigarOps);
  int8_t* refBases = (int8_t*)env->GetByteArrayElements(jrefBases, 0);
  int8_t* bases = (int8_t*)env->GetByteArrayElements(jbases, 0);
  int8_t* baseQuals = (int8_t*)env->GetByteArrayElements(jbaseQuals, 0);
  int8_t* cigarOps = (int8_t*)env->GetByteArrayElements(jcigarOps, 0);
  int* cigarLens = (int*)env->GetIntArrayElements(jcigarLens, 0);

  uint64_t start_sec_ns = getNs();
  int8_t* bqTag = (int8_t*)malloc(readLength*sizeof(int8_t));

  int ret = baq->calculateBAQ(bqTag,
          bases, baseQuals, refBases,
          cigarOps, cigarLens,
          refLength, readLength, numCigar, refOffset);
  elapsed_ns = getNs() - start_sec_ns;

  // release arrays
  env->ReleaseByteArrayElements(jrefBases, refBases, 0);
  env->ReleaseByteArrayElements(jbases, bases, 0);
  env->ReleaseByteArrayElements(jbaseQuals, baseQuals, 0);
  env->ReleaseByteArrayElements(jcigarOps, cigarOps, 0);
  env->ReleaseIntArrayElements(jcigarLens, cigarLens, 0);

  if (ret) {
    free(bqTag);
    return NULL;
  }
  jbyteArray jbqTag = (jbyteArray)env->NewByteArray(readLength);
  env->SetByteArrayRegion(jbqTag, 0, readLength, bqTag);
  // SetIntArrayRegion should release bqTag

  total_baq_time += getNs() - start_ns;
  total_baq_compute_time += elapsed_ns;

  return jbqTag;
}

JNIEXPORT jdoubleArray JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_calculateErrorsNative(
    JNIEnv *env, jobject obj,
    jbyteArray jbases,
    jbyteArray jquals,
    jbyteArray jrefForBAQ,
    jbyteArray jrefBases,
    jbyteArray jcigarOps,
    jintArray  jcigarLens,
    jbyteArray jreadBAQArray,
    jboolean   isExcludeFromBAQ,
    jboolean   isNegativeStrand,
    jint       refOffset)
{
  int readLength = env->GetArrayLength(jbases);
  int numCigar = env->GetArrayLength(jcigarOps);
  int refLength = -1;
  int8_t* refForBAQ = NULL;

  if (jrefForBAQ) { // if read is excluded from BAQ, the ref seq can be null
    refLength = env->GetArrayLength(jrefForBAQ);
    refForBAQ = (int8_t*)env->GetByteArrayElements(jrefForBAQ, 0);
  }
  int8_t* bases = (int8_t*)env->GetByteArrayElements(jbases, 0);
  int8_t* quals = (int8_t*)env->GetByteArrayElements(jquals, 0);
  int8_t* refBases = (int8_t*)env->GetByteArrayElements(jrefBases, 0);
  int8_t* cigarOps = (int8_t*)env->GetByteArrayElements(jcigarOps, 0);
  int* cigarLens = (int*)env->GetIntArrayElements(jcigarLens, 0);

  jdoubleArray ret = (jdoubleArray)env->NewDoubleArray(readLength * 3);
  double* errors = env->GetDoubleArrayElements(ret, 0);

  // NOTE: if baq.excludeReadFromBAQ, and read.getBAQTag != NULL
  // the readBAQArray passed to calculateErrors is not null
  int8_t* readBAQArray = NULL;
  if (jreadBAQArray != NULL) {
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "This read is excluded from BAQ";
    readBAQArray = (int8_t*)env->GetByteArrayElements(jreadBAQArray, 0);
  }

  int result = baq->calculateErrors(readLength, refLength, refOffset,
          bases, quals, refForBAQ, refBases,
          numCigar, cigarOps, cigarLens,
          isNegativeStrand, isExcludeFromBAQ,
          readBAQArray,
          errors, errors+readLength, errors+2*readLength);

  // release arrays
  if (readBAQArray) {
    env->ReleaseByteArrayElements(jreadBAQArray, readBAQArray, 0);
  }

  if (refForBAQ) {
    env->ReleaseByteArrayElements(jrefForBAQ, refForBAQ, 0);
  }
  env->ReleaseByteArrayElements(jbases, bases, 0);
  env->ReleaseByteArrayElements(jquals, quals, 0);
  env->ReleaseByteArrayElements(jrefBases, refBases, 0);
  env->ReleaseByteArrayElements(jcigarOps, cigarOps, 0);
  env->ReleaseIntArrayElements(jcigarLens, cigarLens, 0);
  env->ReleaseDoubleArrayElements(ret, errors, 0);

  if (result) {
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Return NULL from calculateErrors()";
    return NULL;
  }
  else {
    return ret;
  }
}

// private native int updateTableNative
JNIEXPORT int JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_updateTableNative(
    JNIEnv     *env,
    jobject    obj,
    jbyteArray jrefForBAQ,
    jbyteArray jrefBases,
    jbyteArray jbases,
    jbyteArray jbaseQuals,
    jbyteArray jinsertionQuals,
    jbyteArray jdeletionQuals,
    jbyteArray jcigarOps,
    jintArray  jcigarLens,
    jbyteArray jreadBAQArray,
    jstring    jreadGroup,
    jboolean   isNegativeStrand,
    jboolean   isReadPaired,
    jboolean   isSecondOfPair,
    jboolean   isExcludeFromBAQ,
    jint       platformType,
    jint       refOffset,
    jbooleanArray jskips)
{
  uint64_t start_ns = getNs();

  int readLength = env->GetArrayLength(jbases);
  int numCigar = env->GetArrayLength(jcigarOps);
  int refLength = 0;
  int8_t* refForBAQ = NULL;

  if (jrefForBAQ) { // when read is excluded from BAQ, the ref can be null
    refLength = env->GetArrayLength(jrefForBAQ);
    refForBAQ = (int8_t*)env->GetByteArrayElements(jrefForBAQ, 0);
  }

  // JNI array arguments, need to be released
  int8_t* refBases  = (int8_t*)env->GetByteArrayElements(jrefBases, 0);

  int8_t* bases          = (int8_t*)env->GetByteArrayElements(jbases, 0);
  int8_t* baseQuals      = (int8_t*)env->GetByteArrayElements(jbaseQuals, 0);
  int8_t* insertionQuals = (int8_t*)env->GetByteArrayElements(jinsertionQuals, 0);
  int8_t* deletionQuals  = (int8_t*)env->GetByteArrayElements(jdeletionQuals, 0);

  int8_t* cigarOps  = (int8_t*)env->GetByteArrayElements(jcigarOps, 0);
  int*    cigarLens = (int*)env->GetIntArrayElements(jcigarLens, 0);

  const char* readGroup = env->GetStringUTFChars(jreadGroup, NULL);

  // skip array used in update()
  uint8_t* skips = (uint8_t*)env->GetBooleanArrayElements(jskips, 0);

  // NOTE: if baq.excludeReadFromBAQ, and read.getBAQTag != NULL
  // the readBAQArray passed to calculateErrors is not null
  int8_t* readBAQArray = NULL;
  if (jreadBAQArray != NULL) {
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "This read is excluded from BAQ";
    readBAQArray = (int8_t*)env->GetByteArrayElements(jreadBAQArray, 0);
  }

  // first, compute the baq and errors
  // results for the error arrays (snp, indel, deletion)
  uint64_t start_sec_ns = getNs();
  double** isErrors = (double**)malloc(g_numEvents*sizeof(double*));
  isErrors[0] = (double*)malloc(readLength*sizeof(double));
  isErrors[1] = (double*)malloc(readLength*sizeof(double));
  isErrors[2] = (double*)malloc(readLength*sizeof(double));

  int result = baq->calculateErrors(readLength, refLength, refOffset,
          bases, baseQuals, refForBAQ, refBases,
          numCigar, cigarOps, cigarLens,
          isNegativeStrand, isExcludeFromBAQ,
          readBAQArray,
          isErrors[0], isErrors[1], isErrors[2]);

  total_update_baq_time += getNs() - start_sec_ns;

  // release arrays for baq and errors computation
  if (readBAQArray) {
    env->ReleaseByteArrayElements(jreadBAQArray, readBAQArray, 0);
  }
  if (refForBAQ) {
    env->ReleaseByteArrayElements(jrefForBAQ, refForBAQ, 0);
  }
  env->ReleaseByteArrayElements(jrefBases, refBases, 0);
  env->ReleaseByteArrayElements(jcigarOps, cigarOps, 0);
  env->ReleaseIntArrayElements(jcigarLens, cigarLens, 0);

  if (result) {
    // this means a null BAQ array is a result,
    // and the rest of the computation needs to be skipped
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "baqArray is null, skip update table";

    // release all the JNI arrays
    env->ReleaseByteArrayElements(jbases, bases, 0);
    env->ReleaseByteArrayElements(jbaseQuals, baseQuals, 0);
    env->ReleaseByteArrayElements(jinsertionQuals, insertionQuals, 0);
    env->ReleaseByteArrayElements(jdeletionQuals, deletionQuals, 0);
    env->ReleaseStringUTFChars(jreadGroup, readGroup);
    env->ReleaseBooleanArrayElements(jskips, skips, 0);

    // release native pointers
    free(isErrors[0]);
    free(isErrors[1]);
    free(isErrors[2]);
    free(isErrors);

    uint64_t total_elapsed_ns = getNs() - start_ns;
    total_update_total_time += total_elapsed_ns;

    // following the application convention letting it
    // know that we didn't update the table
    return 0;
  }

  // second, compute the covariates
  start_sec_ns = getNs();
  int* keys = (int*)malloc(g_numEvents*readLength*g_numCovariates*sizeof(int));

  try {
    cov->compute(keys, readLength, std::string(readGroup),
          bases, baseQuals, insertionQuals, deletionQuals,
          platformType,
          isNegativeStrand, isReadPaired, isSecondOfPair);
  } catch (std::runtime_error &e) {
    throwAccError(env, e.what());
    return 0;
  }

  total_update_covariate_time += getNs() - start_sec_ns;

  env->ReleaseByteArrayElements(jbases, bases, 0);
  env->ReleaseByteArrayElements(jbaseQuals, baseQuals, 0);
  env->ReleaseByteArrayElements(jinsertionQuals, insertionQuals, 0);
  env->ReleaseByteArrayElements(jdeletionQuals, deletionQuals, 0);
  env->ReleaseStringUTFChars(jreadGroup, readGroup);

  start_sec_ns = getNs();

  // finally, perform the update of recal tables
  table->update(readLength, keys, skips, isErrors);

  total_update_compute_time += getNs() - start_sec_ns;

  // free JNI array finally
  env->ReleaseBooleanArrayElements(jskips, skips, 0);

  // free temp pointers
  free(keys);
  free(isErrors[0]);
  free(isErrors[1]);
  free(isErrors[2]);
  free(isErrors);

  uint64_t total_elapsed_ns = getNs() - start_ns;
  total_update_total_time += total_elapsed_ns;
  total_update_num_calls += 1;

  return 1;
}

// private native byte[][] recalibrateNative
JNIEXPORT jobjectArray JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_recalibrateNative(
    JNIEnv *env, jobject obj,
    jbyteArray jbases,
    jbyteArray jbaseQuals,
    jbyteArray jinsertionQuals,
    jbyteArray jdeletionQuals,
    jstring  jreadGroup,
    jboolean isNegativeStrand,
    jboolean isReadPaired,
    jboolean isSecondOfPair,
    jint platformType)
{
  //{

    //PLACE_TIMER;
  //}
  timer.start(0);
  int readLength = env->GetArrayLength(jbases);
  int8_t* bases = (int8_t*)env->GetByteArrayElements(jbases, 0);
  int8_t** quals = (int8_t**)malloc(g_numEvents*sizeof(int8_t*));
  quals[0] = (int8_t*)env->GetByteArrayElements(jbaseQuals, 0);
  quals[1] = (int8_t*)env->GetByteArrayElements(jinsertionQuals, 0);
  quals[2] = (int8_t*)env->GetByteArrayElements(jdeletionQuals, 0);

  const char* readGroup = env->GetStringUTFChars(jreadGroup, NULL);

  timer.start(1);
  // compute covariates
  int* keys = (int*)malloc(g_numEvents*readLength*g_numCovariates*sizeof(int));
  try {
    cov->compute(keys, readLength, std::string(readGroup),
          bases, quals[0], quals[1], quals[2],
          platformType,
          isNegativeStrand, isReadPaired, isSecondOfPair);
  } catch (std::runtime_error &e) {
    throwAccError(env, e.what());
    return 0;
  }

  timer.stop(1);

  // recalibrate
  int8_t** recal_quals = (int8_t**)malloc(g_numEvents*sizeof(int8_t*));
  for (int i = 0; i < g_numEvents; i++) {
    // the recal table should not need to be freed since it's
    // returned to java
    recal_quals[i] = (int8_t*)malloc(readLength*sizeof(int8_t));
  }

  timer.start(2);
  table->recalibrate(readLength, keys, quals, recal_quals);
  timer.stop(2);

  // return value is byte[][] --> B]]
  jclass byte_2d_cls = env->FindClass("[B");
  if (!byte_2d_cls) DLOG(ERROR) << "cannot find byte array class";

  env->ReleaseByteArrayElements(jbases, bases, 0);
  env->ReleaseByteArrayElements(jbaseQuals, quals[0], 0);
  env->ReleaseByteArrayElements(jinsertionQuals, quals[1], 0);
  env->ReleaseByteArrayElements(jdeletionQuals, quals[2], 0);
  env->ReleaseStringUTFChars(jreadGroup, readGroup);
  free(quals);
  free(keys);

  jobjectArray ret = env->NewObjectArray(g_numEvents, byte_2d_cls, NULL);

  for (int i = 0; i < g_numEvents; i++) {
    jbyteArray array = env->NewByteArray(readLength);
    env->SetByteArrayRegion(array, 0, readLength, recal_quals[i]);
    env->SetObjectArrayElement(ret, i, array);
    env->ReleaseByteArrayElements(array, recal_quals[i], 0);
  }

  free(recal_quals);
  timer.stop(0);

  //DLOG_EVERY_N(INFO, TIMER_SAMPLE_RATE) << "Time for " << timer.get_laps(0)
  //           << " calls on compute covariate: "
  //           << timer.get_ms(1) << " ms";
  //DLOG_EVERY_N(INFO, TIMER_SAMPLE_RATE) << "Time for " << timer.get_laps(0)
  //          << " calls on recalibrate: "
  //          << timer.get_ms(2) << " ms";
  //DLOG_EVERY_N(INFO, TIMER_SAMPLE_RATE) << "Time for " << timer.get_laps(0)
  //          << " calls on jni communication: "
  //          << timer.get_ms(0) - timer.get_ms(1) - timer.get_ms(2) << " ms";

  return ret;
}

JNIEXPORT jobjectArray JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_getTableNative(
    JNIEnv* env, jobject obj) {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Return values from the RecalibrationTables";

  // TODO: need to double check on these signatures
  jclass ret_cls = env->FindClass("com/falconcomputing/genomics/bqsr/FalconRecalibrationEngine$RecalDatumTable");
  if (!ret_cls) {
    DLOG(ERROR) << "Cannot find class for RecalDatumTable";
    return NULL;
  }
  jmethodID ret_init = env->GetMethodID(ret_cls, "<init>",
                        "(Lcom/falconcomputing/genomics/bqsr/FalconRecalibrationEngine;II)V");
  if (!ret_init) {
    DLOG(ERROR) << "Cannot find constructor for RecalDatumTable";
    return NULL;
  }

  // first create the object for return
  jobjectArray ret = (jobjectArray)env->NewObjectArray(g_numCovariates - 1, ret_cls, NULL);

  // start filling the return array with object instances
  for (int i = 1; i < g_numCovariates; i++) {
    std::vector<int> table_dims = table->getTableDimensions(i);
    int table_size = table->getTableSize(i);
    int num_dimensions = table_dims.size();

    // construct one table in the return array
    jobject ret_obj = env->NewObject(ret_cls, ret_init, obj, table_size, num_dimensions);
    if (!ret_obj) {
      DLOG(ERROR) << "Cannot construct object of RecalDatumTable";
      return NULL;
    }

    env->SetObjectArrayElement(ret, i-1, ret_obj);

    // get reference for each member fields
    jlongArray numOccur = (jlongArray)env->GetObjectField(ret_obj,
                             env->GetFieldID(ret_cls, "numOccurance", "[J"));
    jdoubleArray numError = (jdoubleArray)env->GetObjectField(ret_obj,
                                  env->GetFieldID(ret_cls, "numMismatches", "[D"));
    jintArray tableDimensions = (jintArray)env->GetObjectField(ret_obj,
                                  env->GetFieldID(ret_cls, "tableDimensions", "[I"));

    if (!numOccur || !numError || !tableDimensions) {
      DLOG(ERROR) << "Cannot reference to members of RecalDatumTable";
      return NULL;
    }

    // copy over the results
    jlong* num_occur = env->GetLongArrayElements(numOccur, NULL);
    jdouble* num_error = env->GetDoubleArrayElements(numError, NULL);
    jint* table_dimensions = env->GetIntArrayElements(tableDimensions, NULL);

    if (table->getTable(i)) { // DatumTable can be empty
      Datum* tableContents = table->getTable(i);
      for (int k = 0; k < table_size; k++) {
        num_occur[k] = tableContents[k].numOccurance;
        num_error[k] = tableContents[k].numMismatches;
      }
    }
    else if (table->getFullTable(i)) {
      RecalDatum* tableContents = table->getFullTable(i);
      for (int k = 0; k < table_size; k++) {
        num_occur[k] = tableContents[k].numOccurance;
        num_error[k] = tableContents[k].numMismatches;
      }
    }
    else {
      DLOG(ERROR) << "recalibration table is not initialized";
    }

    for (size_t k = 0; k < table_dims.size(); k++) {
      table_dimensions[k] = table_dims[k];
      //DLOG(INFO) << "Table #" << i << ", "
      //           << "dim-" << k << " = " << table_dimensions[k];
    }

    env->ReleaseLongArrayElements(numOccur, num_occur, 0);
    env->ReleaseDoubleArrayElements(numError, num_error, 0);
    env->ReleaseIntArrayElements(tableDimensions, table_dimensions, 0);
  }
  return ret;
}

JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_updateReadGroupCovariatesNative(
    JNIEnv *env, jobject obj, jobject readGroupCovariates)
{
  jclass cls = env->GetObjectClass(readGroupCovariates);
  jmethodID method = env->GetMethodID(cls, "keyForReadGroup", "(Ljava/lang/String;)I");
  if (!method) {
    DLOG(ERROR) << "Failed to get the method";
    return;
  }

  for (auto rg : cov->getRGKeys()) {
    jstring jreadGroup = env->NewStringUTF(rg.c_str());
    jint idx = env->CallIntMethod(readGroupCovariates, method, jreadGroup);
    DLOG(INFO) << "Set rg: " << rg << " with idx = " << idx;
  }
}

JNIEXPORT jdouble JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_log10QempPriorNative(
    JNIEnv *env, jobject obj,
    jdouble qemp, jdouble qreported) {

  MathUtils util;
  return util.log10QempPrior(qemp, qreported);
}

JNIEXPORT jdouble JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_log10QempLikelihoodNative(
    JNIEnv *env, jobject obj,
    jdouble qemp, jlong numObservations, jlong numMismatches) {

  MathUtils util;
  return util.log10QempLikelihood(qemp, numObservations, numMismatches);
}

JNIEXPORT jdouble JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_bayesianEstimateOfEmpiricalQualityNative(
    JNIEnv *env, jobject obj,
    jlong numObservations, jlong numMismatches, jdouble qreported) {

  MathUtils util;
  return util.bayesianEstimateOfEmpiricalQuality(numObservations, numMismatches, qreported);
}

JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_bqsr_FalconRecalibrationEngine_finalizeNative(
    JNIEnv *env, jobject obj) {

  if (total_update_num_calls > 0) {
    uint64_t total_update_jni_time = total_update_total_time -
    total_update_covariate_time -
    total_update_compute_time -
    total_update_baq_time;

    // output performance counter
    DLOG(INFO) << "Total number of calls to update(): " << total_update_num_calls;
    DLOG(INFO) << "Total time on BAQ compute: " << (double)total_update_baq_time / 1e6 << " ms";
    DLOG(INFO) << "Total time on compute covariates: " << (double)total_update_covariate_time / 1e6 << " ms";
    DLOG(INFO) << "Total time on table update: " << (double)total_update_compute_time / 1e6 << " ms";
    DLOG(INFO) << "Total time on jni communication: " << (double)total_update_jni_time / 1e6 << " ms";
  }

  // free up resource
  if (baq) delete baq;
  if (cov) delete cov;
  if (table) delete table;

  baq = NULL;
  cov = NULL;
  table = NULL;
}
