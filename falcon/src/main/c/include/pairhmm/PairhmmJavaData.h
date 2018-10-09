#ifndef PAIRHMMJAVADATA_H
#define PAIRHMMJAVADATA_H

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PROFILING

class JavaData {
 public:
  // cache field ids
  static void init(JNIEnv *env, jclass readDataHolder, jclass haplotypeDataHolder) {
    m_readBasesFid = getFieldId(env, readDataHolder, "readBases", "[B");
    m_readQualsFid = getFieldId(env, readDataHolder, "readQuals", "[B");
    m_insertionGopFid = getFieldId(env, readDataHolder, "insertionGOP", "[B");
    m_deletionGopFid = getFieldId(env, readDataHolder, "deletionGOP", "[B");
    m_overallGcpFid = getFieldId(env, readDataHolder, "overallGCP", "[B");
    m_haplotypeBasesFid = getFieldId(env, haplotypeDataHolder, "haplotypeBases", "[B");
  }

  void getData(JNIEnv *env, 
      jobjectArray& readDataArray, 
      jobjectArray& haplotypeDataArray, 
      std::vector<read_t>& reads, 
      std::vector<hap_t>& haps, 
      int numRead, int numHap) 
  {
    float cells = 0.0;
    reads.clear();
    haps.clear();
    // get haplotypes
    for (int i = 0; i < numHap; i++) {
      int curHapLen;
      char* hapBases = getCharArray(env, haplotypeDataArray, i, m_haplotypeBasesFid, curHapLen);
      hap_t curHap;
      curHap.len = curHapLen;
      curHap._b  = hapBases;
      haps.push_back(curHap);
    }

    // get reads 
    for (int r = 0; r < numRead; r++) {
      int curReadLen;
      char* readBases = getCharArray(env, readDataArray, r, m_readBasesFid, curReadLen);
      char* insGops = getCharArray(env, readDataArray, r, m_insertionGopFid, curReadLen);
      char* delGops = getCharArray(env, readDataArray, r, m_deletionGopFid, curReadLen);
      char* gapConts = getCharArray(env, readDataArray, r, m_overallGcpFid, curReadLen);
      char* readQuals = getCharArray(env, readDataArray, r, m_readQualsFid, curReadLen);
      read_t curRead;
      curRead.len = curReadLen;
      curRead._b  = readBases;
      curRead._i  = insGops;
      curRead._d  = delGops;
      curRead._c  = gapConts;
      curRead._q  = readQuals;
      reads.push_back(curRead);
    }
    for(int i = 0; i < numRead; i++){
      for(int j = 0; j < numHap; j++){
        cells = cells + reads[i].len * haps[j].len;
      }
    }
  }

  double* getOutputArray(JNIEnv *env, jdoubleArray array) {
    return getDoubleArray(env, array);
  }

  void releaseData(JNIEnv *env) {
    for (int i = 0; (size_t)i < m_byteArrays.size(); i++) {
      env->ReleaseByteArrayElements(m_byteArrays[i].first, m_byteArrays[i].second, 0);
    }
    m_byteArrays.clear();
    for (int i = 0; (size_t)i < m_doubleArrays.size(); i++) {
      env->ReleaseDoubleArrayElements(m_doubleArrays[i].first, m_doubleArrays[i].second, 0);
    }
    m_doubleArrays.clear();
  }

  static jfieldID getFieldId(JNIEnv *env, jclass clazz, const char *name, const char *sig) {
    jfieldID id = env->GetFieldID(clazz, name, sig);
    if (id == NULL) {
      env->ThrowNew(env->FindClass("java/lang/Exception"), "Unable to get field ID"); 
    }
    return id;
  }

  char* getCharArray(JNIEnv* env, jobjectArray& array, int index, jfieldID fieldId, int& length) {
    jobject object = env->GetObjectArrayElement(array, index);
    jbyteArray byteArray = (jbyteArray)env->GetObjectField(object, fieldId);
    jbyte* primArray = (jbyte*)env->GetByteArrayElements(byteArray, NULL);
    if (primArray == NULL) {
      env->ThrowNew(env->FindClass("java/lang/OutOfMemoryError"), "Unable to access jbyteArray");
    }
    length = env->GetArrayLength(byteArray);
    m_byteArrays.push_back(std::make_pair(byteArray, primArray));
    return (char*)primArray;
  }

  double* getDoubleArray(JNIEnv *env, jdoubleArray array) {
    jdouble* primArray = (jdouble*)env->GetDoubleArrayElements(array, NULL);
    if (primArray == NULL) {
      env->ThrowNew(env->FindClass("java/lang/OutOfMemoryError"), "Unable to access jdoubleArray");
    }
    m_doubleArrays.push_back(std::make_pair(array, primArray));
    return (double*)primArray;
  }

  std::vector<std::pair<jbyteArray, jbyte*> > m_byteArrays;
  std::vector<std::pair<jdoubleArray, jdouble*> > m_doubleArrays;

  static jfieldID m_readBasesFid;
  static jfieldID m_readQualsFid;
  static jfieldID m_insertionGopFid;
  static jfieldID m_deletionGopFid;
  static jfieldID m_overallGcpFid;
  static jfieldID m_haplotypeBasesFid;  
};

jfieldID JavaData::m_readBasesFid;
jfieldID JavaData::m_readQualsFid;
jfieldID JavaData::m_insertionGopFid;
jfieldID JavaData::m_deletionGopFid;
jfieldID JavaData::m_overallGcpFid;
jfieldID JavaData::m_haplotypeBasesFid;  

#endif
