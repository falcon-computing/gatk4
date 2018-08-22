#ifndef PAIRHMMJAVADATA_H
#define PAIRHMMJAVADATA_H

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pairhmm/common/pairhmm_common.h"

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

    void getData(JNIEnv *env, jobjectArray& readDataArray, jobjectArray& haplotypeDataArray, pairhmmInput* ret, int numRead, int numHap, bool& useFPGA){
        float cells = 0.0;
        ret->reads.clear();
        ret->haps.clear();
        // get haplotypes
        for (int i = 0; i < numHap; i++) {
            int curHapLen;
            char* haps = getCharArray(env, haplotypeDataArray, i, m_haplotypeBasesFid, curHapLen);
            Hap curHap;
            for(int j = 0; j < curHapLen; j++){
                curHap.bases.push_back(haps[j]);
            }
            ret->haps.push_back(curHap);
            if(curHapLen > MAX_HAP_LEN)
                useFPGA = false;
        }

        // get reads 
        for (int r = 0; r < numRead; r++) {
            int curReadLen;
            char* reads = getCharArray(env, readDataArray, r, m_readBasesFid, curReadLen);
            char* insGops = getCharArray(env, readDataArray, r, m_insertionGopFid, curReadLen);
            char* delGops = getCharArray(env, readDataArray, r, m_deletionGopFid, curReadLen);
            char* gapConts = getCharArray(env, readDataArray, r, m_overallGcpFid, curReadLen);
            char* readQuals = getCharArray(env, readDataArray, r, m_readQualsFid, curReadLen);
            Read curRead;
            for(int j = 0; j < curReadLen; j++){
                curRead.bases.push_back(reads[j]);
                curRead._i.push_back(insGops[j]);
                curRead._d.push_back(delGops[j]);
                curRead._c.push_back(gapConts[j]);
                curRead._q.push_back(readQuals[j]);
            }
            ret->reads.push_back(curRead);
            if(curReadLen > MAX_READ_LEN)
                useFPGA = false;
        }
        for(int i = 0; i < numRead; i++){
            for(int j = 0; j < numHap; j++){
                cells = cells + ret->reads[i].bases.size() * ret->haps[j].bases.size();
            }
        }
        float FPGA_time = cells / FPGA_perf + DIE_NUM * 1e6 + DIE_NUM * sizeof(FPGAInput) / 3.0 + numRead * numHap * sizeof(float) / 3.0;
        float AVX_time = cells / AVX_perf;
        if(AVX_time <= FPGA_time)
            useFPGA = false;
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
