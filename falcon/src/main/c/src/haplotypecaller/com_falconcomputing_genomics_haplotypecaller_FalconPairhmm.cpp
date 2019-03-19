#include <algorithm>
#include <assert.h>
#include <map>
#include <math.h>
#include <mutex>
#include <stdio.h>
#include <stdlib.h>
#include <thread>
#include <unordered_map>
#include <vector>
#include <xmmintrin.h>

#include "Common.h"
#include "com_falconcomputing_genomics_haplotypecaller_FalconPairhmm.h"
#include "gkl-pairhmm/avx_impl.h"
#ifndef NO_PROFILE
#include "ksight/tools.h"
#endif
#include "PairHMMClient.h"
#include "PairHMMWorker.h"
#include "PairhmmJavaData.h"

std::mutex pmm_client_map_mutex;
std::unordered_map<std::thread::id, PairHMMClient*> pmm_client_map;

PairHMMClient* get_pairhmm_client() {
    std::thread::id this_id = std::this_thread::get_id();
    PairHMMClient* client;
    pmm_client_map_mutex.lock();
    std::unordered_map<std::thread::id, PairHMMClient*>::iterator iter = pmm_client_map.find(this_id);
    if(iter == pmm_client_map.end()){
        client = new PairHMMClient();
        pmm_client_map[this_id] = client;
    }
    else{
        client = iter->second;
    }
    pmm_client_map_mutex.unlock();
    return client;
}

void release_pairhmm_client() {
    std::thread::id this_id = std::this_thread::get_id();
    pmm_client_map_mutex.lock();
    std::unordered_map<std::thread::id, PairHMMClient*>::iterator iter = pmm_client_map.find(this_id);
    if (iter != pmm_client_map.end()) {
        PairHMMClient* client = iter->second;
        delete client;
        pmm_client_map.erase(this_id);   
    }
    pmm_client_map_mutex.unlock();
}

JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_haplotypecaller_FalconPairhmm_initNative
(JNIEnv* env, jclass cls, jclass readDataHolder, jclass haplotypeDataHolder, jboolean use_double, jint max_threads, jboolean use_fpga){

    PairHMMClient* client = get_pairhmm_client();
    JavaData javaData;
    javaData.init(env, readDataHolder, haplotypeDataHolder);

    // enable FTZ
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

    // init convert char table
    ConvertChar::init();
}


JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_haplotypecaller_FalconPairhmm_computeLikelihoodsNative
(JNIEnv* env, jobject obj, jobjectArray readDataArray, jobjectArray haplotypeDataArray, jdoubleArray likelihoodArray) {
#ifndef NO_PROFILE
  PLACE_TIMER;
#endif

  // get Java data
  int num_read = env->GetArrayLength(readDataArray);
  int num_hap  = env->GetArrayLength(haplotypeDataArray);

  std::vector<read_t> reads;
  std::vector<hap_t>  haps;

  JavaData javaData;
  javaData.getData(env, readDataArray, haplotypeDataArray, reads, haps, num_read, num_hap);
  double* likelihoods = javaData.getOutputArray(env, likelihoodArray);


  // call PairHMM worker
  PairHMMClient* client = get_pairhmm_client();
  PairHMMWorker worker(client, num_read, num_hap, &reads[0], &haps[0]);
  worker.run();

  // copy output to Java data
  worker.getOutput(likelihoods);

  // release Java data
  javaData.releaseData(env);
}


JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_haplotypecaller_FalconPairhmm_doneNative(JNIEnv* env, jobject obj) {
#ifndef NO_PROFILE
  ksight::ksight.print_total(); 
#endif
  release_pairhmm_client();
}
