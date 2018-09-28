#include <vector>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "com_falconcomputing_genomics_haplotypecaller_FalconPairhmm.h"
#include "gkl-pairhmm/Context.h"
#include "gkl-pairhmm/avx_impl.h"
#include "ksight/tools.h"
#include "pairhmm/client/PairhmmClient.h"
#include "pairhmm/client/PairhmmWorker.h"
#include "pairhmm/common/PairhmmJavaData.h"
#include <xmmintrin.h>

JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_haplotypecaller_FalconPairhmm_initNative
(JNIEnv* env, jclass cls, jclass readDataHolder, jclass haplotypeDataHolder, jboolean use_double, jint max_threads, jboolean use_fpga){

    PairhmmContext* context = get_pairhmm_context();
    PairHMMClient* client = get_pairhmm_client();
    JavaData javaData;
    javaData.init(env, readDataHolder, haplotypeDataHolder);
    context->g_use_double = use_double;

    // enable FTZ
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

    // init convert char table
    ConvertChar::init();
    for(int i = 0; i < 3; i++)
        inputBank[i].first = "bankID";
#ifdef DEPLOY_aws
    inputBank[0].second = 3; inputBank[1].second = 1; inputBank[2].second = 0;
#else
    inputBank[0].second = 0; inputBank[1].second = 1; inputBank[2].second = 3;
#endif
}


JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_haplotypecaller_FalconPairhmm_computeLikelihoodsNative
(JNIEnv* env, jobject obj, jobjectArray readDataArray, jobjectArray haplotypeDataArray, jdoubleArray likelihoodArray){
    PLACE_TIMER;
    //==================================================================
    // get Java data
    PairhmmContext* context = get_pairhmm_context();
    int numRead = env->GetArrayLength(readDataArray);
    int numHap = env->GetArrayLength(haplotypeDataArray);
    bool useFPGA = !context->g_use_double;
    JavaData javaData;
    javaData.getData(env, readDataArray, haplotypeDataArray, context->reads, context->haps, numRead, numHap, useFPGA);
    context->likelihoods = javaData.getOutputArray(env, likelihoodArray);

    //==================================================================
    // call PairHMM worker
    context->pmm_count++;
    PairHMMClient* client = get_pairhmm_client();
    PairHMMWorker worker(client, numRead, numHap, &context->reads[0], &context->haps[0]);
    worker.run();
    
    //==================================================================
    // copy output to Java data
    worker.getOutput(context->likelihoods);

    //==================================================================
    // release Java data
    javaData.releaseData(env);
}



JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_haplotypecaller_FalconPairhmm_doneNative
(JNIEnv* env, jobject obj)
{
#ifndef NO_PROFILE
    ksight::ksight.print_total(); 
#endif

    release_pairhmm_context();
    release_pairhmm_client();
}
