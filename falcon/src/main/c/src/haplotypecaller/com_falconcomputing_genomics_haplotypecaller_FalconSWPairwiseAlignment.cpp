#include <jni.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <emmintrin.h>
#include <immintrin.h>
#include <xmmintrin.h>

#include "com_falconcomputing_genomics_haplotypecaller_FalconSWPairwiseAlignment.h"
#include "Common.h"
#include "smithwaterman/avx2_impl.h"

jfieldID state_list_fd;
jfieldID length_list_fd;
jfieldID cigar_list_size_fd;
jfieldID alignment_offset_fd;
jfieldID ret_code_fd;

struct SWPairJNIObj{
    jclass results_class;
    jobject ret_object;
    jintArray length_list;
    jintArray state_list;
};


std::mutex obj_map_mutex;
std::map<std::thread::id, struct SWPairJNIObj> SWPairJNIObj_map;

struct SWPairJNIObj init_JNI_obj(JNIEnv* env, jobject& obj){
    std::thread::id this_id = std::this_thread::get_id();
    std::map<std::thread::id, struct SWPairJNIObj>::iterator JNIObj_mapIter = SWPairJNIObj_map.find(this_id);
    if(JNIObj_mapIter == SWPairJNIObj_map.end()){
        struct SWPairJNIObj objs;
        jintArray length_list_global = env->NewIntArray(MAX_SEQ_LENGTH);
        if (env->ExceptionCheck()) { 
            throwAccError(env, "failed new jintArray length_list");
        }
        objs.length_list = (jintArray) env->NewGlobalRef(length_list_global);
        jintArray state_list_global = env->NewIntArray(MAX_SEQ_LENGTH);
        if (env->ExceptionCheck()) { 
            throwAccError(env, "failed new jintArray state_list");
        }
        objs.state_list = (jintArray) env->NewGlobalRef(state_list_global);
    
        jclass results_class_global = env->FindClass("Lcom/falconcomputing/genomics/haplotypecaller/FalconSWPairwiseAlignment$results_native;");
        if(results_class_global == NULL){
            throwAccError(env, "Failed to FindClass results_native");
        }
        objs.results_class = (jclass) env->NewGlobalRef(results_class_global);
 
        jmethodID results_class_cnstrctr = env->GetMethodID(results_class_global, "<init>", "(Lcom/falconcomputing/genomics/haplotypecaller/FalconSWPairwiseAlignment;[I[IIII)V");
        if(results_class_cnstrctr == NULL){
            throwAccError(env, "Failed to Get constructor of results_native");
        }

        jobject ret_object_global = env->NewObject(results_class_global, results_class_cnstrctr, obj, state_list_global, length_list_global, 0, 0, (jint)0);
        objs.ret_object = (jclass) env->NewGlobalRef(ret_object_global);

        obj_map_mutex.lock();
        SWPairJNIObj_map[this_id] = objs;
        obj_map_mutex.unlock();
        return objs;
    }
    else{
        return JNIObj_mapIter->second;
    }
}

JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_haplotypecaller_FalconSWPairwiseAlignment_init_1native(JNIEnv *env, jobject obj){
    if (license_verify() != 0) {
        throwAccError(env, "license check failed");
        return;
    }
    struct SWPairJNIObj JNIObjs = init_JNI_obj(env, obj);
   
    state_list_fd = env->GetFieldID(JNIObjs.results_class, "state_list", "[I");
    length_list_fd = env->GetFieldID(JNIObjs.results_class, "length_list", "[I");
    cigar_list_size_fd = env->GetFieldID(JNIObjs.results_class, "cigar_list_size", "I");
    alignment_offset_fd = env->GetFieldID(JNIObjs.results_class, "alignment_offset", "I");
    ret_code_fd = env->GetFieldID(JNIObjs.results_class, "ret_code", "I");
}

JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_haplotypecaller_FalconSWPairwiseAlignment_done_1native(JNIEnv *env, jobject obj){
    /*if(results_class != NULL){
        env->NewGlobalRef(results_class);
    }
    if(Exception != NULL){
        env->NewGlobalRef(Exception);
    }*/
}

JNIEXPORT jobject JNICALL Java_com_falconcomputing_genomics_haplotypecaller_FalconSWPairwiseAlignment_align_1native(JNIEnv *env, jobject thisObj, jbyteArray reference, jbyteArray alternate, jint w_match, jint w_mismatch, jint w_open, jint w_extend, jint overhang_strategy, jboolean cutoff){

    struct SWPairJNIObj JNIObjs = init_JNI_obj(env, thisObj);
    jobject ret_object = JNIObjs.ret_object;
    jintArray length_list = JNIObjs.length_list;
    jintArray state_list = JNIObjs.state_list;

    int ret_code = 0;
    int n = (int)(env->GetArrayLength(reference)) + 1;
    int m = (int)(env->GetArrayLength(alternate)) + 1;
       
    int overhang_strategy_native = (int)overhang_strategy; //0 is SOFT_CLIP, 1 is INDEL, 2 is LEADING_INDEL, 3 is IGNORE
    if(overhang_strategy < 0 || overhang_strategy > 3){
        //fprintf(stderr, "Wrong strategy\n");
        ret_code = -1;
        throwAccError(env, "Illegal overhang_strategy");
        return NULL;
    }
   

    //now initialize the native reference and native alternate arrays
    jboolean isCopy_reference;
    jboolean isCopy_alternate;
    char* reference_native = (char*)env->GetByteArrayElements(reference, &isCopy_reference);
    if(!reference_native){
        ret_code = -1;
        throwAccError(env, "Failed to get native reference[]");
        return NULL;
    }
    
    char* alternate_native = (char*)env->GetByteArrayElements(alternate, &isCopy_alternate);
    if(!alternate_native){
        ret_code = -1;
        throwAccError(env, "Failed to get native alternate[]");
        return NULL;
    }

    int w_match_native = (int)w_match;
    int w_mismatch_native = (int)w_mismatch;
    int w_open_native = (int)w_open;
    int w_extend_native = (int)w_extend;
    int* length_list_native = NULL;
    int* state_list_native = NULL;
    int Cigar_list_size;
    int alignment_offset_native;


    SWPairwiseAlignmentOnceGKL(reference_native, n - 1, alternate_native, m - 1, w_match_native, w_mismatch_native, w_open_native, w_extend_native, overhang_strategy_native, &length_list_native, &state_list_native, Cigar_list_size, alignment_offset_native); 
    /*calculate_matrix_ret = calculateMatrixRowWiseSIMDUnroll4x(reference_native, alternate_native, sw, m, n, btrack, overhang_strategy_native, cutoff_native, w_match_native, w_mismatch_native, w_open_native, w_extend_native);

    if(calculate_matrix_ret < 0){
        return NULL;
    }

    calculate_cigar_ret = calculateCigar(env, sw, btrack, n, m, overhang_strategy_native, &length_list_native, &state_list_native, &Cigar_list_size, &alignment_offset_native);
    if(calculate_cigar_ret < 0){
        return NULL;
    }*/
        
///now covert the length_list_native and state_list_native and alignment_offset_native to java vars
    env->SetIntArrayRegion(length_list, 0, Cigar_list_size, length_list_native);
    if (env->ExceptionCheck()) { 
        throwAccError(env, "Failed to setIntArrayRegion length_list");
        return NULL;
    }

    env->SetIntArrayRegion(state_list, 0, Cigar_list_size, state_list_native);
    if (env->ExceptionCheck()) { 
        throwAccError(env, "Failed to setIntArrayRegion state_list");
        return NULL;
    }

    jint alignment_offset = (jint)alignment_offset_native;
 
    jint ret_code_java = (jint)ret_code;

    env->SetObjectField(ret_object, state_list_fd, state_list);
    env->SetObjectField(ret_object, length_list_fd, length_list);
    env->SetIntField(ret_object, cigar_list_size_fd, Cigar_list_size);
    env->SetIntField(ret_object, alignment_offset_fd, alignment_offset);
    env->SetIntField(ret_object, ret_code_fd, ret_code_java);

    if (env->ExceptionCheck()) { 
        throwAccError(env, "Failed to new return object");
        return NULL;
    }

    free(length_list_native);
    free(state_list_native);
    return ret_object;
  
}
/*
int calculateMatrix(JNIEnv *env, char* reference, char* alternate, int** sw, int ncol, int nrow, int** btrack, int overhang_strategy_native, int cutoff, int w_match, int w_mismatch, int w_open, int w_extend){
    int matrix_min_cutoff;
    if(cutoff) matrix_min_cutoff = 0;
    else matrix_min_cutoff = (int)(-1e8);
    int i = 0;
    signed int lowInitValue = -1073741824;
    int* best_gap_v;
    int* best_gap_h;
    int* gap_size_v;
    int* gap_size_h;
    if(!(best_gap_v = (int*)_mm_malloc((ncol + 1) * sizeof(int), 64))){
        throwAccError(env, "failed to malloc best_gap_v[]");
        //env->ThrowNew(Exception_SWPair, "failed to malloc best_gap_v[]");
        return -1;
    }
    if(!(best_gap_h = (int*)_mm_malloc((nrow + 1) * sizeof(int), 64))){
        throwAccError(env, "failed to malloc best_gap_h[]");
        //env->ThrowNew(Exception_SWPair, "failed to malloc best_gap_h[]");
        return -1;
    }
    if(!(gap_size_v = (int*)_mm_malloc((ncol + 1) * sizeof(int), 64))){
        throwAccError(env, "failed to malloc gap_size_v[]");
        //env->ThrowNew(Exception_SWPair, "failed to malloc gap_size_v[]");
        return -1;
    }
    if(!(gap_size_h = (int*)_mm_malloc((nrow + 1) * sizeof(int), 64))){
        throwAccError(env, "failed to malloc gap_size_h[]");
        //env->ThrowNew(Exception_SWPair, "failed to malloc gap_size_h[]");
        return -1;
    }

    for(i = 0; i < ncol + 1; ++i){
        best_gap_v[i] = lowInitValue;
        gap_size_v[i] = 0;
    }
    for(i = 0; i < nrow + 1; ++i){
        best_gap_h[i] = lowInitValue;
        gap_size_h[i] = 0;
    }
    
    if(overhang_strategy_native == OVERHANG_STRATEGY_INDEL || overhang_strategy_native == OVERHANG_STRATEGY_LEADING_INDEL){
        sw[0][1] = w_open;
        int currentValue = w_open;
        for(i = 2; i < ncol; i++){
            currentValue += w_extend;
            sw[0][i] = currentValue;
        }
        sw[1][0] = w_open;
        currentValue = w_open;
        for(i = 2; i < nrow; i++){
            currentValue += w_extend;
            sw[i][0] = currentValue;
        }
        
    }


    //build smith-waterman matrix and keep backtrack info
    int curRow_id = 0;
    int lastRow_id = 0;
    int curBackTrackRow_id = 0;
    int j = 0;
    char a_base = 0;
    char b_base = 0;
    int step_diag = 0;
    int prev_gap = 0;
    int step_down = 0;
    int kd = 0;
    int step_right = 0;
    int ki = 0;
    int diagHighestOrEqual = 0; 
    for(i = 1; i < nrow; ++i){
        a_base = reference[i - 1];
        lastRow_id = curRow_id;
        curRow_id = i;
        curBackTrackRow_id = i;
        for(j = 1; j< ncol; ++j){
            b_base = alternate[j - 1];
            
            step_diag = sw[lastRow_id][j - 1] + wd(a_base, b_base, w_match, w_mismatch);
            prev_gap = sw[lastRow_id][j] + w_open;
            best_gap_v[j] += w_extend;
            if(prev_gap > best_gap_v[j]){
                best_gap_v[j] = prev_gap;
                gap_size_v[j] = 1;
            }
            else{
                gap_size_v[j]++;
            }

            step_down = best_gap_v[j];
            kd = gap_size_v[j];
            
            prev_gap = sw[curRow_id][j - 1] + w_open;
            best_gap_h[i] += w_extend;
            if(prev_gap > best_gap_h[i]){
                best_gap_h[i] = prev_gap;
                gap_size_h[i] = 1;
            }
            else{
                gap_size_h[i]++;
            }
            step_right = best_gap_h[i];
            ki = gap_size_h[i];
            //priority here will be step diagonal, step righ, step down
            diagHighestOrEqual = (step_diag >= step_down) && (step_diag >= step_right);
            
            if(diagHighestOrEqual){
                sw[curRow_id][j] = max(matrix_min_cutoff, step_diag);
                btrack[curBackTrackRow_id][j] = 0;
            }
            else if(step_right >= step_down){
                sw[curRow_id][j] = max(matrix_min_cutoff, step_right);
                btrack[curBackTrackRow_id][j] = -ki;
            }
            else{
                sw[curRow_id][j] = max(matrix_min_cutoff, step_down);
                btrack[curBackTrackRow_id][j] = kd;
            }
        }
    }

    _mm_free(best_gap_h);
    _mm_free(best_gap_v);
    _mm_free(gap_size_h);
    _mm_free(gap_size_v);
    return 0;
}

int calculateCigar(JNIEnv *env, int** sw, int** btrack, int nrow, int ncol, int overhang_strategy, int** length_list_ptr, int** state_list_ptr, int* Cigar_list_length, int* alignment_offset_ptr){
    int p1 = 0;
    int p2 = 0;
    int refLength = nrow - 1;
    int altLength = ncol - 1;
    int maxscore = -2147483648;
    int segment_length = 0;
    int i = 0;
    int j = 0;
    int curScore = 0;
    if(overhang_strategy == OVERHANG_STRATEGY_INDEL){
        p1 = refLength;
        p2 = altLength;
    }
    else{
        p2 = altLength;

        for(i = 1; i < nrow; ++i){
            curScore = sw[i][altLength];
            if(curScore >= maxscore){
                p1 = i;
                maxscore = curScore;
            }
        }
        //now look for a larger score on the bottom-most row
        if(overhang_strategy != OVERHANG_STRATEGY_LEADING_INDEL){
            for(j = 1; j < ncol; ++j){
                curScore = sw[refLength][j];
                if(curScore > maxscore || (curScore == maxscore && (abs(refLength - j) < abs(p1 - p2)))){
                    p1 = refLength;
                    p2 = j;
                    maxscore = curScore;
                    segment_length = altLength - j;
                }
            }
        }
    }
    int lce_length = 0;
    struct CigarList* lce_head_node = NULL;
            
    if(segment_length > 0 && overhang_strategy == OVERHANG_STRATEGY_SOFTCLIP){
        if(!(lce_head_node = (struct CigarList*) malloc(sizeof(struct CigarList)))){
            throwAccError(env, "failed to malloc lce_head_node");
            //env->ThrowNew(Exception_SWPair,"failed to malloc lce_head_node");
            return -1;
        }
        lce_head_node->state = STATE_CLIP;
        lce_head_node->length = segment_length;
        lce_head_node->next = NULL;
        lce_length++;
        segment_length = 0;
    }

    int state = STATE_MATCH;
    int btr = 0;
    int new_state =0;
    int step_length = 0;
    struct CigarList* lce_new_node = NULL;
    do{
        btr = btrack[p1][p2];
        step_length = 1;
        if(btr > 0){
            new_state = STATE_DELETION;
            step_length = btr;
        }
        else if(btr < 0){
            new_state = STATE_INSERTION;
            step_length = (-btr);
        }
        else
            new_state = STATE_MATCH;

        switch(new_state){
            case STATE_MATCH: p1--; p2--; break;
            case STATE_INSERTION: p2 -= step_length; break; 
            case STATE_DELETION: p1 -= step_length; break;
        }

        if(new_state == state) segment_length += step_length;
        else{
            lce_new_node = NULL;
            if(!(lce_new_node = (struct CigarList*) malloc(sizeof(struct CigarList)))){
                throwAccError(env, "failed to malloc lce_new_node");
                //env->ThrowNew(Exception_SWPair, "failed to malloc lce_new_node");
                return -1;
            }
            lce_new_node->state = state;
            lce_new_node->length = segment_length;
            lce_new_node->next = lce_head_node;
            lce_head_node = lce_new_node;
            lce_length++;
            segment_length = step_length;
            state = new_state;
        }
    }while(p1 > 0 && p2 > 0);

    
    if(overhang_strategy == OVERHANG_STRATEGY_SOFTCLIP){
        if(!(lce_new_node = (struct CigarList*) malloc(sizeof(struct CigarList)))){
            throwAccError(env, "failed to malloc lce_new_node");
            //env->ThrowNew(Exception_SWPair, "failed to malloc lce_new_node");
            return -1;
        }   
        lce_new_node->state = state;
        lce_new_node->length = segment_length;
        lce_new_node->next = lce_head_node;
        lce_head_node = lce_new_node;
        lce_length++;

        if(p2 > 0){
            if(!(lce_new_node = (struct CigarList*) malloc(sizeof(struct CigarList)))){
                //fprintf(stderr, "No mem space for CigarList node\n");
                throwAccError(env, "failed to malloc lce_new_node");
                //env->ThrowNew(Exception_SWPair, "failed to malloc lce_new_node");
                return -1;
            } 
            lce_new_node->state = STATE_CLIP;
            lce_new_node->length = p2;
            lce_new_node->next = lce_head_node;
            lce_head_node = lce_new_node;
            lce_length++;
        }

        *alignment_offset_ptr = p1;
    }
    else if(overhang_strategy == OVERHANG_STRATEGY_IGNORE){
        if(!(lce_new_node = (struct CigarList*) malloc(sizeof(struct CigarList)))){
            throwAccError(env, "failed to malloc lce_new_node");
            //env->ThrowNew(Exception_SWPair, "failed to malloc lce_new_node");
            return -1;
        }   
        lce_new_node->state = state;
        lce_new_node->length = segment_length + p2;
        lce_new_node->next = lce_head_node;
        lce_head_node = lce_new_node;
        lce_length++;

        *alignment_offset_ptr = p1 - p2;
    }
    else{
        if(!(lce_new_node = (struct CigarList*) malloc(sizeof(struct CigarList)))){
            throwAccError(env, "failed to malloc lce_new_node");
            //env->ThrowNew(Exception_SWPair, "failed to malloc lce_new_node");
            return -1;
        }   
        lce_new_node->state = state;
        lce_new_node->length = segment_length;
        lce_new_node->next = lce_head_node;
        lce_head_node = lce_new_node;
        lce_length++;
        
        if(p1 > 0){
          if(!(lce_new_node = (struct CigarList*) malloc(sizeof(struct CigarList)))){
                throwAccError(env, "failed to malloc lce_new_node");
                //env->ThrowNew(Exception_SWPair, "failed to malloc lce_new_node");
                return -1;
            } 
            lce_new_node->state = STATE_DELETION;
            lce_new_node->length = p1;
            lce_new_node->next = lce_head_node;
            lce_head_node = lce_new_node;
            lce_length++;
        }
        else if(p2 > 0){
            if(!(lce_new_node = (struct CigarList*) malloc(sizeof(struct CigarList)))){
                throwAccError(env, "failed to malloc lce_new_node");
                //env->ThrowNew(Exception_SWPair, "failed to malloc lce_new_node");
                return -1;
            } 
            lce_new_node->state = STATE_INSERTION;
            lce_new_node->length = p2;
            lce_new_node->next = lce_head_node;
            lce_head_node = lce_new_node;
            lce_length++;
        }
        *alignment_offset_ptr = 0;
    }
    //lce_new_node is already a reversed linked list
    //now cast it into a int array
    if(!(*length_list_ptr = (int*)malloc(lce_length*sizeof(int)))){
        throwAccError(env, "failed to malloc length_list_ptr");
        //env->ThrowNew(Exception_SWPair, "failed to malloc length_list_ptr");
        return -1;
    }
    if(!(*state_list_ptr = (int*)malloc(lce_length*sizeof(int)))){
        throwAccError(env, "failed to malloc state_list_ptr");
        //env->ThrowNew(Exception_SWPair, "failed to malloc state_list_ptr");
        return -1;
    }
    struct CigarList* lce_list_ptr = lce_head_node;
    for(i = 0; i < lce_length; ++i){
        (*length_list_ptr)[i] = lce_list_ptr->length;
        (*state_list_ptr)[i] = lce_list_ptr->state;
        lce_list_ptr = lce_list_ptr->next;
    }
    *Cigar_list_length = lce_length;

    //free the lce list
    struct CigarList* lce_list_cur = lce_head_node;
    while(lce_list_cur){
       lce_list_ptr = lce_list_cur->next;
       free(lce_list_cur);
       lce_list_cur = lce_list_ptr;
    }
    return 0;
}*/

