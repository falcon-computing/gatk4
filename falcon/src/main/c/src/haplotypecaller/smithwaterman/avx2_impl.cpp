#include "avx2_impl.h"
#include "avx2-smithwaterman.h"
#include "Common.h"

std::mutex ptrs_map_mutex;
std::map<std::thread::id, struct sw_ptrs> sw_ptrs_map;
std::mutex cigar_ptrs_map_mutex;
std::map<std::thread::id, struct sw_ptrs> cigar_ptrs_map;

int32_t (*runSWOnePairBT_fp_avx2)(int32_t match, int32_t mismatch, int32_t open, int32_t extend,uint8_t *seq1, uint8_t *seq2, int32_t len1, int32_t len2, int8_t overhangStrategy, struct Cigar* cigarRet, struct sw_ptrs ptrs)= &runSWOnePairBT_avx2;

struct sw_ptrs init_ptrs(){
    std::thread::id this_id = std::this_thread::get_id();
    std::map<std::thread::id, struct sw_ptrs>::iterator iter = sw_ptrs_map.find(this_id);
    if(iter == sw_ptrs_map.end()){
#define AVX_LENGTH 8
        //fprintf(stderr,"thread id%d needs malloc matrices\n", this_id);
        struct sw_ptrs matrices;
        matrices.E_  = (int32_t *)_mm_malloc((6 * (MAX_SEQ_LEN+ AVX_LENGTH)) * sizeof(int32_t), 64);
        matrices.backTrack_ = (int16_t *)_mm_malloc((2 * MAX_SEQ_LEN * MAX_SEQ_LEN + 2 * AVX_LENGTH) * sizeof(int16_t), 64);
        matrices.cigarBuf_  = (int16_t *)_mm_malloc(4 * MAX_SEQ_LEN * sizeof(int16_t), 64);
        ptrs_map_mutex.lock();
        sw_ptrs_map[this_id] = matrices;
        ptrs_map_mutex.unlock();
        return matrices;
    }
    else{
        //fprintf(stderr,"thread id%d already have matrices\n", this_id);
        return iter->second;
    }
}

struct sw_ptrs init_cigar_ptrs(){
    std::thread::id this_id = std::this_thread::get_id();
    std::map<std::thread::id, struct sw_ptrs>::iterator iter = cigar_ptrs_map.find(this_id);
    if(iter == cigar_ptrs_map.end()){
#define AVX_LENGTH 8
        //fprintf(stderr,"thread id%d needs malloc matrices\n", this_id);
        struct sw_ptrs matrices;
        matrices.E_  = (int32_t *)_mm_malloc((6 * (MAX_SEQ_LEN+ AVX_LENGTH)) * sizeof(int32_t), 64);
        matrices.backTrack_ = (int16_t *)_mm_malloc((2 * MAX_SEQ_LEN * MAX_SEQ_LEN + 2 * AVX_LENGTH) * sizeof(int16_t), 64);
        matrices.cigarBuf_  = (int16_t *)_mm_malloc(4 * MAX_SEQ_LEN * sizeof(int16_t), 64);
        cigar_ptrs_map_mutex.lock();
        cigar_ptrs_map[this_id] = matrices;
        cigar_ptrs_map_mutex.unlock();
        return matrices;
    }
    else{
        return iter->second;
    }
}


int SWPairwiseAlignmentGKL(char* ref, int refLength, char alts[][MAX_SEQ_LENGTH], int batchSize, int* altLengths, struct Cigar* cigarResults, int* alignmentOffsets, int option){
    int i = 0;
    int overhang_strategy = 0;
    struct sw_ptrs ptrs = init_cigar_ptrs();
    for(i = 0; i < batchSize; ++i){
        uint8_t* ref_ptr = (uint8_t*)ref;
        uint8_t* alt_ptr = (uint8_t*)alts[i];
        alignmentOffsets[i] = runSWOnePairBT_fp_avx2(W_MATCH, W_MISMATCH, W_OPEN, W_EXTEND, ref_ptr, alt_ptr, refLength, altLengths[i], overhang_strategy, &cigarResults[i], ptrs);
            
    }
    return 0;
}


int SWPairwiseAlignmentOnceGKL(char* ref, int refLength, char* alt, int altLength, int w_match, int w_mismatch, int w_open, int w_extend, int overhang_strategy, int** length_list, int** state_list, int& Cigar_list_size, int& alignment_offset){
    uint8_t* ref_ptr = (uint8_t*)ref;
    uint8_t* alt_ptr = (uint8_t*)alt;
    struct Cigar cigar_ret;
    struct sw_ptrs ptrs = init_ptrs();

    alignment_offset = runSWOnePairBT_fp_avx2(w_match, w_mismatch, w_open, w_extend, ref_ptr, alt_ptr, refLength, altLength, overhang_strategy, &cigar_ret, ptrs);

    Cigar_list_size = cigar_ret.CigarElementNum;
    *length_list = (int*)malloc(Cigar_list_size * sizeof(int));
    *state_list = (int*)malloc(Cigar_list_size * sizeof(int));
    for(int i = 0; i < Cigar_list_size; i++){
        (*length_list)[i] = cigar_ret.cigarElements[i].length;
        (*state_list)[i] = cigar_ret.cigarElements[i].state;
    }

    return 0;

}


