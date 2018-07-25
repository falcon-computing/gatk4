#include "pairhmm/client/PairhmmClient.h"

Context<float> g_ctxf;
Context<double> g_ctxd;

uint8_t ConvertChar::conversionTable[255];


std::pair<std::string, int> inputBank[3];
float m2m_table[128 * 128] = M2M_INIT;
std::mutex pmm_ctx_map_mutex;
std::map<std::thread::id, struct PairhmmContext*> pmm_map;

struct timespec pmm_diff_time(struct timespec start, struct timespec end)
{
    struct timespec temp;
    if ((end.tv_nsec-start.tv_nsec)<0) {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
}
struct PairhmmContext* init_pairhmm_context(){
    std::thread::id this_id = std::this_thread::get_id();
    std::map<std::thread::id, struct PairhmmContext*>::iterator iter = pmm_map.find(this_id);
    if(iter == pmm_map.end()){
        struct PairhmmContext* curCtx = new PairhmmContext();
        curCtx->client_input = new pairhmmInput();
        curCtx->host_input[0] = (FPGAInput*)_mm_malloc(sizeof(FPGAInput), 64);
        curCtx->host_input[1] = (FPGAInput*)_mm_malloc(sizeof(FPGAInput), 64);
        curCtx->host_input[2] = (FPGAInput*)_mm_malloc(sizeof(FPGAInput), 64);
        curCtx->overflow_avx_time = 0;
        curCtx->direct_avx_time = 0;
        curCtx->fpga_prepare_data_time = 0;
        curCtx->fpga_time = 0;
        curCtx->reject_time = 0;
        curCtx->total_time = 0;
        curCtx->pmm_count = 0;
        pmm_ctx_map_mutex.lock();
        pmm_map[this_id] = curCtx;
        pmm_ctx_map_mutex.unlock();
        return curCtx;
    }
    else{
        return iter->second;
    }

}

void PairhmmClient::compute(){
    struct PairhmmContext* curCtx = init_pairhmm_context();
    struct timespec time1, time2, diff;
    clock_gettime(CLOCK_REALTIME, &time1);
    curCtx->rejected = true;
    testcase testCase;
    for(int r = curCtx->read_base_index; r < curCtx->read_base_index + curCtx->cur_numRead; r++){
        testCase.rslen = curCtx->client_input->reads[r].bases.size();
        testCase.rs = curCtx->client_input->reads[r].bases.c_str(); 
        testCase.i = curCtx->client_input->reads[r]._i.c_str();
        testCase.d = curCtx->client_input->reads[r]._d.c_str();
        testCase.c = curCtx->client_input->reads[r]._c.c_str();
        testCase.q = curCtx->client_input->reads[r]._q.c_str();
        for(int j = curCtx->hap_base_index; j < curCtx->hap_base_index + curCtx->cur_numHap; j++){
            testCase.haplen = curCtx->client_input->haps[j].bases.size();
            testCase.hap = curCtx->client_input->haps[j].bases.c_str();
            curCtx->host_raw_output[r * curCtx->client_input->haps.size() + j] = compute_fp_avxs(&testCase);
        }
    }
    clock_gettime(CLOCK_REALTIME, &time2);
    diff = pmm_diff_time(time1, time2);
    curCtx->reject_time += diff.tv_sec * 1e9 + diff.tv_nsec;

}
    



