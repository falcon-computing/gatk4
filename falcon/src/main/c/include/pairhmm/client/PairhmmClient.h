#ifndef PAIRHMMCLIENT_H
#define PAIRHMMCLIENT_H

#define LOG_HEADER "PairhmmClient"
#include <glog/logging.h>

#include "Common.h"
#include "blaze/Client.h"
#include "pairhmm/client/m2m.h"
#include "pairhmm/common/pairhmm_common.h"
#include "pairhmm/common/Context.h"
#include "pairhmm/avx/avx_impl.h"

#define TOTAL_PE_NUM (SLR0_PE_NUM + SLR1_PE_NUM + SLR2_PE_NUM)

struct PairhmmContext{
    bool g_use_double;
    int g_max_threads;
    bool g_use_fpga;
    pairhmmInput* client_input;
    FPGAInput* host_input[3];
    std::vector<float> host_raw_output;
    double overflow_avx_time;
    double direct_avx_time;
    double fpga_prepare_data_time;
    double fpga_time;
    double reject_time;
    double total_time;
    int pmm_count;
    bool rejected;
    int read_base_index;
    int hap_base_index;
    int cur_numRead;
    int cur_numHap;
    double* javaResults;
};



extern Context<float> g_ctxf;
extern Context<double> g_ctxd;
extern float m2m_table[128 * 128];
extern std::mutex pmm_ctx_map_mutex;
extern std::map<std::thread::id, struct PairhmmContext*> pmm_map;
extern uint8_t converstionTable[255];

extern std::pair<std::string, int> inputBank[3];


struct timespec pmm_diff_time(struct timespec start, struct timespec end);
struct PairhmmContext* init_pairhmm_context();
using namespace blaze;
class PairhmmClient: public Client {
    public:
    PairhmmClient(): Client("PairhmmTest", DIE_NUM * 3 + 1, DIE_NUM) {
        //ConvertChar::init();
        //_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    }
    void compute();
};


#endif
