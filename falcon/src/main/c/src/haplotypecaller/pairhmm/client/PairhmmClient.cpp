#include <glog/logging.h>
#include <unordered_map>
#include <vector>

#include "blaze/Client.h"
#include "gkl-pairhmm/avx_impl.h"
#include "pairhmm/client/PairhmmWorker.h"
#include "pairhmm/client/PairhmmHostInterface.h"

Context<float> g_ctxf;
Context<double> g_ctxd;

uint8_t ConvertChar::conversionTable[255];


std::pair<std::string, int> inputBank[3];
float m2m_table[128 * 128] = M2M_INIT;
std::mutex pmm_ctx_map_mutex;
std::unordered_map<std::thread::id, PairhmmContext*> pmm_map;
std::mutex pmm_client_map_mutex;
std::unordered_map<std::thread::id, PairHMMClient*> pmm_client_map;

PairhmmContext* get_pairhmm_context(){
    std::thread::id this_id = std::this_thread::get_id();
    PairhmmContext* context;
    pmm_ctx_map_mutex.lock();
    std::unordered_map<std::thread::id, PairhmmContext*>::iterator iter = pmm_map.find(this_id);
    if(iter == pmm_map.end()){
        context = new PairhmmContext();
        context->pmm_count = 0;
        context->likelihoods = NULL;
        pmm_map[this_id] = context;
    }
    else{
        context = iter->second;
    }
    pmm_ctx_map_mutex.unlock();
    return context;
}

void release_pairhmm_context() {
    std::thread::id this_id = std::this_thread::get_id();
    pmm_ctx_map_mutex.lock();
    std::unordered_map<std::thread::id, PairhmmContext*>::iterator iter = pmm_map.find(this_id);
    if (iter != pmm_map.end()) {
        PairhmmContext* context = iter->second;
        delete context;
        pmm_map.erase(this_id);   
    }
    pmm_ctx_map_mutex.unlock();
}


using namespace blaze;

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

PairHMMClient::PairHMMClient(): 
  blaze::Client("PairHMM", 3, 1)
{
  createInput(1, 1, 4 + MAX_RSDATA_NUM*(4 + 5*MAX_READ_LEN), 1);
  createInput(2, 1, 4 + MAX_HAPDATA_NUM*(4 + MAX_HAP_LEN), 1);
}

inline testcase convert_avx_input(
    read_t & read,
    hap_t  & hap) 
{
  testcase ret;
  ret.rslen  = read.len; 
  ret.rs     = read._b; 
  ret.i      = read._i;
  ret.d      = read._d;
  ret.c      = read._c;
  ret.q      = read._q;
  ret.haplen = hap.len;
  ret.hap    = hap._b;
  return ret;
}

void PairHMMClient::setup(
    read_t* reads, int num_read, 
    hap_t*  haps,  int num_hap) 
{
  PLACE_TIMER;
  num_read_ = num_read;
  num_hap_  = num_hap;
  reads_ = reads;
  haps_  = haps;

  uint64_t total_rl = 0;
  uint64_t total_hl = 0;
  for (int i = 0; i < num_read; ++i) {
    total_rl += reads[i].len;
  }
  for (int i = 0; i < num_hap; ++i) {
    total_hl += haps[i].len;
  }
  num_cell_ = total_rl * total_hl;

  setInput(0, &num_cell_, 1, 1, sizeof(uint64_t));

  // serialize
  uint64_t read_size = serialize(getInputPtr(1), reads, num_read);
  uint64_t hap_size  = serialize(getInputPtr(2), haps, num_hap);

  // just use it to update block size
  createInput(1, 1, read_size, 1);
  createInput(2, 1, hap_size, 1);
}

// load balance compute routine
void PairHMMClient::compute() {
  ksight::AutoTimer __timer("compute on client cpu");
  ksight::ksight.add("num cells on cpu", num_cell_);

  float* output = (float*)createOutput(0, 
      1, num_read_ * num_hap_, sizeof(float));

  for (int i = 0; i < num_read_; ++i) {
    for (int j = 0; j < num_hap_; ++j) {
      // use GKL routines to perform computation
      testcase avx_v = convert_avx_input(reads_[i], haps_[j]);
      output[i*num_hap_ + j] = compute_fp_avxs(&avx_v);
    }
  }
}

