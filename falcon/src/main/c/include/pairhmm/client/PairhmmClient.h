#ifndef PairHMMCLIENT_H
#define PairHMMCLIENT_H

#include <stdexcept>
#include <vector>
#include <unordered_map>

#include "Common.h"
#include "blaze/Client.h"
#include "pairhmm/client/m2m.h"
#include "pairhmm/client/PairhmmHostInterface.h"
#include "gkl-pairhmm/Context.h"
#include "gkl-pairhmm/avx_impl.h"

#define TOTAL_PE_NUM (SLR0_PE_NUM + SLR1_PE_NUM + SLR2_PE_NUM)


typedef struct {
    bool g_use_double;
    int  g_max_threads; /*never used*/
    bool g_use_fpga;    /*never used*/
    int  pmm_count;
    std::vector<read_t> reads;
    std::vector<hap_t>  haps;
    double             *likelihoods;
} PairhmmContext;

extern Context<float> g_ctxf;
extern Context<double> g_ctxd;
extern float m2m_table[128 * 128];
extern std::mutex pmm_ctx_map_mutex;
extern std::unordered_map<std::thread::id, PairhmmContext*> pmm_map;
extern uint8_t converstionTable[255];
extern std::pair<std::string, int> inputBank[3];

PairhmmContext* get_pairhmm_context();
void release_pairhmm_context();

class PairHMMClient : public blaze::Client {

 public:
  PairHMMClient();

  void setup(read_t* reads, int num_read, 
             hap_t*  haps,  int num_hap);

  void compute();

  int getMaxInputNumItems(int idx);
  int getMaxInputItemLength(int idx);
  int getMaxInputDataWidth(int idx);

 private:
  // used to perform easy calculation on cpu
  int     num_read_;
  int     num_hap_;
  read_t* reads_;
  hap_t*  haps_;
};

PairHMMClient* get_pairhmm_client();
void release_pairhmm_client();

#endif
