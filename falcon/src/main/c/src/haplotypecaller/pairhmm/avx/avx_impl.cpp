#include "pairhmm/avx/avx_impl.h"
#include "pairhmm/avx/avx-pairhmm.h"
#include "pairhmm/common/pairhmm_common.h"

float (*compute_fp_avxs)(testcase*) = &compute_full_prob_avxs;
double (*compute_fp_avxd)(testcase*) = &compute_full_prob_avxd;

