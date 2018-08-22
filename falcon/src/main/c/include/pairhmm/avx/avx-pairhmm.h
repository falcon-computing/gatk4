#include <stdint.h>
#include "pairhmm/common/pairhmm_common.h"
#include "pairhmm/common/Context.h"

#include "pairhmm/avx/avx-types.h"

#undef SIMD_ENGINE
#define SIMD_ENGINE avx

#include "pairhmm/avx/avx-functions-float.h"
#include "pairhmm/avx/avx-vector-shift.h"
#include "pairhmm/avx/avx-pairhmm-template.h"

#include "pairhmm/avx/avx-functions-double.h"
#include "pairhmm/avx/avx-vector-shift.h"
#include "pairhmm/avx/avx-pairhmm-template.h"
