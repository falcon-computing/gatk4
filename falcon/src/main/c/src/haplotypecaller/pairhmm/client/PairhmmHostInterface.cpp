#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <vector>

//#include "blaze/Timer.h"
//#define PLACE_TIMER ;
//#define PLACE_TIMER1(s) ;

#include "pairhmm/client/PairhmmHostInterface.h"

// Encode a scalar value to serialized data
template <typename T>
static inline void putT(void* buf, uint64_t & buf_idx, T value) {
  int len = sizeof(value);
  memcpy((char*)buf + buf_idx, reinterpret_cast<char*>(&value), len);
  buf_idx += len;
}

// Store a string with its length to serialized data
static inline void putStr(void* buf, uint64_t & buf_idx, const char* str, size_t len) {
  memcpy((char*)buf + buf_idx, str, len);
  buf_idx += len;
}

uint64_t serialize(void* buf, const read_t * reads, int num) {
  //PLACE_TIMER1("reads");

  uint64_t buf_idx = 0;
  putT(buf, buf_idx, num);

  for (int i = 0; i < num; i++) {
    int len = reads[i].len; 
    putT(buf, buf_idx, len);

    putStr(buf, buf_idx, reads[i]._b, len);
    putStr(buf, buf_idx, reads[i]._q, len);
    putStr(buf, buf_idx, reads[i]._i, len);
    putStr(buf, buf_idx, reads[i]._d, len);
    putStr(buf, buf_idx, reads[i]._c, len);
  }
  return buf_idx;
}

uint64_t serialize(void* buf, const hap_t * haps, int num) {
  //PLACE_TIMER1("haps");

  uint64_t buf_idx = 0;
  putT(buf, buf_idx, num);

  for (int i = 0; i < num; i++) {
    int len = haps[i].len;
    putT(buf, buf_idx, len);
    putStr(buf, buf_idx, haps[i]._b, len);
  }
  return buf_idx;
}


