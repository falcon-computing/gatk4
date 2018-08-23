#ifndef FALCONGENOMICS_COMMON_H
#define FALCONGENOMICS_COMMON_H

#include <jni.h>
#include <sys/syscall.h>
#include <sys/time.h>
#include <time.h>
#include <mutex>
#include <thread>
#include <map>

#include <stdexcept>
#include <string>
#define LOG_HEADER "Falcon Genomics Library"
#include <glog/logging.h>

#define TIMER_SAMPLE_RATE 1000000

inline uint64_t getUs() {
  struct timespec tr;
  clock_gettime(CLOCK_REALTIME, &tr);

  return (uint64_t)tr.tv_sec*1e6 + tr.tv_nsec/1e3;
}

inline uint64_t getNs() {
  struct timespec tr;
  clock_gettime(CLOCK_REALTIME, &tr);

  return (uint64_t)tr.tv_sec*1e9 + tr.tv_nsec;
}

inline jint throwAccError(JNIEnv *env, const char* message) {
  jclass cls_exception = env->FindClass("com/falconcomputing/genomics/AccelerationException");
  if (cls_exception) {
    return env->ThrowNew(cls_exception, message);
  }
  else {
    return -1;
  }
}

#endif

