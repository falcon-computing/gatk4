#include "com_falconcomputing_genomics_NativeLibraryLoader.h"
#include "falcon-lic/genome.h"

JNIEXPORT jboolean JNICALL Java_com_falconcomputing_genomics_NativeLibraryLoader_license_1verify
  (JNIEnv *env, jclass cls) 
{
  return (license_verify() == 0);
}
