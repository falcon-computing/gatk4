/**
 * Copyright: Falcon Computing Solutions, Inc.
 * TODO: figure out the correct license text
 */
package com.falconcomputing.genomics.haplotypecaller;
import com.falconcomputing.genomics.NativeLibraryLoader;
import com.falconcomputing.genomics.AccelerationException;
/*
import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
*/
import com.intel.gkl.IntelGKLUtils;

import org.broadinstitute.gatk.nativebindings.NativeLibrary;
import org.apache.log4j.Logger;

import org.broadinstitute.gatk.nativebindings.pairhmm.HaplotypeDataHolder;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeBinding;
import org.broadinstitute.gatk.nativebindings.pairhmm.ReadDataHolder;

import java.io.File;

public class FalconPairhmm implements PairHMMNativeBinding{

    private final static String NATIVE_LIBRARY_NAME = "falcon_genomics";
    private final static Logger logger = Logger.getLogger(FalconPairhmm.class);
    private static boolean initialized = false;
    private static boolean loaded = false;
    private static boolean useFpga = false;
    private boolean profiling = false;

    public FalconPairhmm() {; }

    public synchronized boolean load(File tempDir) {
      if (!loaded) {
        if (!NativeLibraryLoader.load(tempDir, NATIVE_LIBRARY_NAME)) {
          return false;
        }
        loaded = true;
      }
      return true;
    }

    public boolean isInit() {
        return initialized;
    }
    
    public void initialize(PairHMMNativeArguments args) {
        if (!initialized) {
            if (args == null) {
                args = new PairHMMNativeArguments();
                args.useDoublePrecision = false;
                args.maxNumberOfThreads = 1;
            }
        
            initNative(ReadDataHolder.class, HaplotypeDataHolder.class, args.useDoublePrecision, args.maxNumberOfThreads, useFpga);
            // log information about threads
            if (args.maxNumberOfThreads > 1) {
                logger.warn("FalconPairhmm does not support OpenMP");   
            }
            initialized = true;
        }
        initialized = true;
    }

    public void computeLikelihoods(ReadDataHolder[] readDataArray, HaplotypeDataHolder[] haplotypeDataArray, double[] likelihoodArray){
        computeLikelihoodsNative(readDataArray, haplotypeDataArray, likelihoodArray);
    }

    public void done() {
        doneNative();
    }

    private native static void initNative(Class<?> readDataHolderClass, Class<?> haplotypeDataHolderClass, boolean doublePrecision, int maxThreads, boolean useFpga);

    private native void computeLikelihoodsNative(Object[] readDataArray, Object[] haplotypeDataArray, double[] likelihoodArray);

    private native void doneNative();
}
