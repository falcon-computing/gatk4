/**
 * Copyright: Falcon Computing Solutions, Inc.
 * TODO: figure out the correct license text
 */
package com.falconcomputing.genomics.haplotypecaller;
import com.intel.gkl.pairhmm.IntelPairHmm;
import com.falconcomputing.genomics.NativeLibraryLoader;
import com.falconcomputing.genomics.AccelerationException;
import org.broadinstitute.hellbender.utils.pairhmm.*;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.gatk.nativebindings.NativeLibrary;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.broadinstitute.gatk.nativebindings.pairhmm.HaplotypeDataHolder;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeBinding;
import org.broadinstitute.gatk.nativebindings.pairhmm.ReadDataHolder;

import java.io.File;
import java.util.Random;
import org.testng.Assert;
import org.testng.annotations.Test;

public class FalconPairhmmTest{
    private final static Logger logger = LogManager.getLogger(FalconPairhmmTest.class);
    private PairHMMNativeArguments pairHmmNativeArgs = new PairHMMNativeArguments();
    
    private PairHMMNativeBinding golden_pairhmm = new IntelPairHmm();
    private PairHMMNativeBinding target_pairhmm = new FalconPairhmm(); 
    private ReadDataHolder[] reads;
    private ReadDataHolder[] reads_target;
    private HaplotypeDataHolder[] haps;
    private HaplotypeDataHolder[] haps_target;
    private double[] result_target;
    private double[] result_golden;

    
    private void randGenBytes(byte[] src, int length){
        Random randomGen = new Random();
        for(int i = 0; i < length; i++){
            src[i] = (byte)(randomGen.nextInt(50) + 10);
        }
    }
    private byte randBases(){
        Random randomGen = new Random();
        int i = randomGen.nextInt(10001);
        if(i < 2500)
            return 'A';
        else if(i < 5000)
            return 'C';
        else if(i < 7500)
            return 'T';
        else if(i < 10000)
            return 'G';
        else
            return 'N';
    }

    @Test(enabled = true)
    public void testFixedTest(){
        int numRead = 38;
        int numHap = 8;
        reads = new ReadDataHolder[numRead];
        reads_target = new ReadDataHolder[numRead];
        haps = new HaplotypeDataHolder[numHap];
        haps_target = new HaplotypeDataHolder[numHap];
        result_target = new double[numRead * numHap];
        result_golden = new double[numRead * numHap];
        Random randomGen = new Random();
        boolean intel_supported = golden_pairhmm.load(null);
        boolean falcon_supported = target_pairhmm.load(null);
        boolean isAccelerated = false;
        pairHmmNativeArgs.maxNumberOfThreads = 1;
        pairHmmNativeArgs.useDoublePrecision = false;
 
        Assert.assertTrue(intel_supported);
        if (!falcon_supported) return;

        golden_pairhmm.initialize(pairHmmNativeArgs);
        target_pairhmm.initialize(pairHmmNativeArgs);
        
        logger.info("reads size is " + reads.length + "haps size is " + haps.length);
        for(int i = 0; i < numRead; i++){
            int length = 100;
            reads[i] = new ReadDataHolder();
            reads[i].readBases = new byte[length];
            reads[i].readQuals = new byte[length];
            reads[i].insertionGOP = new byte[length];
            reads[i].deletionGOP = new byte[length];
            reads[i].overallGCP = new byte[length];
            for(int j = 0; j < length; j++){
                reads[i].readBases[j] = randBases();
                reads[i].readQuals[j] = (byte)(randomGen.nextInt(10) + 29);
                reads[i].insertionGOP[j] = (byte)(randomGen.nextInt(20) + 30);
                reads[i].deletionGOP[j] = (byte)(randomGen.nextInt(20) + 30);
                reads[i].overallGCP[j] = 10;
            }
            reads_target[i] = new ReadDataHolder();
            reads_target[i].readBases = new byte[length];
            reads_target[i].readQuals = new byte[length];
            reads_target[i].insertionGOP = new byte[length];
            reads_target[i].deletionGOP = new byte[length];
            reads_target[i].overallGCP = new byte[length];
            for(int j = 0; j < length; j++){
                reads_target[i].readBases[j] = reads[i].readBases[j];
                reads_target[i].readQuals[j] = reads[i].readQuals[j];
                reads_target[i].insertionGOP[j] = reads[i].insertionGOP[j];
                reads_target[i].deletionGOP[j] = reads[i].deletionGOP[j];
                reads_target[i].overallGCP[j] = reads[i].overallGCP[j];
            }

        }
        for(int i = 0; i < numHap; i++){
            int length = 100;
            haps[i] = new HaplotypeDataHolder();
            haps[i].haplotypeBases = new byte[length];
            for(int j = 0; j < length; j++){
                haps[i].haplotypeBases[j] = randBases();
            }
            haps_target[i] = new HaplotypeDataHolder();
            haps_target[i].haplotypeBases = new byte[length];
            for(int j = 0; j < length; j++){
                haps_target[i].haplotypeBases[j] = haps[i].haplotypeBases[j];
            }
        }
        double cells = 0;
        for(int i = 0; i < numRead; i++){
            for(int j = 0; j < numHap; j++){
                cells += reads[i].readBases.length * haps[j].haplotypeBases.length;
            }
        }   
        logger.info("done init inputs");
        long start_time = System.nanoTime();
        golden_pairhmm.computeLikelihoods(reads, haps, result_golden);
        long avx_time = System.nanoTime() - start_time;
        logger.info("done golden compute with " + avx_time * (1e-9) + "secs");
        logger.info("Intel AVX GCUPS is " + cells / avx_time); 

        start_time = System.nanoTime();
        logger.info("prepare for native");
        target_pairhmm.computeLikelihoods(reads_target, haps_target, result_target);
        long falcon_time = System.nanoTime() - start_time;
        logger.info("done target compute");
        logger.info("done target compute with " + falcon_time * (1e-9) + "secs");
        logger.info("Falcon GCUPS is " + cells / falcon_time); 
        target_pairhmm.done();
        logger.info("tore down target");
        golden_pairhmm.done();
        logger.info("tore down golden");
        double count_err = 0;
        for(int i = 0; i < numRead * numHap; i++){
            double diff = Math.abs((result_target[i] - result_golden[i]) / result_golden[i]);
            //logger.info("reads id: " + i / numHap + "haps id: " + i % numHap + "golden = " + result_golden[i] + " target = " + result_target[i]);
            //logger.info("reads length " + reads[i / numHap].readBases.length + " haps length " + haps[i % numHap].haplotypeBases.length);
            if(result_golden[i] == Double.POSITIVE_INFINITY){
                if(result_target[i] != Double.POSITIVE_INFINITY)
                    count_err = count_err + 1;
                Assert.assertTrue(result_target[i] == Double.POSITIVE_INFINITY);
            }
            else if(result_golden[i] == Double.NEGATIVE_INFINITY){
                if(result_target[i] != Double.NEGATIVE_INFINITY)
                    count_err = count_err + 1;
                Assert.assertTrue(result_target[i] == Double.NEGATIVE_INFINITY);
            }
            else{
                if(diff >= 2e-3)
                    count_err = count_err + 1;
            }
            Assert.assertTrue(diff < 2e-3);
        }
        logger.info(count_err + " out of " + numRead * numHap + " runs is wrong");
    }

    @Test(enabled = true)
      public void testHugeTest() {

        int numRead = 2050;
        int numHap = 32;
        reads = new ReadDataHolder[numRead];
        reads_target = new ReadDataHolder[numRead];
        haps = new HaplotypeDataHolder[numHap];
        haps_target = new HaplotypeDataHolder[numHap];
        result_target = new double[numRead * numHap];
        result_golden = new double[numRead * numHap];
        Random randomGen = new Random();
        boolean intel_supported = golden_pairhmm.load(null);
        boolean falcon_supported = target_pairhmm.load(null);
        boolean isAccelerated = false;
        pairHmmNativeArgs.maxNumberOfThreads = 1;
        pairHmmNativeArgs.useDoublePrecision = false;

        Assert.assertTrue(intel_supported);
        if (!falcon_supported) return;

        golden_pairhmm.initialize(pairHmmNativeArgs);
        target_pairhmm.initialize(pairHmmNativeArgs);

        logger.info("reads size is " + reads.length + "haps size is " + haps.length);
        for(int i = 0; i < numRead; i++){
          int length = randomGen.nextInt(191) + 1;
          reads[i] = new ReadDataHolder();
          reads[i].readBases = new byte[length];
          reads[i].readQuals = new byte[length];
          reads[i].insertionGOP = new byte[length];
          reads[i].deletionGOP = new byte[length];
          reads[i].overallGCP = new byte[length];
          for(int j = 0; j < length; j++){
            reads[i].readBases[j] = randBases();
            reads[i].readQuals[j] = (byte)(randomGen.nextInt(10) + 29);
            reads[i].insertionGOP[j] = (byte)(randomGen.nextInt(20) + 30);
            reads[i].deletionGOP[j] = (byte)(randomGen.nextInt(20) + 30);
            reads[i].overallGCP[j] = 10;
          }
          reads_target[i] = new ReadDataHolder();
          reads_target[i].readBases = new byte[length];
          reads_target[i].readQuals = new byte[length];
          reads_target[i].insertionGOP = new byte[length];
          reads_target[i].deletionGOP = new byte[length];
          reads_target[i].overallGCP = new byte[length];
          for(int j = 0; j < length; j++){
            reads_target[i].readBases[j] = reads[i].readBases[j];
            reads_target[i].readQuals[j] = reads[i].readQuals[j];
            reads_target[i].insertionGOP[j] = reads[i].insertionGOP[j];
            reads_target[i].deletionGOP[j] = reads[i].deletionGOP[j];
            reads_target[i].overallGCP[j] = reads[i].overallGCP[j];
          }

        }
        for(int i = 0; i < numHap; i++){
          int length = randomGen.nextInt(1023) + 1;
          haps[i] = new HaplotypeDataHolder();
          haps[i].haplotypeBases = new byte[length];
          for(int j = 0; j < length; j++){
            haps[i].haplotypeBases[j] = randBases();
          }
          haps_target[i] = new HaplotypeDataHolder();
          haps_target[i].haplotypeBases = new byte[length];
          for(int j = 0; j < length; j++){
            haps_target[i].haplotypeBases[j] = haps[i].haplotypeBases[j];
          }
        }
        double cells = 0;
        for(int i = 0; i < numRead; i++){
          for(int j = 0; j < numHap; j++){
            cells += reads[i].readBases.length * haps[j].haplotypeBases.length;
          }
        }   
        logger.info("done init inputs");
        long start_time = System.nanoTime();
        golden_pairhmm.computeLikelihoods(reads, haps, result_golden);
        long avx_time = System.nanoTime() - start_time;
        logger.info("done golden compute with " + avx_time * (1e-9) + "secs");
        logger.info("Intel AVX GCUPS is " + cells / avx_time); 

        start_time = System.nanoTime();
        target_pairhmm.computeLikelihoods(reads_target, haps_target, result_target);
        long falcon_time = System.nanoTime() - start_time;
        logger.info("done target compute");
        logger.info("done target compute with " + falcon_time * (1e-9) + "secs");
        logger.info("Falcon GCUPS is " + cells / falcon_time); 
        target_pairhmm.done();
        golden_pairhmm.done();
        double count_err = 0;
        for(int i = 0; i < numRead * numHap; i++){
          double diff = Math.abs((result_target[i] - result_golden[i]) / result_golden[i]);
          //logger.info("reads id: " + i / numHap + "haps id: " + i % numHap + "golden = " + result_golden[i] + " target = " + result_target[i]);
          //logger.info("reads length " + reads[i / numHap].readBases.length + " haps length " + haps[i % numHap].haplotypeBases.length);
          if(result_golden[i] == Double.POSITIVE_INFINITY){
            if(result_target[i] != Double.POSITIVE_INFINITY)
              count_err = count_err + 1;
            Assert.assertTrue(result_target[i] == Double.POSITIVE_INFINITY);
          }
          else if(result_golden[i] == Double.NEGATIVE_INFINITY){
            if(result_target[i] != Double.NEGATIVE_INFINITY)
              count_err = count_err + 1;
            Assert.assertTrue(result_target[i] == Double.NEGATIVE_INFINITY);
          }
          else{
            if(diff >= 2e-3)
              count_err = count_err + 1;
          }
          Assert.assertTrue(diff < 2e-3);
        }
        logger.info(count_err + " out of " + numRead * numHap + " runs is wrong");
      }


}


