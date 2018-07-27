/**
 * Copyright (C) 2015-2017 Falcon Computing Solutions, Inc.
 * All Rights Reserved.
 * All information contained herein is considered trade secrets and
 * confidential. Dissemination or reproduction, in whole or in part, is
 * strictly forbidden unless A prior written permission is obtained from
 * Falcon Computing Solutions.
 */

package com.falconcomputing.genomics.bqsr;

import htsjdk.samtools.SAMUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.broadinstitute.hellbender.utils.recalibration.*;
import org.apache.commons.math.optimization.fitting.GaussianFunction;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;

import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.AfterMethod;

public class RecalDatumTest extends RecalDatum {

  private final static Logger logger = LogManager.getLogger(RecalDatumTest.class);
  private final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

  // this part is initialized every time we run a test
  private FalconRecalibrationEngine engine;

  public RecalDatumTest() {
    super(1L, 0.1, (byte)64);
  }

  @Test(enabled = true)
  public void TestLog10Prior() {
    final double[] qemp = {20, 30, 40, 50, 60, 70, 80, 90};
    final double[] qreported = {20, 30, 40, 50, 60, 70, 80, 90};
    for (double e : qemp) {
      for (double r : qreported) {
        final double gatk_res = RecalDatum.log10QempPrior(e, r);
        //System.out.println(gatk_res);
        final double falcon_res = engine.log10QempPriorNative(e, r);
        //System.out.printf("java out: %f , native out: %f\n", gatk_res, falcon_res);
        Assert.assertEquals(falcon_res, gatk_res, 1e-8);
      }
    }
  }

  @Test(enabled = true)
  public void TestLog10Likelihood() {
    final double[] qemp = {20, 30, 40, 50, 60, 70, 80, 90};
    final long[] numObservations = {1, 2, 3, 10};
    final long[] numMismatches = {0, 1, 2, 5};
    for (double e : qemp) {
      for (long o : numObservations) {
        for (long m : numMismatches) {
          if (m >= o) continue;
          final double gatk_res = RecalDatum.log10QempLikelihood(e, o, m);
          final double falcon_res = engine.log10QempLikelihoodNative(e, o, m);
          if (gatk_res != falcon_res) {
            logger.info(String.format("error happens at %f, %d, %d with magnitude of %f", e, o, m, Math.abs(gatk_res - falcon_res)));
          }
          Assert.assertEquals(falcon_res, gatk_res, 1e-8);
        }
      }
    }
  }

 // @Test(enabled = true)
 // public void TestBayesianEstimate() {
 //   final long[] numObservations = {1, 2, 3, 10};
 //   final long[] numMismatches = {0, 1, 2, 5};
 //   final double[] qReported = {20, 30, 40, 50, 60, 70, 80, 90};

 //   for (double q : qReported) {
 //     for (long o : numObservations) {
 //       for (long m : numMismatches) {
 //         if (m >= o) continue;
 //         final double gatk_res = RecalDatum.bayesianEstimateOfEmpiricalQuality(o, m, q);
 //         final double falcon_res = engine.bayesianEstimateOfEmpiricalQualityNative(o, m, q);
 //         Assert.assertEquals(falcon_res, gatk_res, 1e-8);
 //       }
 //     }
 //   }
 // }

  @BeforeMethod
  public void setUp() {
    engine = new FalconRecalibrationEngine(RAC, null);
    final boolean isLoaded = engine.load(null);
    Assert.assertTrue(isLoaded);
  }

  @AfterMethod
  public void tearDown() {
    engine.finalizeData();
    engine = null;
  }
}
