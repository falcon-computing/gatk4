/**
 * Copyright Falcon Computing Solutions, Inc.
 * TODO: figure out the correct license text
 */
package com.falconcomputing.genomics.bqsr;

import com.falconcomputing.genomics.NativeLibraryLoader;
import com.falconcomputing.genomics.AccelerationException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.tribble.Feature;

//import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.hellbender.engine.filters.*;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.recalibration.*;
import org.broadinstitute.hellbender.utils.recalibration.covariates.*;

//import org.broadinstitute.gatk.engine.datasources.reference.*;
import org.broadinstitute.hellbender.engine.datasources.*;

//import org.broadinstitute.gatk.nativebindings.NativeLibrary;
import org.broadinstitute.gatk.nativebindings.*;

//import org.broadinstitute.gatk.tools.walkers.bqsr.*;
import org.broadinstitute.hellbender.tools.walkers.bqsr.*;

//import org.broadinstitute.gatk.utils.commandline.Advanced;

//import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.barclay.argparser.Argument;


//import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.barclay.argparser.ArgumentCollection;


//import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
//import org.broadinstitute.gatk.utils.contexts.ReferenceContext.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceContext.*;

//import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
//TODO

//import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;

//import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.BaseUtils;

//import org.broadinstitute.gatk.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;

//import org.broadinstitute.gatk.utils.collections.Pair;
//TODO
import org.apache.commons.lang3.tuple.Pair;

//import org.broadinstitute.gatk.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;

//import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
//TODO

//import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

//import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;

//import org.broadinstitute.gatk.utils.recalibration.*;
// already imported

//import org.broadinstitute.gatk.utils.sam.*;



import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.AfterMethod;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.*;
import java.util.*;
import java.net.*;

public class FalconRecalibrationEngineTest {

  private final static Logger logger = LogManager.getLogger(FalconRecalibrationEngineTest.class);
  private final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();

  //TODO
  //private TestHelper helper;

  // this part is initialized every time we run a test
  private FalconRecalibrationEngine engine;

  // file input for testing
  private String refFilename = "ref.fa";
  private String bamFilename = "input-1.bam";
  private String grpFilename = "input-1.grp";

  private Path refPath; // path to reference genome
  private Path bamPath; // path to input bam
  private Path grpPath; // path to bqsr report

  public FalconRecalibrationEngineTest() {

    // read system properties for input filenames
    if (System.getProperty("bqsr.input_bam") != null) {
      this.bamFilename = System.getProperty("bqsr.input_bam");
      logger.info("Use bam file: " + this.bamFilename);
    }
    if (System.getProperty("bqsr.reference") != null) {
      this.refFilename = System.getProperty("bqsr.reference");
      logger.info("Use reference: " + this.refFilename);
    }
    if (System.getProperty("bqsr.input_grp") != null) {
      this.grpFilename = System.getProperty("bqsr.input_grp");
      logger.info("Use table: " + this.grpFilename);
    }

    // get path for input data
    try {
      URL bam_url = this.getClass().getResource(bamFilename);
      URL ref_url = this.getClass().getResource(refFilename);
      URL grp_url = this.getClass().getResource(grpFilename);

      bamPath = Paths.get(bam_url.toURI());
      refPath = Paths.get(ref_url.toURI());
      grpPath = Paths.get(grp_url.toURI());
    }
    catch (URISyntaxException e) {
      e.printStackTrace();
      return;
    }
    //TODO
    //helper = new TestHelper(refPath.toFile());
  }

  @Test(enabled = true, groups = {"bqsr"})
  public void TestLicense() {
    ;
  }

  @Test(enabled = true, groups = {"bqsr"})
  public void TestBasicInterface() {

    //final Covariate[] covariates = getCovariates();
    final StandardCovariateList covariates = getCovariates();
    final int numReadGroups = 2;
    //int numCovariates = covariates.length;
    int numCovariates = covariates.size();
    int numEvents = EventType.values().length;
    //int qualLength = covariates[1].maximumKeyValue()+1;
    int qualLength = covariates.get(1).maximumKeyValue()+1;

    // Initialize the table
    try {
      engine.init(covariates, numReadGroups);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }

    final FalconRecalibrationEngine.RecalDatumTable[] recalTables = engine.getTables();
    Assert.assertEquals(recalTables.length, numCovariates-1);

    // Test the dimensions of the table to make sure the interface is correct
    for (int i = 1; i < numCovariates; i++) { // skip the first table
      //int covLength = covariates[i].maximumKeyValue()+1;
      int covLength = covariates.get(1).maximumKeyValue()+1;
      int expectedTableSize = numReadGroups*numEvents;
      expectedTableSize = expectedTableSize*qualLength;
      if (i > 1)
        expectedTableSize = expectedTableSize*covLength;

      // check the dimensions of the returned table
      Assert.assertEquals(recalTables[i-1].numOccurance.length, expectedTableSize);
      Assert.assertEquals(recalTables[i-1].numMismatches.length, expectedTableSize);
      Assert.assertEquals(recalTables[i-1].tableSize, expectedTableSize);
      if (i == 1) {
        Assert.assertEquals(recalTables[i-1].tableDimensionsSize, 3);
      }
      else {
        Assert.assertEquals(recalTables[i-1].tableDimensionsSize, 4);
      }
      Assert.assertEquals(recalTables[i-1].tableDimensions[0], numEvents);
      Assert.assertEquals(recalTables[i-1].tableDimensions[1], numReadGroups);
      Assert.assertEquals(recalTables[i-1].tableDimensions[2], qualLength);
      if (i > 1) {
        Assert.assertEquals(recalTables[i-1].tableDimensions[3], covLength);
      }
    }
  }
  /*
  @Test(enabled = true, groups = {"bqsr"})
  public void TestCycleCovariatesOnSynthesizedData() {

    final int readLength = 10;
    final byte[] bases = new byte[readLength];
    final byte[] baseQuals = new byte[readLength];
    final byte[] insertionQuals = new byte[readLength];
    final byte[] deletionQuals = new byte[readLength];
    final boolean[] skips = new boolean[readLength];
    final double[] snpErrors = new double[readLength];
    final double[] insertionErrors = new double[readLength];
    final double[] deletionsErrors = new double[readLength];

    for ( int i = 0; i < readLength; i++ ) {
      bases[i] = 'A';
      baseQuals[i] = (byte)(i % SAMUtils.MAX_PHRED_SCORE);
      insertionQuals[i] = (byte)((i+1) % SAMUtils.MAX_PHRED_SCORE);
      deletionQuals[i] = (byte)((i+2) % SAMUtils.MAX_PHRED_SCORE);
      skips[i] = i % 2 == 0;
      snpErrors[i] = 1.0 / (i+1);
      insertionErrors[i] = 0.5 / (i+1);
      deletionsErrors[i] = 0.3 / (i+1);
    }

    final EnumMap<EventType, double[]> errors = new EnumMap<EventType, double[]>(EventType.class);
    errors.put(EventType.BASE_SUBSTITUTION, snpErrors);
    errors.put(EventType.BASE_INSERTION, insertionErrors);
    errors.put(EventType.BASE_DELETION, deletionsErrors);

    final EnumMap<EventType, byte[]> quals = new EnumMap<EventType, byte[]>(EventType.class);
    quals.put(EventType.BASE_SUBSTITUTION, baseQuals);
    quals.put(EventType.BASE_INSERTION, insertionQuals);
    quals.put(EventType.BASE_DELETION, deletionQuals);

    final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(bases, baseQuals, readLength + "M");
    read.setBaseQualities(insertionQuals, EventType.BASE_INSERTION);
    read.setBaseQualities(deletionQuals, EventType.BASE_DELETION);

    // set read group to test unsupported platform
    final GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord("@RG\tID:test\tSM:test\tPL:bgiseq");
    read.setReadGroup(rg);

    try {
      final int[][][] falcon_keys = engine.computeCycleCovariates(read);
    }
    catch (AccelerationException e) {
      logger.info("caught exception for unsupported platform");
      return;
    }
    Assert.fail("should have caught exception");
  }

  @Test(enabled = true, groups = {"bqsr"})
  public void TestCycleCovariates() {
    final Covariate[] covariates = getCovariates();
    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();
    final int numReadGroups = header.getReadGroups().size();
    try {
      engine.init(covariates, numReadGroups);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }

    for (SAMRecord record : reader) {
      final GATKSAMRecord read = new GATKSAMRecord(record);
      final ReadCovariates cov = RecalUtils.computeCovariates(read, covariates);
      try {
        final int[][][] falcon_keys = engine.computeCycleCovariates(read);
        final int contextCovIdx = 3;

        for (EventType event : EventType.values()) {
          //System.out.println(Arrays.toString(falcon_keys[event.ordinal()]));
          final int[][] gatk_keys = cov.getKeySet(event);
          for (int i = 0; i < read.getReadBases().length; i++) {
            //System.out.println(Arrays.toString(gatk_keys[i]));
            Assert.assertEquals(
                  falcon_keys[event.ordinal()][i][contextCovIdx],
                  gatk_keys[i][contextCovIdx]);
          }
        }
      }
      catch (AccelerationException e) {
        Assert.fail("should not caught exception here: " + e.getMessage());
      }
    }
  }

  @Test(enabled = true, groups = {"bqsr"})
  public void TestContextCovariates() {
    final Covariate[] covariates = getCovariates();
    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();
    final int numReadGroups = header.getReadGroups().size();
    logger.info(String.format("number of read groups: %d", numReadGroups));

    try {
      engine.init(covariates, numReadGroups);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }

    for (SAMRecord record : reader) {
      final GATKSAMRecord read = new GATKSAMRecord(record);
      final ReadCovariates cov = RecalUtils.computeCovariates(read, covariates);

      final int[][][] falcon_keys = engine.computeContextCovariates(read);

      final int contextCovIdx = 2;
      for (EventType event : EventType.values()) {
        //System.out.println(Arrays.toString(falcon_keys[event.ordinal()]));
        final int[][] gatk_keys = cov.getKeySet(event);
        for (int i = 0; i < read.getReadBases().length; i++) {
          //System.out.println(Arrays.toString(gatk_keys[i]));
          Assert.assertEquals(
              falcon_keys[event.ordinal()][i][contextCovIdx],
              gatk_keys[i][contextCovIdx]);
        }
      }
    }
  }

  @Test(enabled = true, groups = {"bqsr"})
  public void TestCovariates() {
    final Covariate[] covariates = getCovariates();
    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();
    final int numReadGroups = header.getReadGroups().size();

    try {
      engine.init(covariates, numReadGroups);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }

    int numCovariates = covariates.length;

    for (SAMRecord record : reader) {
      final GATKSAMRecord read = new GATKSAMRecord(record);
      final ReadCovariates cov = RecalUtils.computeCovariates(read, covariates);
      final int[] falcon_keys = engine.computeCovariates(read);
      int readLength = read.getReadBases().length;

      int idx = 0;
      for (int i = 0; i < readLength; i++) {
        for (int j = 0; j < numCovariates; j++) {
          for (EventType event : EventType.values()) {
            final int gatk_key = cov.getKeySet(event)[i][j];
            Assert.assertEquals(falcon_keys[idx++], gatk_key);
          }
        }
      }
    }
  }

  @Test(enabled = true, groups = {"bqsr"})
  public void TestReadGroupCovariate() {
    final Covariate[] covariates = getCovariates();
    // create a separate covariate class to compare
    final Covariate[] gatk_covariates = getCovariates();

    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();
    final int numReadGroups = header.getReadGroups().size();
    final int numCovariates = covariates.length;

    try {
      engine.init(covariates, numReadGroups);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }

    int numRecords = 0;
    for (SAMRecord record : reader) {
      final GATKSAMRecord read = new GATKSAMRecord(record);
      final ReadCovariates cov = RecalUtils.computeCovariates(read, gatk_covariates);
      final int[] falcon_keys = engine.computeCovariates(read);
    }
    engine.updateReadGroupCovariates();

    // check read group covariates
    final ReadGroupCovariate falcon_rgCov = (ReadGroupCovariate)covariates[0];
    final ReadGroupCovariate gatk_rgCov = (ReadGroupCovariate)gatk_covariates[0];
    Assert.assertEquals(gatk_rgCov.maximumKeyValue(), falcon_rgCov.maximumKeyValue());

    for (int i = 0; i < gatk_rgCov.maximumKeyValue(); i++) {
      final String gatk_rg = gatk_rgCov.formatKey(i);
      final String falcon_rg = falcon_rgCov.formatKey(i);
      Assert.assertEquals(gatk_rg, falcon_rg);
    }
  }

  @Test(enabled = true, groups = {"bqsr"})
  public void TestBAQCalculation() {
    final Covariate[] covariates = getCovariates();
    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();
    final int numReadGroups = header.getReadGroups().size();
    final int numCovariates = covariates.length;

    try {
      engine.init(covariates, numReadGroups);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }

    RecalibrationEngine recalibrationEngine = new RecalibrationEngine(covariates, numReadGroups, RAC.RECAL_TABLE_UPDATE_LOG, false);

    int numRecords = 0;
    for (SAMRecord record : reader) {
      //if (numRecords > 1) break;
      final GATKSAMRecord org_read = new GATKSAMRecord(record);
      final GATKSAMRecord read = ReadClipper.hardClipSoftClippedBases(ReadClipper.hardClipAdaptorSequence(org_read));
      final ReferenceContext ref = helper.getRefContext(org_read);

      final int[] isSNP = helper.falconCalculateIsSNP(read, ref, org_read);
      final int[] isInsertion = helper.falconCalculateIsIndel(read, EventType.BASE_INSERTION);
      final int[] isDeletion = helper.falconCalculateIsIndel(read, EventType.BASE_DELETION);
      final int nErrors = helper.falconNumEvents(isSNP, isInsertion, isDeletion);

      if (nErrors == 0) continue;

      final byte[] gatk_baqArray = helper.falconCalculateBAQArray(read);
      final byte[] falcon_baqArray = engine.calculateBAQArray(read);

      if (gatk_baqArray == null) {
        Assert.assertEquals(falcon_baqArray, null);
      }
      else {
        Assert.assertEquals(falcon_baqArray.length, gatk_baqArray.length);
        //System.out.println(Arrays.toString(falcon_baqArray));
        //System.out.println(Arrays.toString(gatk_baqArray));
        for (int i = 0; i < gatk_baqArray.length; i++) {
          Assert.assertEquals(falcon_baqArray[i], gatk_baqArray[i]);
        }
      }
      //System.out.println(String.format("finish read %d", numRecords));
      numRecords++;
    }
  }

  @Test(enabled = true, groups = {"bqsr"})
  public void TestFractionalErrors() {
    final Covariate[] covariates = getCovariates();
    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();
    final int numReadGroups = header.getReadGroups().size();
    final int numCovariates = covariates.length;

    try {
      engine.init(covariates, numReadGroups);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }

    RecalibrationEngine recalibrationEngine = new RecalibrationEngine(covariates, numReadGroups, RAC.RECAL_TABLE_UPDATE_LOG, false);

    int numRecords = 0;
    for (SAMRecord record : reader) {
      //if (numRecords > 1) break;
      final GATKSAMRecord org_read = new GATKSAMRecord(record);
      final GATKSAMRecord read = ReadClipper.hardClipSoftClippedBases(ReadClipper.hardClipAdaptorSequence(org_read));
      final ReferenceContext ref = helper.getRefContext(org_read);

      int readLength = read.getReadBases().length;

      final int[] isSNP = helper.falconCalculateIsSNP(read, ref, org_read);
      final int[] isInsertion = helper.falconCalculateIsIndel(read, EventType.BASE_INSERTION);
      final int[] isDeletion = helper.falconCalculateIsIndel(read, EventType.BASE_DELETION);
      final int nErrors = helper.falconNumEvents(isSNP, isInsertion, isDeletion);

      final byte[] baqArray = nErrors == 0 ? helper.falconFlatBAQArray(read) : helper.falconCalculateBAQArray(read);

      //logger.info("isSNP = " + Arrays.toString(isSNP));
      //logger.info("isInsertion = " + Arrays.toString(isInsertion));
      //logger.info("isDeletion = " + Arrays.toString(isDeletion));
      //logger.info("nErrors = " + Integer.toString(nErrors));

      final double[][] errors = engine.calculateFractionalErrorArray(read, org_read, ref);

      if (baqArray == null) { // some reads just can't be BAQ'ed
        Assert.assertEquals(null, errors);
      }
      else {
        //logger.info("baqArray = " + Arrays.toString(baqArray));
        //logger.info(String.format("read %d has %d errors, and baqArray is not null", numRecords, nErrors));

        //final boolean[] skip = calculateSkipArray(read, metaDataTracker); // skip known sites of variation as well as low quality and non-regular bases
        final double[] snpErrors = helper.falconCalculateFractionalErrorArray(isSNP, baqArray);
        final double[] insertionErrors = helper.falconCalculateFractionalErrorArray(isInsertion, baqArray);
        final double[] deletionErrors = helper.falconCalculateFractionalErrorArray(isDeletion, baqArray);

        //logger.info("snpErrors = " + Arrays.toString(snpErrors));
        //logger.info("insertionErrors = " + Arrays.toString(insertionErrors));
        //logger.info("deletionErrors = " + Arrays.toString(deletionErrors));
        Assert.assertNotNull(errors);

        for (int i = 0; i < readLength; i++) {
          Assert.assertEquals(snpErrors[i], errors[0][i]);
          Assert.assertEquals(insertionErrors[i], errors[1][i]);
          Assert.assertEquals(deletionErrors[i], errors[2][i]);
        }
      }

      //System.out.println(String.format("finish read %d", numRecords));
      numRecords++;
    }
  }

  @Test(enabled = true, groups = {"bqsr"})
  public void TestTableUpdateWithRealData() {
    final Covariate[] covariates = getCovariates();
    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();
    final int numReadGroups = header.getReadGroups().size();
    final int numCovariates = covariates.length;

    try {
      engine.init(covariates, numReadGroups);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }

    RecalibrationEngine recalibrationEngine = new RecalibrationEngine(covariates, numReadGroups, RAC.RECAL_TABLE_UPDATE_LOG, false);

    int numRecords = 0;
    for (SAMRecord record : reader) {
      final GATKSAMRecord org_read = new GATKSAMRecord(record);
      final GATKSAMRecord read = ReadClipper.hardClipSoftClippedBases(ReadClipper.hardClipAdaptorSequence(org_read));
      final ReferenceContext ref = helper.getRefContext(org_read);

      final int readLength = read.getReadBases().length;
      final boolean[] skip = new boolean[readLength];
      Arrays.fill(skip, false);

      // perform falcon table update
      int ret = 0;
      try {
        ret = engine.update(read, org_read, ref, skip);
      }
      catch (AccelerationException e) {
        logger.error("exception caught in init(): "+ e.getMessage());
        return;
      }

      final int[] isSNP = helper.falconCalculateIsSNP(read, ref, org_read);
      final int[] isInsertion = helper.falconCalculateIsIndel(read, EventType.BASE_INSERTION);
      final int[] isDeletion = helper.falconCalculateIsIndel(read, EventType.BASE_DELETION);
      final int nErrors = helper.falconNumEvents(isSNP, isInsertion, isDeletion);
      final byte[] baqArray = nErrors == 0 ? helper.falconFlatBAQArray(read) : helper.falconCalculateBAQArray(read);

      if (baqArray != null) { // some reads just can't be BAQ'ed
        //final boolean[] skip = calculateSkipArray(read, metaDataTracker); // skip known sites of variation as well as low quality and non-regular bases

        final double[] snpErrors = helper.falconCalculateFractionalErrorArray(isSNP, baqArray);
        final double[] insertionErrors = helper.falconCalculateFractionalErrorArray(isInsertion, baqArray);
        final double[] deletionErrors = helper.falconCalculateFractionalErrorArray(isDeletion, baqArray);

        final ReadCovariates cov = RecalUtils.computeCovariates(read, covariates);
        // aggregate all of the info into our info object, and update the data
        final ReadRecalibrationInfo info = new ReadRecalibrationInfo(read, cov, skip, snpErrors, insertionErrors, deletionErrors);
        recalibrationEngine.updateDataForRead(info);
      }
      else {
        Assert.assertEquals(ret, 0);
      }
      numRecords++;
    }

    // get table results
    recalibrationEngine.finalizeData();

    RecalibrationTables gatk_table = recalibrationEngine.getFinalRecalibrationTables();
    RecalibrationTables our_table = engine.getFinalRecalibrationTables();

    // compare all tables
    for (int i = 0; i < numCovariates; i++) {
      List<RecalDatum> gatk_table_contents = gatk_table.getTable(i).getAllValues();
      List<RecalDatum> our_table_contents = our_table.getTable(i).getAllValues();
      Assert.assertEquals(our_table_contents.size(), gatk_table_contents.size());
      //System.out.println(String.format("%d: %d == %d", i, our_table_contents.size(), gatk_table_contents.size()));
      for (int k = 0; k < gatk_table_contents.size(); k++) {
        //System.out.println(String.format("[%d-O] %d == %d", k, our_table_contents.get(k).getNumObservations(), gatk_table_contents.get(k).getNumObservations()));
        //System.out.println(String.format("[%d-M] %f == %f", k, our_table_contents.get(k).getNumMismatches(), gatk_table_contents.get(k).getNumMismatches()));
        compareRecalDatum(our_table_contents.get(k), gatk_table_contents.get(k));
      }
    }
  }

  @Test(enabled = true, groups = {"pr"})
  public void TestInitForRecalibrate() {
    final RecalibrationReport report = getRecalReport();

    final RecalibrationTables gatk_tables = report.getRecalibrationTables();
    final Covariate[] requestedCovariates = report.getRequestedCovariates();
    final QuantizationInfo quantizationInfo = report.getQuantizationInfo();

    // here assuming staticQuantizedMapping is null
    final List<Byte> quantizedQuals = quantizationInfo.getQuantizedQuals();

    long start_ts = System.nanoTime();

    // use default parameters
    try {
      engine.init(requestedCovariates, gatk_tables,
            quantizedQuals, null,
            false, QualityUtils.MIN_USABLE_Q_SCORE,
            -1.0, false);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }
    logger.info(String.format("init table for PR takes %f ms", (System.nanoTime() - start_ts)/1e6));

    final RecalibrationTables falcon_tables = engine.getFinalRecalibrationTables();

    compareRecalibrationTables(requestedCovariates.length, falcon_tables, gatk_tables);
  }

  @Test(enabled = true, groups = {"pr"})
  public void TestRecalibrate() {
    final boolean disableIndelQuals = false;
    final int preserveQLessThan = QualityUtils.MIN_USABLE_Q_SCORE;
    final double globalQScorePrior = -1.0;
    final boolean emitOriginalQuals = false;

    final RecalibrationReport report = getRecalReport();
    final Covariate[] requestedCovariates = report.getRequestedCovariates();
    // TODO: this covariates need to encode the read groups
    final QuantizationInfo quantizationInfo = report.getQuantizationInfo();
    final RecalibrationTables gatk_tables = report.getRecalibrationTables();

    final int quantizationLevels = 1;

    quantizationInfo.quantizeQualityScores(quantizationLevels);
    final List<Byte> quantizedQuals = quantizationInfo.getQuantizedQuals();

    final BaseRecalibration gatk_engine = new BaseRecalibration(grpPath.toFile(),
                quantizationLevels,
                disableIndelQuals,
                preserveQLessThan,
                emitOriginalQuals,
                globalQScorePrior,
                null, false);

    // use default parameters
    try {
      engine.init(requestedCovariates, gatk_tables,
                quantizedQuals, null,
                disableIndelQuals,
                preserveQLessThan,
                globalQScorePrior,
                emitOriginalQuals);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }

    // NOTE: staticQuantizedMapping is not tested

    final int numReadGroups = gatk_tables.getReadGroupTable().getDimensions()[0];

    final SamReader reader = getInputBamRecords();
    RecalibrationEngine recalibrationEngine = new RecalibrationEngine(requestedCovariates, numReadGroups, RAC.RECAL_TABLE_UPDATE_LOG, false);

    int numRecords = 0;
    for (SAMRecord record : reader) {
      final GATKSAMRecord read = new GATKSAMRecord(record);

      // run recalibrate and then compare the base qualities
      gatk_engine.recalibrateRead(read);

      try {
        final byte[][] quals = engine.recalibrate(read);
        for (final EventType errorModel : EventType.values()) { // recalibrate all three quality strings
          if (disableIndelQuals && errorModel != EventType.BASE_SUBSTITUTION) {
            continue;
          }
          Assert.assertEquals(quals[errorModel.ordinal()], read.getBaseQualities(errorModel));
        }
      }
      catch (AccelerationException e) {
        logger.error("exception caught in init(): "+ e.getMessage());
        return;
      }
      numRecords++;
    }
  }

  @Test(enabled = true, groups = {"bqsr"})
  public void TestFinalize() {
    final Covariate[] covariates = getCovariates();
    final int numReadGroups = 2;
    try {
      engine.init(covariates, numReadGroups);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }
    int numCovariates = covariates.length;
    int numEvents = EventType.values().length;
    int qualLength = covariates[1].maximumKeyValue()+1;

    RecalibrationTables recal_table = engine.getRecalibrationTables();
    engine.finalizeData();
    RecalibrationTables recal_table_1 = engine.getFinalRecalibrationTables();
    RecalibrationTables recal_table_2 = engine.getFinalRecalibrationTables();
  }
  */
  @BeforeMethod
  public void setUp() {
      //TODO
    //engine = new FalconRecalibrationEngine(RAC, helper.getRefReader());
    engine = new FalconRecalibrationEngine(RAC, null);
    final boolean isLoaded = engine.load(null);
    Assert.assertTrue(isLoaded);
  }

  @AfterMethod
  public void tearDown() {
    engine.finalizeData();
    engine = null;
  }
/*
  private final Covariate[] getCovariates() {
    Pair<ArrayList<Covariate>, ArrayList<Covariate>> all_covariates = RecalUtils.initializeCovariates(RAC);
    ArrayList<Covariate> requiredCovariates = all_covariates.getFirst();
    ArrayList<Covariate> optionalCovariates = all_covariates.getSecond();

    final Covariate[] covariates = new Covariate[requiredCovariates.size() + optionalCovariates.size()];
    int cov_idx = 0;
    for (final Covariate covariate : requiredCovariates)
      covariates[cov_idx++] = covariate;

    for (final Covariate covariate : optionalCovariates)
      covariates[cov_idx++] = covariate;

    for (Covariate cov : covariates) { // list all the covariates being used
        cov.initialize(RAC); // initialize any covariate member variables using the shared argument collection
    }
    return covariates;
  }
  */
  private final StandardCovariateList getCovariates() {
      //TODO: second argument needs to be changed
    return new StandardCovariateList(RAC, Collections.singletonList("readGroup"));
  }
  /*
  private final RecalibrationReport getRecalReport() {
    return new RecalibrationReport(grpPath.toFile());
  }

  private final SamReader getInputBamRecords() {

    final SamReaderFactory readerFactory = SamReaderFactory.make();
    readerFactory.validationStringency(ValidationStringency.LENIENT);

    final SamReader reader = readerFactory.open(bamPath);

    return reader;
  }

  private static void compareRecalibrationTables(
        final int numCovariates,
        final RecalibrationTables our_table,
        final RecalibrationTables gatk_table)
  {
    for (int i = 0; i < numCovariates; i++) {
      final List<RecalDatum> gatk_table_contents = gatk_table.getTable(i).getAllValues();
      final List<RecalDatum> our_table_contents = our_table.getTable(i).getAllValues();
      Assert.assertEquals(our_table_contents.size(), gatk_table_contents.size());
      //System.out.println(String.format("%d: %d == %d", i, our_table_contents.size(), gatk_table_contents.size()));
      for (int k = 0; k < gatk_table_contents.size(); k++) {
        //System.out.println(String.format("[%d-O] %d == %d", k, our_table_contents.get(k).getNumObservations(), gatk_table_contents.get(k).getNumObservations()));
        //System.out.println(String.format("[%d-M] %f == %f", k, our_table_contents.get(k).getNumMismatches(), gatk_table_contents.get(k).getNumMismatches()));
        compareRecalDatum(our_table_contents.get(k), gatk_table_contents.get(k));
      }
    }
  }

  private static void compareRecalDatum(final RecalDatum r1, final RecalDatum r2) {
    Assert.assertEquals(r1.getNumObservations(), r2.getNumObservations());
    Assert.assertEquals(r1.getNumMismatches(), r2.getNumMismatches());
    Assert.assertEquals(r1.getEmpiricalErrorRate(), r2.getEmpiricalErrorRate());
    Assert.assertEquals(r1.getEstimatedQReported(), r2.getEstimatedQReported());
    Assert.assertEquals(r1.getEmpiricalQuality(), r2.getEmpiricalQuality());
  }
  */
}
