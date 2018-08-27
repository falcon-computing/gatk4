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

import org.broadinstitute.hellbender.engine.filters.*;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.iterators.SamReaderQueryingIterator;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.recalibration.*;
import org.broadinstitute.hellbender.utils.recalibration.covariates.*;

import org.broadinstitute.hellbender.engine.datasources.*;
import org.broadinstitute.gatk.nativebindings.*;
import org.broadinstitute.hellbender.tools.walkers.bqsr.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceContext.*;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.QualityUtils;

import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMReadGroupRecord;

import org.broadinstitute.hellbender.transformers.BQSRReadTransformer;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.utils.read.ReadUtils;



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
  private TestHelper helper;

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
    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();
    //helper = new TestHelper(refPath.toFile(), RAC, header);
    helper = new TestHelper(refPath, RAC, header);
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
    int numEvents;
    if(RAC.computeIndelBQSRTables){ 
        numEvents = EventType.values().length;
    }
    else{
        numEvents = 1;
    }
    System.out.printf("Peipei Debug, numEvents: %s\n", numEvents);
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
      int covLength = covariates.get(i).maximumKeyValue()+1;
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

    //final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(bases, baseQuals, readLength + "M");

    final SAMRecord read = ArtificialReadUtils.createArtificialSAMRecord(bases, baseQuals, readLength + "M");
    //final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, baseQuals, readLength + "M");
    final GATKRead togatkread = new SAMRecordToGATKReadAdapter(read);
    //read.setBaseQualities(insertionQuals, EventType.BASE_INSERTION);
    //read.setBaseQualities(deletionQuals, EventType.BASE_DELETION);
    togatkread.setBaseQualities(insertionQuals);
    togatkread.setBaseQualities(deletionQuals);

    // set read group to test unsupported platform

    //final GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord("@RG\tID:test\tSM:test\tPL:bgiseq");
    //read.setReadGroup(rg);
    final SAMReadGroupRecord rg = new SAMReadGroupRecord("@RG\tID:test\tSM:test\tPL:bgiseq");
    //final SAMReadGroupRecord rg = new SAMReadGroupRecord("@RG\tID:test\tSM:test\tPL:ILLUMINA");
    togatkread.setReadGroup(rg.getId());

    try {
      final int[][][] falcon_keys = engine.computeCycleCovariates(togatkread);
      //final int[][][] falcon_keys = engine.computeCycleCovariates(togatkread.getEncapsulatedSamRecord());
    }
    catch (AccelerationException e) {
      logger.info("caught exception for unsupported platform");
      return;
    }
    Assert.fail("should have caught exception");
  }



  @Test(enabled = true, groups = {"bqsr"})
  public void TestCycleCovariates() {
    //final Covariate[] covariates = getCovariates();
    //final  StandardCovariateList covariates = getCovariates();
    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();
    final StandardCovariateList covariates = new StandardCovariateList(RAC, header);
    final int numReadGroups = header.getReadGroups().size();
    try {
      engine.init(covariates, numReadGroups);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }

    final CovariateKeyCache keyCache= new CovariateKeyCache();
    for (SAMRecord record : reader) {
      final GATKRead read = new SAMRecordToGATKReadAdapter(record);
      final ReadCovariates cov = RecalUtils.computeCovariates(read, header, covariates, true, keyCache);
      try {
        final int[][][] falcon_keys = engine.computeCycleCovariates(read, header);
        final int contextCovIdx = 3;

        if(RAC.computeIndelBQSRTables){
            for (EventType event : EventType.values()) {
                //System.out.println(Arrays.toString(falcon_keys[event.ordinal()]));
                final int[][] gatk_keys = cov.getKeySet(event);
                //for (int i = 0; i < read.getReadBases().length; i++) {
                for (int i = 0; i < read.getBases().length; i++) {
                    //System.out.println(Arrays.toString(gatk_keys[i]));
                    Assert.assertEquals(
                            falcon_keys[event.ordinal()][i][contextCovIdx],
                            gatk_keys[i][contextCovIdx]);
                }
                }
            }
            else{
                final int[][] gatk_keys = cov.getKeySet(EventType.BASE_SUBSTITUTION);
                for (int i = 0; i < read.getBases().length; i++) {
                    Assert.assertEquals(
                            falcon_keys[0][i][contextCovIdx],
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
    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();
    final StandardCovariateList covariates = new StandardCovariateList(RAC, header);
    final int numReadGroups = header.getReadGroups().size();
    logger.info(String.format("number of read groups: %d", numReadGroups));

    try {
      engine.init(covariates, numReadGroups);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }

    final CovariateKeyCache keyCache= new CovariateKeyCache();
    for (SAMRecord record : reader) {
      final GATKRead read = new SAMRecordToGATKReadAdapter(record);
      final ReadCovariates cov = RecalUtils.computeCovariates(read, header, covariates, true, keyCache);

      final int[][][] falcon_keys = engine.computeContextCovariates(read);

      final int contextCovIdx = 2;
      if(RAC.computeIndelBQSRTables){
          for (EventType event : EventType.values()) {
              final int[][] gatk_keys = cov.getKeySet(event);
              for (int i = 0; i < read.getBases().length; i++) {
                  Assert.assertEquals(
                          falcon_keys[event.ordinal()][i][contextCovIdx],
                          gatk_keys[i][contextCovIdx]);
              }
          }
      }
      else{
          final int[][] gatk_keys0 = cov.getKeySet(EventType.BASE_SUBSTITUTION);
          final int[][] gatk_keys1 = cov.getKeySet(EventType.BASE_INSERTION);
          final int[][] gatk_keys2 = cov.getKeySet(EventType.BASE_DELETION);
          for (int i = 0; i < read.getBases().length; i++) {
              Assert.assertEquals(
                      falcon_keys[0][i][contextCovIdx],
                      gatk_keys0[i][contextCovIdx]);
          }
      }
    }
  }

  @Test(enabled = true, groups = {"bqsr"})
  public void TestCovariates() {
    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();
    final StandardCovariateList covariates = new StandardCovariateList(RAC, header);
    final int numReadGroups = header.getReadGroups().size();

    try {
      engine.init(covariates, numReadGroups, header);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }

    int numCovariates = covariates.size();

    final CovariateKeyCache keyCache= new CovariateKeyCache();
    for (SAMRecord record : reader) {
      final GATKRead read = new SAMRecordToGATKReadAdapter(record);
      final ReadCovariates cov = RecalUtils.computeCovariates(read, header, covariates, true, keyCache);
      final int[] falcon_keys = engine.computeCovariates(read, header);
      int readLength = read.getBases().length;

      int idx = 0;
      for (int i = 0; i < readLength; i++) {
          for (int j = 0; j < numCovariates; j++) {
              if(RAC.computeIndelBQSRTables){
                  for (EventType event : EventType.values()) {
                      final int gatk_key = cov.getKeySet(event)[i][j];
                      Assert.assertEquals(falcon_keys[idx++], gatk_key);
                  }
              }
              else{
                  final int gatk_key = cov.getKeySet(EventType.BASE_SUBSTITUTION)[i][j];
                  //if(falcon_keys[idx]!=gatk_key){
                  //    System.out.printf("Peipei Debug i : %d, j: %d \n", i, j);
                  //}
                  //idx++;
                  Assert.assertEquals(falcon_keys[idx++], gatk_key);
              }
          }
      }
    }
  }

  @Test(enabled = true, groups = {"bqsr"})
  public void TestReadGroupCovariate() {

    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();

    final StandardCovariateList covariates = new StandardCovariateList(RAC, header);
    final StandardCovariateList gatk_covariates = new StandardCovariateList(RAC, header);

    final int numReadGroups = header.getReadGroups().size();
    final int numCovariates = covariates.size();

    try {
      engine.init(covariates, numReadGroups);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }

    int numRecords = 0;
    final CovariateKeyCache keyCache= new CovariateKeyCache();
    for (SAMRecord record : reader) {
      final GATKRead read = new SAMRecordToGATKReadAdapter(record);
      final ReadCovariates cov = RecalUtils.computeCovariates(read, header, gatk_covariates, true, keyCache);
      final int[] falcon_keys = engine.computeCovariates(read, header);
    }
    engine.updateReadGroupCovariates();

    // check read group covariates
    final ReadGroupCovariate falcon_rgCov = (ReadGroupCovariate)covariates.get(0);
    final ReadGroupCovariate gatk_rgCov = (ReadGroupCovariate)gatk_covariates.get(0);
    Assert.assertEquals(gatk_rgCov.maximumKeyValue(), falcon_rgCov.maximumKeyValue());

    for (int i = 0; i < gatk_rgCov.maximumKeyValue(); i++) {
      final String gatk_rg = gatk_rgCov.formatKey(i);
      final String falcon_rg = falcon_rgCov.formatKey(i);
      Assert.assertEquals(gatk_rg, falcon_rg);
    }
  }

  @Test(enabled = true, groups = {"bqsr"})
  public void TestBAQCalculation() {
    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();
    final StandardCovariateList covariates = new StandardCovariateList(RAC, header);
    final int numReadGroups = header.getReadGroups().size();
    final int numCovariates = covariates.size();

    try {
      engine.init(covariates, numReadGroups, header);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }

    //BaseRecalibrationEngine recalibrationEngine = new BaseRecalibrationEngine(covariates, numReadGroups, RAC.RECAL_TABLE_UPDATE_LOG, false);

    int numRecords = 0;
    for (SAMRecord record : reader) {
      //if (numRecords > 1) break;
      final GATKRead org_read = new SAMRecordToGATKReadAdapter(record);
      final GATKRead read = ReadClipper.hardClipSoftClippedBases(ReadClipper.hardClipAdaptorSequence(org_read));
      //final ReferenceContext ref = helper.getRefContext(org_read);
      int[] isSNP = new int[read.getLength()];
      int[] isInsertion = new int[isSNP.length];
      int[] isDeletion = new int[isSNP.length];

      BaseRecalibrationEngine br = new BaseRecalibrationEngine(RAC, header);

      final int nErrors = BaseRecalibrationEngine.calculateIsSNPOrIndel(read, helper.getRefDataSource(), isSNP, isInsertion, isDeletion);


      if (nErrors == 0) continue;

      final byte[] gatk_baqArray = helper.falconCalculateBAQArray(read);
      final byte[] falcon_baqArray = engine.calculateBAQArray(read, helper.getRefDataSource());

      if (gatk_baqArray == null) {
        Assert.assertEquals(falcon_baqArray, null);
      }
      else {
        Assert.assertEquals(falcon_baqArray.length, gatk_baqArray.length);
        for (int i = 0; i < gatk_baqArray.length; i++) {
          Assert.assertEquals(falcon_baqArray[i], gatk_baqArray[i]);
        }
      }
      numRecords++;
    }
  }


  @Test(enabled = true, groups = {"bqsr"})
  public void TestFractionalErrors() {
    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();
    final StandardCovariateList covariates = new StandardCovariateList(RAC, header);
    final int numReadGroups = header.getReadGroups().size();
    final int numCovariates = covariates.size();

    try {
      engine.init(covariates, numReadGroups, header);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }

    //RecalibrationEngine recalibrationEngine = new RecalibrationEngine(covariates, numReadGroups, RAC.RECAL_TABLE_UPDATE_LOG, false);

    int numRecords = 0;
    for (SAMRecord record : reader) {
      //if (numRecords > 1) break;
         final GATKRead org_read = new SAMRecordToGATKReadAdapter(record);
      final GATKRead read = ReadClipper.hardClipSoftClippedBases(ReadClipper.hardClipAdaptorSequence(org_read));

      int readLength = read.getLength();
      int[] isSNP = new int[read.getLength()];
      int[] isInsertion = new int[isSNP.length];
      int[] isDeletion = new int[isSNP.length];

      BaseRecalibrationEngine br = new BaseRecalibrationEngine(RAC, header);

      final int nErrors = BaseRecalibrationEngine.calculateIsSNPOrIndel(read, helper.getRefDataSource(), isSNP, isInsertion, isDeletion);

      //final byte[] baqArray = nErrors == 0 ? helper.falconFlatBAQArray(read) : helper.falconCalculateBAQArray(read);
      final boolean enableBAQ = RAC.enableBAQ;
      final byte[] baqArray = (nErrors == 0 || !enableBAQ) ? helper.falconFlatBAQArray(read) : helper.falconCalculateBAQArray(read);



      final double[][] errors = engine.calculateFractionalErrorArray(read, org_read, helper.getRefDataSource());


      if (baqArray == null) { // some reads just can't be BAQ'ed
        Assert.assertEquals(null, errors);
      }
      else {

        //final boolean[] skip = calculateSkipArray(read, metaDataTracker); // skip known sites of variation as well as low quality and non-regular bases
        final double[] snpErrors = helper.falconCalculateFractionalErrorArray(isSNP, baqArray);
        final double[] insertionErrors = helper.falconCalculateFractionalErrorArray(isInsertion, baqArray);
        final double[] deletionErrors = helper.falconCalculateFractionalErrorArray(isDeletion, baqArray);


        Assert.assertNotNull(errors);
        if(RAC.computeIndelBQSRTables){
            for (int i = 0; i < readLength; i++) {
                Assert.assertEquals(snpErrors[i], errors[0][i]);
                Assert.assertEquals(insertionErrors[i], errors[1][i]);
                Assert.assertEquals(deletionErrors[i], errors[2][i]);
            }
        }
        else{
            for (int i = 0; i < readLength; i++) {
                Assert.assertEquals(snpErrors[i], errors[0][i]);
                //Assert.assertEquals(insertionErrors[i], errors[1][i]);
                //Assert.assertEquals(deletionErrors[i], errors[2][i]);
            }
        }
      }

      //System.out.println(String.format("finish read %d", numRecords));
      numRecords++;
    }
  }

  @Test(enabled = true, groups = {"bqsr"})
  public void TestFinalize() {
    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();
    final StandardCovariateList covariates = new StandardCovariateList(RAC, header);
    final int numReadGroups = 2;
    try {
      //engine.init(covariates, numReadGroups);
      engine.init(covariates, numReadGroups, header);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }
    int numCovariates = covariates.size();
    int numEvents = EventType.values().length;
    int qualLength = covariates.get(1).maximumKeyValue()+1;


    RecalibrationTables recal_table = engine.getRecalibrationTables();

    System.out.println("my test finalize");
    for (int i = 0; i < numCovariates; i++) {
      //List<RecalDatum> gatk_table_contents = recal_table.getTable(i).getAllValues();
      List<RecalDatum> our_table_contents = recal_table.getTable(i).getAllValues();
      //if (our_table_contents.size() != gatk_table_contents.size()) {
      System.out.printf("%d: ours: %d\n", i, our_table_contents.size());
      //}
    }
    engine.finalizeData();
    RecalibrationTables recal_table_1 = engine.getFinalRecalibrationTables();
    RecalibrationTables recal_table_2 = engine.getFinalRecalibrationTables();
    //RecalibrationTables our_table = engine.getFinalRecalibrationTables();

    for (int i = 0; i < numCovariates; i++) {
      //List<RecalDatum> gatk_table_contents = recal_table.getTable(i).getAllValues();
      List<RecalDatum> our_table_contents = recal_table_1.getTable(i).getAllValues();
      //if (our_table_contents.size() != gatk_table_contents.size()) {
        System.out.printf("%d: ours: %d\n", i, our_table_contents.size());
      //}
    }
  }



  @Test(enabled = true, groups = {"bqsr"})
  public void TestTableUpdateWithRealData() {
    //RAC.computeIndelBQSRTables = true;
    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();
    final StandardCovariateList covariates = new StandardCovariateList(RAC, header);
    final int numReadGroups = header.getReadGroups().size();
    final int numCovariates = covariates.size();

    try {
      //engine.init(covariates, numReadGroups);
      engine.init(covariates, numReadGroups, header);
    }
    catch (AccelerationException e) {
      logger.error("exception caught in init(): "+ e.getMessage());
      return;
    }

    //RecalibrationEngine recalibrationEngine = new RecalibrationEngine(covariates, numReadGroups, RAC.RECAL_TABLE_UPDATE_LOG, false);
    BaseRecalibrationEngine recalibrationEngine = new BaseRecalibrationEngine(RAC, header);

    int numRecords = 0;
    for (SAMRecord record : reader) {
      final GATKRead org_read = new SAMRecordToGATKReadAdapter(record);
      //final GATKRead read = ReadClipper.hardClipSoftClippedBases(ReadClipper.hardClipAdaptorSequence(org_read));
      final ReadTransformer transform = recalibrationEngine.makeReadTransform();
      final GATKRead read = transform.apply(org_read);
      RecalUtils.parsePlatformForRead(read, header, RAC);

      //final ReferenceContext ref = helper.getRefContext(org_read);

      //final int readLength = read.getReadBases().length;
      final int readLength = read.getLength();
      final boolean[] skip = new boolean[readLength];
      Arrays.fill(skip, false);

      //System.out.println("@@@ before update");
      //for (int i = 0; i < 4; i++){
          //System.out.println(engine.getDebugTable().getTable(i).getAllValues().size());
      //}
      // perform falcon table update
      int ret = 0;
      try {
        //ret = engine.update(read, org_read, ref, skip);
        ret = engine.update(read, org_read, helper.getRefDataSource(), header, skip);
      }
      catch (AccelerationException e) {
        logger.error("exception caught in init(): "+ e.getMessage());
        return;
      }
      //

      int[] isSNP = new int[read.getLength()];
      int[] isInsertion = new int[isSNP.length];
      int[] isDeletion = new int[isSNP.length];
      final int nErrors = BaseRecalibrationEngine.calculateIsSNPOrIndel(read, helper.getRefDataSource(), isSNP, isInsertion, isDeletion);
      //final byte[] baqArray = nErrors == 0 ? helper.falconFlatBAQArray(read) : helper.falconCalculateBAQArray(read);
      final boolean enableBAQ = RAC.enableBAQ;
      final byte[] baqArray = (nErrors == 0 || !enableBAQ) ? helper.falconFlatBAQArray(read) : helper.falconCalculateBAQArray(read);

      if (baqArray != null) { // some reads just can't be BAQ'ed
        //final boolean[] skip = calculateSkipArray(read, metaDataTracker); // skip known sites of variation as well as low quality and non-regular bases

        final double[] snpErrors = helper.falconCalculateFractionalErrorArray(isSNP, baqArray);
        final double[] insertionErrors = helper.falconCalculateFractionalErrorArray(isInsertion, baqArray);
        final double[] deletionErrors = helper.falconCalculateFractionalErrorArray(isDeletion, baqArray);

        //final ReadCovariates cov = RecalUtils.computeCovariates(read, covariates);
        final CovariateKeyCache keyCache= new CovariateKeyCache();
        final ReadCovariates cov = RecalUtils.computeCovariates(read, header, covariates, true, keyCache);

        // aggregate all of the info into our info object, and update the data
        final ReadRecalibrationInfo info = new ReadRecalibrationInfo(read, cov, skip, snpErrors, insertionErrors, deletionErrors);
        //recalibrationEngine.updateDataForRead(info);
        recalibrationEngine.updateRecalTablesForRead(info);
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

   
    //System.out.println("test gatk score");
    int counter=0;
    final NestedIntegerArray<RecalDatum> byQualTable = gatk_table.getQualityScoreTable();
    for ( final NestedIntegerArray.Leaf<RecalDatum> leaf : byQualTable.getAllLeaves() ) {

        final int rgKey = leaf.keys[0];
        final int eventIndex = leaf.keys[2];
        final RecalDatum qualDatum = leaf.value;
        // create a copy of qualDatum, and initialize byReadGroup table with it
        //System.out.printf("null branch  rgKey: %d, eventIndex: %d , qualDatum: %s\n", rgKey, eventIndex, qualDatum.toString());
        counter+=1;
    }    
    //System.out.printf("@@@ counter is : %d\n",counter);


    // compare all tables
    //
    for (int i = 0; i < numCovariates; i++) {
      List<RecalDatum> gatk_table_contents = gatk_table.getTable(i).getAllValues();
      List<RecalDatum> our_table_contents = our_table.getTable(i).getAllValues();
        //System.out.printf("%d: gatk: %d, ours: %d\n", i, gatk_table_contents.size(), our_table_contents.size());
      //// starts
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
    //final Covariate[] requestedCovariates = report.getRequestedCovariates();
    StandardCovariateList requestedCovariates = report.getCovariates();
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

    compareRecalibrationTables(requestedCovariates.size(), falcon_tables, gatk_tables);
  }


  @Test(enabled = true, groups = {"pr"})
  public void TestRecalibrate() {
    //final boolean disableIndelQuals = false;
    final boolean disableIndelQuals = true;
    final int preserveQLessThan = QualityUtils.MIN_USABLE_Q_SCORE;
    final double globalQScorePrior = -1.0;
    final boolean emitOriginalQuals = false;
    final ApplyBQSRArgumentCollection bqsrArgs = new ApplyBQSRArgumentCollection();

    final RecalibrationReport report = getRecalReport();
    //final Covariate[] requestedCovariates = report.getRequestedCovariates();
    StandardCovariateList requestedCovariates = report.getCovariates();
    // TODO: this covariates need to encode the read groups
    final QuantizationInfo quantizationInfo = report.getQuantizationInfo();
    final RecalibrationTables gatk_tables = report.getRecalibrationTables();

    //final int quantizationLevels = 1;

    //quantizationInfo.quantizeQualityScores(quantizationLevels);
    if (bqsrArgs.quantizationLevels == 0) { // quantizationLevels == 0 means no quantization, preserve the quality scores
        quantizationInfo.noQuantization();
    } else if (bqsrArgs.quantizationLevels > 0 && bqsrArgs.quantizationLevels != quantizationInfo.getQuantizationLevels()) { // any other positive value means, we want a different quantization than the one pre-calculated in the recalibration report. Negative values mean the user did not provide a quantization argument, and just wants to use what's in the report.
        quantizationInfo.quantizeQualityScores(bqsrArgs.quantizationLevels);
    }


    
    final List<Byte> quantizedQuals = quantizationInfo.getQuantizedQuals();
    System.out.printf("ingatk quantizedQuals size is %d, array is %s\n", quantizedQuals.size(), Arrays.toString(quantizedQuals.toArray()));

    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();

    final BQSRReadTransformer gatk_engine = new BQSRReadTransformer(header, grpPath.toFile(), bqsrArgs);


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

    //RecalibrationEngine recalibrationEngine = new RecalibrationEngine(requestedCovariates, numReadGroups, RAC.RECAL_TABLE_UPDATE_LOG, false);

    int numRecords = 0;
    final SamReader reader1 = getInputBamRecords();
    for (SAMRecord record : reader1) {
      final GATKRead originalRead1 = new SAMRecordToGATKReadAdapter(record);
      //bqsrArgs.useOriginalBaseQualities = true;
      //bqsrArgs.useOriginalBaseQualities default is false;
      final GATKRead toFalconread = bqsrArgs.useOriginalBaseQualities ? ReadUtils.resetOriginalBaseQualities(originalRead1) : originalRead1;
      //System.out.printf("bqsrArgs.useOriginalBaseQualities is : ");
      //System.out.println(bqsrArgs.useOriginalBaseQualities);
      //System.out.printf("read num: %d\n", numRecords);
      //System.out.printf("before falc apply, readnum: %d\n", numRecords);
      //System.out.printf("before falc apply: read %s\n", Arrays.toString(toFalconread.getBaseQualities()));
      try {
        final byte[][] quals = engine.recalibrate(toFalconread, header);
        //GATKRead new_gatkRead = gatk_engine.apply(toFalconread); // new_gatkRead is the same as toFalconread as toFalconread is returned
        gatk_engine.apply(toFalconread);
        for (final EventType errorModel : EventType.values()) { // recalibrate all three quality strings
          if (disableIndelQuals && errorModel != EventType.BASE_SUBSTITUTION) {
            continue;
          }
          //System.out.printf("after  falc apply: read qual[] %s\n", Arrays.toString(quals[errorModel.ordinal()]));
          //System.out.printf("after  gatk apply: read Base() %s\n", Arrays.toString(toFalconread.getBaseQualities()));
          Assert.assertEquals(quals[errorModel.ordinal()], toFalconread.getBaseQualities());
        }
      }
      catch (AccelerationException e) {
        logger.error("exception caught in init(): "+ e.getMessage());
        return;
      }
      numRecords++;
    }
  }

  @Test(enabled = true, groups = {"pr"})
  public void TestRecalibrateUseOriginalBaseQualities() {
    //final boolean disableIndelQuals = false;
    final boolean disableIndelQuals = true;
    final int preserveQLessThan = QualityUtils.MIN_USABLE_Q_SCORE;
    final double globalQScorePrior = -1.0;
    final boolean emitOriginalQuals = false;
    final ApplyBQSRArgumentCollection bqsrArgs = new ApplyBQSRArgumentCollection();
    bqsrArgs.useOriginalBaseQualities = true;
    System.out.printf("bqsrArgs.useOriginalBaseQualities is : ");
    System.out.println(bqsrArgs.useOriginalBaseQualities);

    final RecalibrationReport report = getRecalReport();
    //final Covariate[] requestedCovariates = report.getRequestedCovariates();
    StandardCovariateList requestedCovariates = report.getCovariates();
    // TODO: this covariates need to encode the read groups
    final QuantizationInfo quantizationInfo = report.getQuantizationInfo();
    final RecalibrationTables gatk_tables = report.getRecalibrationTables();

    //final int quantizationLevels = 1;

    //quantizationInfo.quantizeQualityScores(quantizationLevels);
    if (bqsrArgs.quantizationLevels == 0) { // quantizationLevels == 0 means no quantization, preserve the quality scores
        quantizationInfo.noQuantization();
    } else if (bqsrArgs.quantizationLevels > 0 && bqsrArgs.quantizationLevels != quantizationInfo.getQuantizationLevels()) { // any other positive value means, we want a different quantization than the one pre-calculated in the recalibration report. Negative values mean the user did not provide a quantization argument, and just wants to use what's in the report.
        quantizationInfo.quantizeQualityScores(bqsrArgs.quantizationLevels);
    }


    
    final List<Byte> quantizedQuals = quantizationInfo.getQuantizedQuals();
    System.out.printf("ingatk quantizedQuals size is %d, array is %s\n", quantizedQuals.size(), Arrays.toString(quantizedQuals.toArray()));

    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();

    final BQSRReadTransformer gatk_engine = new BQSRReadTransformer(header, grpPath.toFile(), bqsrArgs);




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

    //RecalibrationEngine recalibrationEngine = new RecalibrationEngine(requestedCovariates, numReadGroups, RAC.RECAL_TABLE_UPDATE_LOG, false);

    int numRecords = 0;
    final SamReader reader1 = getInputBamRecords();
    for (SAMRecord record : reader1) {
      final GATKRead originalRead1 = new SAMRecordToGATKReadAdapter(record);
      final GATKRead toFalconread = bqsrArgs.useOriginalBaseQualities ? ReadUtils.resetOriginalBaseQualities(originalRead1) : originalRead1;
      try {
        final byte[][] quals = engine.recalibrate(toFalconread, header);
        //GATKRead new_gatkRead = gatk_engine.apply(toFalconread); // new_gatkRead is the same as toFalconread as toFalconread is returned
        gatk_engine.apply(toFalconread);
        for (final EventType errorModel : EventType.values()) { // recalibrate all three quality strings
          if (disableIndelQuals && errorModel != EventType.BASE_SUBSTITUTION) {
            continue;
          }
          //System.out.printf("after  falc apply: read qual[] %s\n", Arrays.toString(quals[errorModel.ordinal()]));
          //System.out.printf("after  gatk apply: read Base() %s\n", Arrays.toString(toFalconread.getBaseQualities()));
          Assert.assertEquals(quals[errorModel.ordinal()], toFalconread.getBaseQualities());
        }
      }
      catch (AccelerationException e) {
        logger.error("exception caught in init(): "+ e.getMessage());
        return;
      }
      numRecords++;
    }
  }

  @Test(enabled = true, groups = {"pr"})
  public void TestRecalibrateWithSQQ() {
    final int hey;
    //final boolean disableIndelQuals = false;
    final boolean disableIndelQuals = true;
    final int preserveQLessThan = QualityUtils.MIN_USABLE_Q_SCORE;
    final double globalQScorePrior = -1.0;
    final boolean emitOriginalQuals = false;
    final ApplyBQSRArgumentCollection bqsrArgs = new ApplyBQSRArgumentCollection();

    final RecalibrationReport report = getRecalReport();
    //final Covariate[] requestedCovariates = report.getRequestedCovariates();
    StandardCovariateList requestedCovariates = report.getCovariates();
    // TODO: this covariates need to encode the read groups
    final QuantizationInfo quantizationInfo = report.getQuantizationInfo();
    final RecalibrationTables gatk_tables = report.getRecalibrationTables();

    bqsrArgs.staticQuantizationQuals.add(10);
    bqsrArgs.staticQuantizationQuals.add(20);
    bqsrArgs.staticQuantizationQuals.add(30);

    byte[] staticQuantizedMapping;
    if(bqsrArgs.staticQuantizationQuals != null && !bqsrArgs.staticQuantizationQuals.isEmpty()) {
      staticQuantizedMapping = BQSRReadTransformer.constructStaticQuantizedMapping(bqsrArgs.staticQuantizationQuals, bqsrArgs.roundDown);
    }
    else
    {
      staticQuantizedMapping = null;
    }
    //final int quantizationLevels = 1;

    //quantizationInfo.quantizeQualityScores(quantizationLevels);
    if (bqsrArgs.quantizationLevels == 0) { // quantizationLevels == 0 means no quantization, preserve the quality scores
      quantizationInfo.noQuantization();
    } else if (bqsrArgs.quantizationLevels > 0 && bqsrArgs.quantizationLevels != quantizationInfo.getQuantizationLevels()) { // any other positive value means, we want a different quantization than the one pre-calculated in the recalibration report. Negative values mean the user did not provide a quantization argument, and just wants to use what's in the report.
      quantizationInfo.quantizeQualityScores(bqsrArgs.quantizationLevels);
    }


    //
    final List<Byte> quantizedQuals = quantizationInfo.getQuantizedQuals();
    System.out.printf("ingatk quantizedQuals size is %d, array is %s\n", quantizedQuals.size(), Arrays.toString(quantizedQuals.toArray()));

    final SamReader reader = getInputBamRecords();
    final SAMFileHeader header = reader.getFileHeader();

    final BQSRReadTransformer gatk_engine = new BQSRReadTransformer(header, grpPath.toFile(), bqsrArgs);



    // use default parameters
    try {
      engine.init(requestedCovariates, gatk_tables,
              quantizedQuals, staticQuantizedMapping,
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

    //RecalibrationEngine recalibrationEngine = new RecalibrationEngine(requestedCovariates, numReadGroups, RAC.RECAL_TABLE_UPDATE_LOG, false);

    int numRecords = 0;
    final SamReader reader1 = getInputBamRecords();
    for (SAMRecord record : reader1) {
      final GATKRead originalRead1 = new SAMRecordToGATKReadAdapter(record);
      final GATKRead toFalconread = bqsrArgs.useOriginalBaseQualities ? ReadUtils.resetOriginalBaseQualities(originalRead1) : originalRead1;

      //System.out.printf("read num: %d\n", numRecords);
      //System.out.printf("before falc apply, readnum: %d\n", numRecords);
      //System.out.printf("before falc apply: read %s\n", Arrays.toString(toFalconread.getBaseQualities()));
      try {
        final byte[][] quals = engine.recalibrate(toFalconread, header);
        //GATKRead new_gatkRead = gatk_engine.apply(toFalconread); // new_gatkRead is the same as toFalconread as toFalconread is returned
        gatk_engine.apply(toFalconread);
        for (final EventType errorModel : EventType.values()) { // recalibrate all three quality strings
          if (disableIndelQuals && errorModel != EventType.BASE_SUBSTITUTION) {
            continue;
          }
          //System.out.printf("after  falc apply: read qual[] %s\n", Arrays.toString(quals[errorModel.ordinal()]));
          //System.out.printf("after  gatk apply: read Base() %s\n", Arrays.toString(toFalconread.getBaseQualities()));
          Assert.assertEquals(quals[errorModel.ordinal()], toFalconread.getBaseQualities());
        }
      }
      catch (AccelerationException e) {
        logger.error("exception caught in init(): "+ e.getMessage());
        return;
      }
      numRecords++;
    }
  }


  @BeforeMethod
  public void setUp() {
    //TODO FalconRecalibrationEngine second argument
    //engine = new FalconRecalibrationEngine(RAC, helper.getRefReader());
    //RAC.computeIndelBQSRTables = false;
    int testFlag=0;
    switch (testFlag){
        case 0:
            RAC.computeIndelBQSRTables = false;
            RAC.enableBAQ = false;
            break;

        case 1:
            RAC.computeIndelBQSRTables = false;
            RAC.enableBAQ = true;
            break;

        case 2:
            RAC.computeIndelBQSRTables = true;
            RAC.enableBAQ = false;
            break;

        case 3:
            RAC.computeIndelBQSRTables = true;
            RAC.enableBAQ = true;
            break;


    }
    //RAC.computeIndelBQSRTables = false;
    //RAC.enableBAQ = false;
    //RAC.enableBAQ = true;
    engine = new FalconRecalibrationEngine(RAC, null);
    //engine = new FalconRecalibrationEngine(RAC, helper.getRefReader());
    final boolean isLoaded = engine.load(null);
    Assert.assertTrue(isLoaded);
  }

  @AfterMethod
  public void tearDown() {
    engine.finalizeData();
    engine = null;
  }


  private final StandardCovariateList getCovariates() {
      //TODO: second argument needs to be changed
    return new StandardCovariateList(RAC, Collections.singletonList("readGroup"));
  }

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

}
