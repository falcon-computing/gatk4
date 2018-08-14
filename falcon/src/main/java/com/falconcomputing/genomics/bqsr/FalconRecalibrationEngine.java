/**
 * Copyright (C) 2015-2017 Falcon Computing Solutions, Inc.
 * All Rights Reserved.
 * All information contained herein is considered trade secrets and
 * confidential. Dissemination or reproduction, in whole or in part, is
 * strictly forbidden unless A prior written permission is obtained from
 * Falcon Computing Solutions.
 */

package com.falconcomputing.genomics.bqsr;

import com.falconcomputing.genomics.NativeLibraryLoader;
import com.falconcomputing.genomics.AccelerationException;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

import org.apache.log4j.Logger;

import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.covariates.*;
import org.broadinstitute.hellbender.utils.recalibration.*;
import org.broadinstitute.gatk.nativebindings.NativeLibrary;
import org.broadinstitute.hellbender.utils.baq.BAQ;
import org.broadinstitute.hellbender.utils.collections.NestedIntegerArray;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.NGSPlatform;
import org.broadinstitute.hellbender.utils.recalibration.EventType;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMReadGroupRecord;
//import org.broadinstitute.hellbender.utils.sam.GATKSAMReadGroupRecord;
//import org.broadinstitute.hellbender.utils.sam.GATKSAMRecord;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.SequencerFlowClass;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import htsjdk.samtools.SAMTag;

import java.io.PrintStream;
import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.lang.Exception;
import htsjdk.samtools.SAMFileHeader;

public class FalconRecalibrationEngine implements NativeLibrary {
  private final static Logger logger = Logger.getLogger(FalconRecalibrationEngine.class);
  private static final String NATIVE_LIBRARY_NAME = "falcon_genomics";
  private String nativeLibraryName = "falcon_genomics";
  boolean useFPGA = false;

  private final int numEvents = EventType.values().length;
  private int numCovariates;
  private int numReadGroups;
  //private Covariate[] covariates;
  private StandardCovariateList covariates;
  private RecalibrationTables recalTables;
  private static boolean loaded = false;
  private boolean initialized = false;
  private boolean finalized = false;

  // parameters for ContextCovariate
  private final byte LOW_QUAL_TAIL;
  private final int MISMATCHES_CONTEXT_SIZE;
  private final int INDELS_CONTEXT_SIZE;
  // parameters for CycleCovariate
  private final int MAXIMUM_CYCLE_VALUE;
  private final int CUSHION_FOR_INDEL = 4;
  //private final String FORCE_READGROUP;

  // parameters for BAQ calculation
  private double BAQGOP = BAQ.DEFAULT_GOP; // maybe set by outside?

  // TODO: here may need to allow parameters to be set
  //       for now don't bother sending it to native
  private boolean includeClippedBases = false;
  private static final byte minBaseQual = 4;
  private static final String BAQ_TAG = "BQ";

  private BAQ baq; // BAQ the reads on the fly to generate the alignment uncertainty vector
  private final ReferenceSequenceFile referenceReader; // fasta reference reader for use with BAQ calculation

  public FalconRecalibrationEngine(final RecalibrationArgumentCollection RAC, final ReferenceSequenceFile referenceReader) {
    this.LOW_QUAL_TAIL = RAC.LOW_QUAL_TAIL;
    this.MISMATCHES_CONTEXT_SIZE = RAC.MISMATCHES_CONTEXT_SIZE;
    this.INDELS_CONTEXT_SIZE = RAC.INDELS_CONTEXT_SIZE;
    this.MAXIMUM_CYCLE_VALUE = RAC.MAXIMUM_CYCLE_VALUE;

    //this.FORCE_READGROUP = RAC.FORCE_READGROUP;
    this.referenceReader = referenceReader;

    logger.debug("Created one instance of FalconRecalibrationEngine");
  }

  public synchronized boolean load(File tempDir) {
    if (!NativeLibraryLoader.load(tempDir, NATIVE_LIBRARY_NAME)) {
      logger.warn("Falcon Genomics Acceleration Library is "+
                  "not loaded successfully, switching back to "+
                  "the original implementation.");
      return false;
    }
    if (!loaded) {
      loaded = true;
      logger.info("Falcon Genomics Acceleration Library is loaded.");
    }
    return true;
  }

  // This init() function is used for GATK BaseRecalibrator
  public void init(final StandardCovariateList _covariates,
  //public void init(final Covariate[] _covariates,
                   final int _numReadGroups) throws AccelerationException {

    this.covariates = _covariates;
    this.numCovariates = covariates.size();
    this.numReadGroups = _numReadGroups;

    this.baq = new BAQ(BAQGOP); // setup the BAQ object with the provided gap open penalty

    final int[] covariatesDimensions = new int[numCovariates];

    for (int i = 0; i < covariates.size(); i++) {
      covariatesDimensions[i] = covariates.get(i).maximumKeyValue() + 1;
    }

    // call native method to initialize the covariatesTable
    // TODO: need to make sure thread-safety
    initNative(numReadGroups, numEvents, numCovariates,
               covariatesDimensions,
               LOW_QUAL_TAIL,
               MISMATCHES_CONTEXT_SIZE,
               INDELS_CONTEXT_SIZE,
               MAXIMUM_CYCLE_VALUE,
               CUSHION_FOR_INDEL);

    recalTables = new RecalibrationTables(covariates, numReadGroups);

    if (!initialized) {
      initialized = true;
    }
    logger.debug("Initialized FalconRecalibrationEngine");
  }

  // This init() function is used for GATK BaseRecalibrator with header
  public void init(final StandardCovariateList _covariates,
                   //public void init(final Covariate[] _covariates,
                   final int _numReadGroups, final SAMFileHeader header) throws AccelerationException {

    this.covariates = _covariates;
    this.numCovariates = covariates.size();
    this.numReadGroups = _numReadGroups;
    this.baq = new BAQ(BAQGOP); // setup the BAQ object with the provided gap open penalty

    final int[] covariatesDimensions = new int[numCovariates];

    for (int i = 0; i < covariates.size(); i++) {
      covariatesDimensions[i] = covariates.get(i).maximumKeyValue() + 1;
    }

    // call native method to initialize the covariatesTable
    // TODO: need to make sure thread-safety
    initNative(numReadGroups, numEvents, numCovariates,
            covariatesDimensions,
            LOW_QUAL_TAIL,
            MISMATCHES_CONTEXT_SIZE,
            INDELS_CONTEXT_SIZE,
            MAXIMUM_CYCLE_VALUE,
            CUSHION_FOR_INDEL);

    //after initNative, call bqsrInitReadGroupNative for readGroupCovariate
    final List<String> allReadGroups = ReadGroupCovariate.getReadGroupIDs(header);
    allReadGroups.forEach(
            readGroupId -> {
              //System.out.println(readGroupId);
              bqsrInitReadGroupNative(readGroupId);
            }
    );

    recalTables = new RecalibrationTables(covariates, numReadGroups);

    if (!initialized) {
      initialized = true;
    }
    logger.debug("Initialized FalconRecalibrationEngine");
  }



  // This init() function is used for GATK PrintReads
  public void init(final StandardCovariateList _covariates,
                   final RecalibrationTables _recalTables,
                   final List<Byte> quantizedQuals,
                   final byte[] staticQuantizedMapping,
                   final boolean disableIndelQuals,
                   final int preserveQLessThan,
                   final double globalQScorePrior,
                   final boolean emitOriginalQuals)
    throws AccelerationException
  {
    this.covariates = _covariates;
    this.numCovariates = covariates.size();
    this.recalTables = _recalTables;

    this.numReadGroups = _recalTables.getReadGroupTable().getDimensions()[0];

    final byte[] quantizationTable = new byte[quantizedQuals.size()];
    for (int i = 0; i < quantizedQuals.size(); i++) {
      quantizationTable[i] = quantizedQuals.get(i);
    }
    final Covariate[] toNative = new Covariate[_covariates.size()];
    for(int i = 0; i < _covariates.size(); i++){
      toNative[i] = _covariates.get(i);
    }
    // call native method to initialize the RecalibrationTable
    initNative(numEvents,
               //_covariates,
               toNative,
               quantizationTable,
               staticQuantizedMapping,
               disableIndelQuals,
               preserveQLessThan,
               globalQScorePrior,
               emitOriginalQuals,
               LOW_QUAL_TAIL,
               MISMATCHES_CONTEXT_SIZE,
               INDELS_CONTEXT_SIZE,
               MAXIMUM_CYCLE_VALUE,
               CUSHION_FOR_INDEL);

    // initialize table contents
    for (int i = 0; i < numCovariates; i++) {
      final NestedIntegerArray<RecalDatum> table = _recalTables.getTable(i);
      for (final NestedIntegerArray.Leaf<RecalDatum> leaf : table.getAllLeaves()) {
        setTableNative(leaf.value.getNumObservations(),
            leaf.value.getNumMismatches(),
            leaf.value.getEstimatedQReported(),
            leaf.keys, i);
      }
    }

    if (!initialized) {
      initialized = true;
    }
    logger.debug("Initialized FalconRecalibrationEngine with RecalibrationTables");
  }

  // This function is for unit testing only
  // protected int[][][] computeContextCovariates(final GATKSAMRecord read) {
  //protected int[][][] computeContextCovariates(final SAMRecord read) {
  protected int[][][] computeContextCovariates(final GATKRead read) {

    //final byte[] bases = read.getReadBases();
    final byte[] bases = read.getBases();
    final byte[] quals = read.getBaseQualities();
    //final boolean isNegativeStrand = read.getReadNegativeStrandFlag();
    final boolean isNegativeStrand = read.isReverseStrand();

    final int[] keys = computeContextCovariatesNative(bases, quals, isNegativeStrand);

    int readLength = bases.length;
    int[][][] ret = new int[numEvents][readLength][numCovariates];
    int idx = 0;
    for (int i = 0; i < bases.length; i++) {
      for (int j = 0; j < numCovariates; j++) {
        for (EventType event : EventType.values()) {
          ret[event.ordinal()][i][j] = keys[idx++];
        }
      }
    }
    return ret;
  }

  // This function is for unit testing only
  protected int[][][] computeCycleCovariates(final GATKRead read)
          throws AccelerationException {

    //int readLength = read.getReadBases().length;
    int readLength = read.getBases().length;
    //boolean isNegativeStrand = read.getReadNegativeStrandFlag();
    //boolean isReadPaired = read.getReadPairedFlag();
    boolean isNegativeStrand = read.isReverseStrand();
    boolean isReadPaired = read.isPaired();
    boolean isSecondOfPair;
    try {
      //isSecondOfPair = read.getSecondOfPairFlag();
      isSecondOfPair = read.isSecondOfPair();
    }
    catch (java.lang.IllegalStateException e) {
      isSecondOfPair = false;
    }
    //NGSPlatform ngsPlatform = read.getNGSPlatform();
    //final SAMReadGroupRecord rg = read.getReadGroup();
    final String rg = read.getReadGroup();
    //NGSPlatform ngsPlatform = NGSPlatform.fromReadGroupPL(rg.getPlatform());
    NGSPlatform ngsPlatform = NGSPlatform.fromReadGroupPL(rg);
    int platformType = ngsPlatform.getSequencerType() == SequencerFlowClass.DISCRETE ? 0 : 1;


    final int[] keys = computeCycleCovariatesNative(readLength, platformType,
            isNegativeStrand, isReadPaired, isSecondOfPair);

    int[][][] ret = new int[numEvents][readLength][numCovariates];
    int idx = 0;
    for (int i = 0; i < readLength; i++) {
      for (int j = 0; j < numCovariates; j++) {
        for (EventType event : EventType.values()) {
          ret[event.ordinal()][i][j] = keys[idx++];
        }
      }
    }
    return ret;
  }
  //protected int[][][] computeCycleCovariates(final GATKSAMRecord read)
  //protected int[][][] computeCycleCovariates(final SAMRecord read)
  protected int[][][] computeCycleCovariates(final GATKRead read, final SAMFileHeader header)
    throws AccelerationException {

    //int readLength = read.getReadBases().length;
    int readLength = read.getBases().length;
    //boolean isNegativeStrand = read.getReadNegativeStrandFlag();
    //boolean isReadPaired = read.getReadPairedFlag();
    boolean isNegativeStrand = read.isReverseStrand();
    boolean isReadPaired = read.isPaired();
    boolean isSecondOfPair;
    try {
      //isSecondOfPair = read.getSecondOfPairFlag();
      isSecondOfPair = read.isSecondOfPair();
    }
    catch (java.lang.IllegalStateException e) {
      isSecondOfPair = false;
    }
    //NGSPlatform ngsPlatform = read.getNGSPlatform();
    //final SAMReadGroupRecord rg = read.getReadGroup();
    //read.getReadGroup() is to get individual readID, like SEQ01,
    //header.getReadGroup(read.getReadGroup()))
    //is to use SEQ01 in header, to find the corresponding SAMReadGroupRecord,
    // after getting SAMReadGroupRecord, use .getSAMString() to get the whole line for that readGroup
    //@RG     ID:SEQ01        LB:L1   PL:illumina     SM:SEQ01
    // for a SAMReadGroupRecord, getReadGroupId() or getId() gives SEQ01,not whole line
    final String rg = header.getReadGroup(read.getReadGroup()).getSAMString();
    //NGSPlatform ngsPlatform = NGSPlatform.fromReadGroupPL(rg.getPlatform());
    NGSPlatform ngsPlatform = NGSPlatform.fromReadGroupPL(rg);
    int platformType = ngsPlatform.getSequencerType() == SequencerFlowClass.DISCRETE ? 0 : 1;


   // System.out.println("@@@@@@@@@");
   // System.out.println("!!!!!!!!!!!!!!!");
   // System.out.println(NGSPlatform.isKnown(rg)?"true":"false");
   // System.out.println(rg);
   // System.out.println(ngsPlatform.getDefaultPlatform());
   // System.out.println(platformType);
   // System.out.println("!!!!!!!");

    final int[] keys = computeCycleCovariatesNative(readLength, platformType,
            isNegativeStrand, isReadPaired, isSecondOfPair);

    int[][][] ret = new int[numEvents][readLength][numCovariates];
    int idx = 0;
    for (int i = 0; i < readLength; i++) {
      for (int j = 0; j < numCovariates; j++) {
        for (EventType event : EventType.values()) {
          ret[event.ordinal()][i][j] = keys[idx++];
        }
      }
    }
    return ret;
  }

  // This function is for unit testing only
  //protected int[] computeCovariates(final SAMRecord read) {
  protected int[] computeCovariates(final GATKRead read, final SAMFileHeader header) {

    //final byte[] bases = read.getReadBases();
    final byte[] bases = read.getBases();
    final byte[] baseQuals = read.getBaseQualities();

    //@@peipei

    //final byte[] baseInsertionQuals = ReadUtils.getBaseInsertionQualities(new SAMRecordToGATKReadAdapter(read));
    //final byte[] baseDeletionQuals = ReadUtils.getBaseDeletionQualities(new SAMRecordToGATKReadAdapter(read));
    final byte[] baseInsertionQuals = ReadUtils.getBaseInsertionQualities(read);
    final byte[] baseDeletionQuals = ReadUtils.getBaseDeletionQualities(read);


    //final boolean isNegativeStrand = read.getReadNegativeStrandFlag();
    //final boolean isReadPaired = read.getReadPairedFlag();
    //final boolean isSecondOfPair = read.getSecondOfPairFlag();
    final boolean isNegativeStrand = read.isReverseStrand();
    final boolean isReadPaired = read.isPaired();
    final boolean isSecondOfPair = read.isSecondOfPair();

    //final SAMReadGroupRecord rg = read.getReadGroup();
    //final NGSPlatform ngsPlatform = NGSPlatform.fromReadGroupPL(rg.getPlatform());
    //final String rg = read.getReadGroup();
    //final NGSPlatform ngsPlatform = NGSPlatform.fromReadGroupPL(rg);
    final String rg = header.getReadGroup(read.getReadGroup()).getSAMString();
    // equivalent to
    // final String rg = ReadUtils.getSAMReadGroupRecord(read,header).getSAMString();
    NGSPlatform ngsPlatform = NGSPlatform.fromReadGroupPL(rg);
    final int platformType = ngsPlatform.getSequencerType() == SequencerFlowClass.DISCRETE ? 0 : 1;

    String readGroupId;
    //if (FORCE_READGROUP != null) {
    //  readGroupId = FORCE_READGROUP;
    //}

    //final String platformUnit = rg.getPlatformUnit();
    //final String platformUnit = read.getAttributeAsString(SAMTag.PU.name());
    final String platformUnit = header.getReadGroup(read.getReadGroup()).getPlatformUnit();
    //readGroupId = platformUnit == null ? rg.getId() : platformUnit;
    readGroupId = platformUnit == null ? header.getReadGroup(read.getReadGroup()).getId() : platformUnit;

    final int[] keys = computeCovariatesNative(bases,
            baseQuals, baseInsertionQuals, baseDeletionQuals,
            readGroupId,
            isNegativeStrand, isReadPaired, isSecondOfPair,
            platformType);

    return keys;
  }

  //public byte[] calculateBAQArray(final SAMRecord read) {
  public byte[] calculateBAQArray(final GATKRead read, final ReferenceDataSource refDS) {
    //final GATKRead togatkread = new SAMRecordToGATKReadAdapter(read);
    if (baq.excludeReadFromBAQ(read)) {
      // in this case, simply return the original tag
      return BAQ.getBAQTag(read);
    }

    //int offset = baq.getBandWidth() / 2;
    //long readStart = includeClippedBases ? read.getUnclippedStart() : read.getAlignmentStart();
    //long start = Math.max(readStart - offset - ReadUtils.getFirstInsertionOffset(togatkread), 1);
    //long stop = (includeClippedBases ? read.getUnclippedEnd() : read.getAlignmentEnd()) + offset + ReadUtils.getLastInsertionOffset(togatkread);
    final SimpleInterval referenceWindow = BAQ.getReferenceWindowForRead(read, baq.getBandWidth());

    //if (stop > referenceReader.getSequenceDictionary().getSequence(read.getReferenceName()).getSequenceLength()) {
    if (referenceWindow.getEnd() > refDS.getSequenceDictionary().getSequence(read.getContig()).getSequenceLength()){
      return null;
    }
    else {
      //ReferenceSequence refSeq = referenceReader.getSubsequenceAt(read.getReferenceName(), start, stop);
      final ReferenceSequence refSeq = refDS.queryAndPrefetch(referenceWindow.getContig(), referenceWindow.getStart(), referenceWindow.getEnd());
      // preparing input arguments for native calls
      byte[] refBases = refSeq.getBases();
      //byte[] bases = read.getReadBases();
      byte[] bases = read.getBases();
      byte[] quals = read.getBaseQualities();      // in general we are overwriting quals, so just get a pointer to them

      // prepare cigar arrays
      final List<CigarElement> cigarElements = read.getCigar().getCigarElements();
      final int numCigarElements = cigarElements.size();
      byte[] cigarOps = new byte[numCigarElements];
      int[] cigarLens = new int[numCigarElements];
      int idx = 0;
      for (CigarElement elt : cigarElements) {
        cigarOps[idx] = CigarOperator.enumToCharacter(elt.getOperator());
        cigarLens[idx] = elt.getLength();
        idx++;
      }

      //System.out.println(Arrays.toString(refBases));
      //System.out.println(Arrays.toString(bases));
      //System.out.println(Arrays.toString(quals));
      //System.out.println(Arrays.toString(cigarOps));
      //System.out.println(Arrays.toString(cigarLens));


      byte[] bqTag = calculateBAQArrayNative(
          refBases, bases, quals,
          cigarOps, cigarLens,
          //(int)(start - readStart));
          (int)(referenceWindow.getStart() - read.getStart()));

      //byte[] bqTag2 = calculateBAQArrayNative_mock(
      //    refBases, bases, quals,
      //    cigarOps, cigarLens,
      //    (int)(start - readStart));

      if (bqTag == null) {
        final boolean readHasBAQTag = BAQ.hasBAQTag(read);
        //if (read.getStringAttribute(BAQ_TAG) != null) {
          // clear the attribute
          //read.setAttribute(BAQ_TAG, null);
        //}
        if(readHasBAQTag){
          read.clearAttribute(BAQ_TAG);
        }
        return null;
      }
      else {

        //boolean isError = false;
        //for (int i = 0; i < bqTag.length; i++) {
        //  if (bqTag[i] != bqTag2[i]) {
        //    logger.error("bqTag mismatch");
        //    isError = true;
        //    break;
        //  }
        //}
        //if (isError) {
        //  System.out.println("FALCON: " + Arrays.toString(bqTag));
        //  System.out.println("GATK: " + Arrays.toString(bqTag2));
        //}

        // TODO: is this really necessary?
        String baqStr = new String(bqTag);
        read.setAttribute(BAQ_TAG, baqStr);

        return bqTag;
      }
    }
  }

  //public double[][] calculateFractionalErrorArray(final SAMRecord read,
  //                    final SAMRecord org_read,
  //                    final ReferenceContext ref) {
  public double[][] calculateFractionalErrorArray(final GATKRead read,
                                                  final GATKRead org_read,
                                                  //final ReferenceContext ref) {
                                                  final ReferenceDataSource refDS) {

    // put the check of baq.isExcludeFromBAQ() outside, since we need
    // to pass the read BAQ tag to native if it is available
    //final GATKRead togatkread = new SAMRecordToGATKReadAdapter(read);
    //final boolean isExcludeFromBAQ = baq.excludeReadFromBAQ(togatkread);
    final boolean isExcludeFromBAQ = baq.excludeReadFromBAQ(read);
    //final byte[] readBAQArray = isExcludeFromBAQ ? BAQ.getBAQTag(togatkread) : null;
    final byte[] readBAQArray = isExcludeFromBAQ ? BAQ.getBAQTag(read) : null;

    // preparation for BAQ calculation
    int refOffset = -1;
    byte[] refForBAQ = null;

    // skip calcBAQFromHMM() and directly return readBAQ if it's not null
    if (!isExcludeFromBAQ) {
      //final int offset = baq.getBandWidth() / 2;
      //final long readStart = includeClippedBases ? read.getUnclippedStart() : read.getAlignmentStart();
      //final long start = Math.max(readStart - offset - ReadUtils.getFirstInsertionOffset(togatkread), 1);
      //final long stop = (includeClippedBases ? read.getUnclippedEnd() : read.getAlignmentEnd()) + offset + ReadUtils.getLastInsertionOffset(togatkread);
      final SimpleInterval referenceWindow = BAQ.getReferenceWindowForRead(read, baq.getBandWidth());

      //if (stop > referenceReader.getSequenceDictionary().getSequence(read.getReferenceName()).getSequenceLength()) {
      if (referenceWindow.getEnd() > refDS.getSequenceDictionary().getSequence(read.getContig()).getSequenceLength()){

        // meaning null return from calculateBAQArray
        // return null;
        // Note: Here should not return null, since if nErrors is zero BAQArray
        // calculation is skipped, native func knows from negative refOffset
        //logger.info("baq array cannot be calculated because stop is out-of-bound");
        ;
      }
      else {
        //refOffset = (int)(start - readStart);
        refOffset = (int)(referenceWindow.getStart() - read.getStart());

        // ref seq for baq calculation
        //ReferenceSequence refSeq = referenceReader.getSubsequenceAt(read.getReferenceName(), start, stop);
        final ReferenceSequence refSeq = refDS.queryAndPrefetch(referenceWindow.getContig(), referenceWindow.getStart(), referenceWindow.getEnd());
        refForBAQ = refSeq.getBases();
      }
    }

    // ref seqs for SNP error calculation
    //final byte[] refBases = Arrays.copyOfRange(ref.getBases(),
                        //read.getAlignmentStart() - org_read.getAlignmentStart(),
                        //ref.getBases().length + read.getAlignmentEnd() - org_read.getAlignmentEnd());
    //final byte[] refBases = Arrays.copyOfRange(refDS.getBases(),
    //                    read.getStart() - org_read.getStart(),
    //                    refDS.getBases().length + read.getEnd() - org_read.getEnd());

    final byte[] refBases = refDS.queryAndPrefetch(read.getContig(), read.getStart(), read.getEnd()).getBases();


    //final byte[] bases = read.getReadBases();
    final byte[] bases = read.getBases();
    final byte[] quals = read.getBaseQualities();      // in general we are overwriting quals, so just get a pointer to them
    final int readLength = bases.length;
    //final boolean isNegativeStrand = read.getReadNegativeStrandFlag();
    final boolean isNegativeStrand = read.isReverseStrand();

    // prepare cigar arrays
    final List<CigarElement> cigarElements = read.getCigar().getCigarElements();
    final int numCigarElements = cigarElements.size();
    byte[] cigarOps = new byte[numCigarElements];
    int[] cigarLens = new int[numCigarElements];
    int idx = 0;
    for (CigarElement elt : cigarElements) {
      cigarOps[idx] = CigarOperator.enumToCharacter(elt.getOperator());
      cigarLens[idx] = elt.getLength();
      idx++;
    }

    final double[] errors = calculateErrorsNative(
                                  bases, quals, refForBAQ, refBases,
                                  cigarOps, cigarLens,
                                  readBAQArray,
                                  isExcludeFromBAQ, isNegativeStrand,
                                  refOffset);

    if (errors == null) {
      // return null baqArray
      return null;
    }

    final double[][] ret = new double[numEvents][readLength];
    idx = 0;
    for (int k = 0; k < numEvents; k++) {
      for (int i = 0; i < readLength; i++) {
        ret[k][i] = errors[idx++];
      }
    }
    return ret;
  }

  // This function is used to update recalibration tables in
  // GATK BaseRecalibrator.map()
  //public int update(final SAMRecord read,
  //                  final SAMRecord org_read,
  //                  final ReferenceContext ref,
  //                  final boolean[] skips) throws AccelerationException
  public int update(final GATKRead read,
                    final GATKRead org_read,
                    final ReferenceDataSource refDS,
                    final SAMFileHeader header,
                    final boolean[] skips) throws AccelerationException
  {
    if (finalized) {
      throw new GATKException("Covariate table already finalized");
    }
    else if (!initialized) {
      throw new GATKException("Covariate table already finalized");
    }
    //final GATKRead togatkread = new SAMRecordToGATKReadAdapter(read);
    // put the check of baq.isExcludeFromBAQ() outside, since we need
    // to pass the read BAQ tag to native if it is available
    final boolean isExcludeFromBAQ = baq.excludeReadFromBAQ(read);
    final byte[] readBAQArray = isExcludeFromBAQ ? BAQ.getBAQTag(read) : null;

    // preparation for BAQ calculation
    int refOffset = 0;
    byte[] refForBAQ = null;

    if (!isExcludeFromBAQ) {
      //final int offset = baq.getBandWidth() / 2;
      //final long readStart = includeClippedBases ? read.getUnclippedStart() : read.getAlignmentStart();
      //final long start = Math.max(readStart - offset - ReadUtils.getFirstInsertionOffset(togatkread), 1);
      //final long stop = (includeClippedBases ? read.getUnclippedEnd() : read.getAlignmentEnd()) + offset + ReadUtils.getLastInsertionOffset(togatkread);
      final SimpleInterval referenceWindow = BAQ.getReferenceWindowForRead(read, baq.getBandWidth());
      if (referenceWindow.getEnd() > refDS.getSequenceDictionary().getSequence(read.getContig()).getSequenceLength()){
        int i;
      //if (stop > referenceReader.getSequenceDictionary().getSequence(read.getReferenceName()).getSequenceLength()) {
        // meaning null baqArray
        //return 0;
        // Note: Here should not return, since if nErrors is zero BAQArray
        // calculation is skipped, native func knows from negative refOffset
        //logger.info("baq array cannot be calculated because stop is out-of-bound");
        ;
      }
      else {
        //refOffset = (int)(start - readStart);
        refOffset = (int)(referenceWindow.getStart() - read.getStart());

        // ref seq for baq calculation
        //ReferenceSequence refSeq = referenceReader.getSubsequenceAt(read.getReferenceName(), start, stop);
        final ReferenceSequence refSeq = refDS.queryAndPrefetch(referenceWindow.getContig(), referenceWindow.getStart(), referenceWindow.getEnd());
        refForBAQ = refSeq.getBases();
      }
    }

    // ref seqs for SNP error calculation
    //final byte[] refBases = Arrays.copyOfRange(ref.getBases(),
    //                    read.getAlignmentStart() - org_read.getAlignmentStart(),
    //                    ref.getBases().length + read.getAlignmentEnd() - org_read.getAlignmentEnd());
    final byte[] refBases = refDS.queryAndPrefetch(read.getContig(), read.getStart(), read.getEnd()).getBases();

    // get inputs from read
    final byte[] bases = read.getBases();
    final byte[] baseQuals = read.getBaseQualities();
    final byte[] baseInsertionQuals = ReadUtils.getBaseInsertionQualities(read);
    final byte[] baseDeletionQuals = ReadUtils.getBaseDeletionQualities(read);

    final boolean isNegativeStrand = read.isReverseStrand();
    final boolean isReadPaired = read.isPaired();
    final boolean isSecondOfPair = read.isSecondOfPair();

    // prepare cigar arrays
    final List<CigarElement> cigarElements = read.getCigar().getCigarElements();
    final int numCigarElements = cigarElements.size();
    final byte[] cigarOps = new byte[numCigarElements];
    final int[] cigarLens = new int[numCigarElements];
    int idx = 0;
    for (CigarElement elt : cigarElements) {
      cigarOps[idx] = CigarOperator.enumToCharacter(elt.getOperator());
      cigarLens[idx] = elt.getLength();
      idx++;
    }

    //final SAMReadGroupRecord rg = read.getReadGroup();
    final String rg = header.getReadGroup(read.getReadGroup()).getSAMString();
    //final NGSPlatform ngsPlatform = NGSPlatform.fromReadGroupPL(rg.getPlatform());
    final NGSPlatform ngsPlatform = NGSPlatform.fromReadGroupPL(rg);

    //final NGSPlatform ngsPlatform = read.getNGSPlatform();
    final int platformType = ngsPlatform.getSequencerType() == SequencerFlowClass.DISCRETE ? 0 : 1;

    String readGroupId;
    //if (FORCE_READGROUP != null) {
    //  readGroupId = FORCE_READGROUP;
    //}
    //final GATKSAMReadGroupRecord rg = read.getReadGroup();
    //final String platformUnit = rg.getPlatformUnit();
    //readGroupId = platformUnit == null ? rg.getId() : platformUnit;
    final String platformUnit = header.getReadGroup(read.getReadGroup()).getPlatformUnit();
    readGroupId = platformUnit == null ? header.getReadGroup(read.getReadGroup()).getId() : platformUnit;


    // call native method to update RecalibrationTables
    return updateTableNative(refForBAQ, refBases,
        bases, baseQuals, baseInsertionQuals, baseDeletionQuals,
        cigarOps, cigarLens, readBAQArray,
        readGroupId, isNegativeStrand, isReadPaired, isSecondOfPair, isExcludeFromBAQ,
        platformType, refOffset,
        skips);
  }

  public byte[][] recalibrate(final GATKRead read, final SAMFileHeader header)
    throws AccelerationException
  {
    if (!initialized) return null;

    try {
    // these parameters are used to compute covariates
      //final GATKRead togatkread = new SAMRecordToGATKReadAdapter(read);
    //final byte[] bases = read.getReadBases();
      final byte[] bases = read.getBases();
    final byte[] baseQuals = read.getBaseQualities();
    final byte[] baseInsertionQuals = ReadUtils.getBaseInsertionQualities(read);
    final byte[] baseDeletionQuals = ReadUtils.getBaseDeletionQualities(read);

    //final boolean isNegativeStrand = read.getReadNegativeStrandFlag();
    //final boolean isReadPaired = read.getReadPairedFlag();
    //final boolean isSecondOfPair = read.getSecondOfPairFlag();
    final boolean isNegativeStrand = read.isReverseStrand();
    final boolean isReadPaired = read.isPaired();
    final boolean isSecondOfPair = read.isSecondOfPair();
    //final NGSPlatform ngsPlatform = read.getNGSPlatform();
    //final SAMReadGroupRecord rg = read.getReadGroup();
    //final NGSPlatform ngsPlatform = NGSPlatform.fromReadGroupPL(rg.getPlatform());
    final String rg = header.getReadGroup(read.getReadGroup()).getSAMString();
      //final NGSPlatform ngsPlatform = NGSPlatform.fromReadGroupPL(rg.getPlatform());
    final NGSPlatform ngsPlatform = NGSPlatform.fromReadGroupPL(rg);
    final int platformType = ngsPlatform.getSequencerType() == SequencerFlowClass.DISCRETE ? 0 : 1;

    String readGroupId;
    //if (FORCE_READGROUP != null) {
    //  readGroupId = FORCE_READGROUP;
    //}
    //final GATKSAMReadGroupRecord rg = read.getReadGroup();
    //final String platformUnit = rg.getPlatformUnit();
    //readGroupId = platformUnit == null ? rg.getId() : platformUnit;
    final String platformUnit = header.getReadGroup(read.getReadGroup()).getPlatformUnit();
    readGroupId = platformUnit == null ? header.getReadGroup(read.getReadGroup()).getId() : platformUnit;

    // return qualities of each events, if it's null it means
    // disableIndelQuals is true
    return recalibrateNative(bases,
            baseQuals, baseInsertionQuals, baseDeletionQuals,
            readGroupId,
            isNegativeStrand, isReadPaired, isSecondOfPair,
            platformType);
    } catch (Exception e) {
      throw new AccelerationException(e.getMessage());
    }
  }

  public StandardCovariateList getCovariates(){
    return covariates;
  }

  protected class RecalDatumTable {
    final public long[]   numOccurance;
    final public double[] numMismatches;
    final public int      tableSize;
    final public int[]    tableDimensions;
    final public int      tableDimensionsSize;

    public RecalDatumTable(int tableSize, int numDimensions) {
      this.numOccurance = new long[tableSize];
      this.numMismatches = new double[tableSize];
      this.tableSize = tableSize;
      this.tableDimensions = new int[numDimensions];
      this.tableDimensionsSize = numDimensions;
    }
  }

  public void updateRecalibrationTables() {
    if (initialized && !finalized) {
      // get data from native
      final RecalDatumTable[] nativeTables = getTableNative();

      // set table contents
      // skip readgroup table, generate it later, as the same as
      // RecalibrationEngine.finalizeData()
      for (int i = 0; i < numCovariates - 1; i++) {
        NestedIntegerArray<RecalDatum> table = recalTables.getTable(i+1);
        for (int k = 0; k < nativeTables[i].tableSize; k++) {
          if (nativeTables[i].numOccurance[k] == 0) continue;
          final long   numOccur = nativeTables[i].numOccurance[k];
          final double numError = nativeTables[i].numMismatches[k];

          int idx = k;
          final int[] keys = new int[nativeTables[i].tableDimensionsSize];

          // keys: event, rg, qual, cov
          // tableDimensions: numEvents, numRG, numQual, numCov
          for (int j = 0; j < nativeTables[i].tableDimensionsSize; j++) {
            keys[j] = idx % nativeTables[i].tableDimensions[j];
            idx = idx / nativeTables[i].tableDimensions[j];
          }
          int eventIndex = keys[0];
          byte qual = (byte)keys[2];

          if (i == 0) { // Qual table
            table.put(new RecalDatum(numOccur, numError, qual),
                      keys[1], keys[2], eventIndex);
          }
          else { // optional tables
            table.put(new RecalDatum(numOccur, numError, qual),
                      keys[1], keys[2], keys[3], eventIndex);
          }
        }
      }
    }
  }

  public RecalDatumTable[] getTables() {
    return getTableNative();
  }
  public RecalibrationTables getDebugTable(){
      return recalTables;
  }


  public RecalibrationTables getRecalibrationTables() {
    if (!finalized) {
      updateRecalibrationTables();
    }
    return recalTables;
  }

  public RecalibrationTables getFinalRecalibrationTables() {

    finalizeData();


    return recalTables;
  }

  public void updateReadGroupCovariates() {
    if (initialized) {
      // update ReadGroupCovariate class to reflect the key-idx mapping
      if (!(covariates.get(0) instanceof ReadGroupCovariate)) {
        throw new GATKException("The first requested covariate should be " +
        "ReadGroupCovariate.");
      }
      updateReadGroupCovariatesNative((ReadGroupCovariate)covariates.get(0));
    }
  }

  public void finalizeData() {
    if (!initialized || finalized) return;


    updateReadGroupCovariates();

    // get latest recal table
    updateRecalibrationTables();


    // here, qualityscore table, context and cycle table are updated
    // then update ReadGroupTable table
    // finalize RecalibrationTables
    // renaming for GATK
    RecalibrationTables finalRecalibrationTables = recalTables;
    // start GATK code from org.broadinstitute.gatk.tools.walkers.bqsr.RecalibrationEngine
    final NestedIntegerArray<RecalDatum> byReadGroupTable = finalRecalibrationTables.getReadGroupTable();
    final NestedIntegerArray<RecalDatum> byQualTable = finalRecalibrationTables.getQualityScoreTable();
    // iterate over all values in the qual table
    int counter=0;
    for ( final NestedIntegerArray.Leaf<RecalDatum> leaf : byQualTable.getAllLeaves() ) {

        final int rgKey = leaf.keys[0];
        final int eventIndex = leaf.keys[2];
        final RecalDatum rgDatum = byReadGroupTable.get(rgKey, eventIndex);
        final RecalDatum qualDatum = leaf.value;
        if ( rgDatum == null ) {
            // create a copy of qualDatum, and initialize byReadGroup table with it
            byReadGroupTable.put(new RecalDatum(qualDatum), rgKey, eventIndex);
            //System.out.printf("null branch  rgKey: %d, eventIndex: %d , qualDatum: %s\n", rgKey, eventIndex, qualDatum.toString());
        } else {
            // combine the qual datum with the existing datum in the byReadGroup table
            //System.out.printf("else branch  rgKey: %d, eventIndex: %d, qualDatum: %s\n", rgKey, eventIndex, qualDatum.toString());
            rgDatum.combine(qualDatum);
        }
        counter+=1;
    }
    System.out.printf("@@@ QualityTable leaves counter is : %d\n",counter);

    logger.debug("Free resource in native space");
    finalizeNative();

    // new in gatk4, comments from gatk4 are
    /* To replicate the results of BQSR whether or not we save tables to disk (which we need in Spark),
     * we need to trim the numbers to a few decimal placed (that's what writing and reading does).
     */
    BaseRecalibrationEngine.roundTableValues(recalTables);
    // end of GATK code
    finalized = true;
  }

  // helper routines for BAQ calculation
  private static boolean stateIsIndel(int state) {
    return (state & 3) != 0;
  }

  /** decode the bit encoded state array values */
  private static int stateAlignedPosition(int state) {
    return state >> 2;
  }

  private static byte capBaseByBAQ(final byte oq,
                                   final byte bq,
                                   final int state,
                                   final int expectedPos)
  {
    byte b;
    boolean isIndel = stateIsIndel(state);
    int pos = stateAlignedPosition(state);
    if (isIndel || pos != expectedPos) // we are an indel or we don't align to our best current position
      b = minBaseQual; // just take b = minBaseQuality
    else
      b = bq < oq ? bq : oq;

    return b;
  }

  // same function as calculateBAQArrayNative()
  public byte[] calculateBAQArrayNative_mock(
              byte[] refBases,
              byte[] bases,
              byte[] basesQuals,
              byte[] cigarOps,
              int[] cigarLens,
              int refOffset)
  {
    int readLength = bases.length;

    // calculate query range (BAQ.calculateQueryRange())
    int queryStart = -1;
    int queryStop = -1;
    int readI = 0;
    // iterate over the cigar elements to determine the start and stop of the read bases for the BAQ calculation
    //for (CigarElement elt : cigarElements) {
    for (int i = 0; i < cigarOps.length; i++) {
      switch (cigarOps[i]) {
        case 'N':
          return null; // cannot handle these
        case 'H':
        case 'P':
        case 'D':
          break; // ignore pads, hard clips, and deletions
        case 'I':
        case 'S':
        case 'M':
        case '=': //case EQ:
        case 'X':
          int prev = readI;
          readI += cigarLens[i];
          if ( includeClippedBases || cigarOps[i] != 'S') {
            if ( queryStart == -1 )
              queryStart = prev;
            queryStop = readI;
          }
          // in the else case we aren't including soft clipped bases, so we don't update
          // queryStart or queryStop
          break;
        default: throw new GATKException("BUG: Unexpected CIGAR element in read");
      }
    }

    if (queryStop == queryStart) {
      // this read is completely clipped away, and yet is present in the file for some reason
      // usually they are flagged as non-PF, but it's possible to push them through the BAM
      //System.err.printf("WARNING -- read is completely clipped away: " + read.format());
      return null;
    }

    //BAQ.BAQCalculationResult baqResult = new BAQ.BAQCalculationResult(query, quals, ref);
    int queryLen = queryStop - queryStart;
    final int[] baqState = new int[basesQuals.length];
    final byte[] baqBq = new byte[basesQuals.length];

    Arrays.fill(baqState, 0);
    Arrays.fill(baqBq, (byte)0);

    baq.hmm_glocal(refBases, bases, queryStart, queryLen, basesQuals, baqState, baqBq);

    readI = 0;
    int refI = 0;
    //for (CigarElement elt : cigarElements) {
    for (int k = 0; k < cigarOps.length; k++) {
      int l = cigarLens[k];
      switch (cigarOps[k]) {
        case 'N': // cannot handle these
          return null;
        case 'H':
        case 'P': // ignore pads and hard clips
          break;
        case 'S':
          refI += l; // move the reference too, in addition to I
        case 'I':
          // todo -- is it really the case that we want to treat I and S the same?
          for (int i = readI; i < readI + l; i++) baqBq[i] = basesQuals[i];
          readI += l;
          break;
        case 'D':
          refI += l;
          break;
        case 'M':
        case '=':
        case 'X':
          for (int i = readI; i < readI + l; i++) {
            int expectedPos = refI - refOffset + (i - readI);
            baqBq[i] = capBaseByBAQ(basesQuals[i], baqBq[i], baqState[i], expectedPos);
          }
          readI += l; refI += l;
          break;
        default:
          throw new GATKException("BUG: Unexpected CIGAR element in read");
      }
    }
    if (readI != readLength) // odd cigar string
      System.arraycopy(basesQuals, 0, baqBq, 0, baqBq.length);

    //BAQ.encodeBQTag(read, hmmResult.bq);

    final byte[] bqTag = new byte[baqBq.length];
    for (int i = 0; i < bqTag.length; i++) {
      final int bq = (int)basesQuals[i] + 64;
      final int baq_i = (int)baqBq[i];
      final int tag = bq - baq_i;
      // problem with the calculation of the correction factor; this is our problem
      if (tag < 0) {
        throw new GATKException("BAQ tag calculation error. BAQ value above base quality");
      }
      // the original quality is too high, almost certainly due to using the wrong encoding in the BAM file
      if (tag > Byte.MAX_VALUE) {
        throw new GATKException("we encountered an extremely high quality score (" + (int)basesQuals[i] + ") with BAQ correction factor of " + baq_i);
      }
      bqTag[i] = (byte)tag;
    }
    return bqTag;
  }

  private native static void bqsrInitReadGroupNative(
          String readGroupId
  );

  private native static void initNative(
      int numReadGroups,
      int numEvents,
      int numCovariates,
      int[] covariatesDimensions,
      byte LOW_QUAL_TAIL,
      int MISMATCHES_CONTEXT_SIZE,
      int INDELS_CONTEXT_SIZE,
      int MAXIMUM_CYCLE_VALUE,
      int CUSHION_FOR_INDEL);

  private native static void initNative(
      int numEvents,
      Covariate[] covariates,
      //StandardCovariateList covariates,
      byte[] quantizationTable,
      byte[] staticQuantizedMapping,
      boolean disableIndelQuals,
      int preserveQLessThan,
      double globalQScorePrior,
      boolean emitOriginalQuals,
      byte LOW_QUAL_TAIL,
      int MISMATCHES_CONTEXT_SIZE,
      int INDELS_CONTEXT_SIZE,
      int MAXIMUM_CYCLE_VALUE,
      int CUSHION_FOR_INDEL);

  private native void setTableNative(
      long numOccurance,
      double numMismatches,
      double estimatedQReported,
      int[] keys,
      int  cov_idx);

  // These routines are used for unit tests, in real use
  // all of them will be combined into updateTableNative()
  // - computeContextCovariates()
  // - computeCycleCovariates()
  // - computeCovariates()
  // - calculateBAQArray()
  // - calculateErrors()
  private native int[] computeContextCovariatesNative(
      byte[] bases,
      byte[] quals,
      boolean isNegativeStrand);

  private native int[] computeCycleCovariatesNative(
      int readLength,
      int platformType,
      boolean isNegativeStrand,
      boolean isReadPaired,
      boolean isSecondOfPair);

  private native int[] computeCovariatesNative(
      byte[] bases,
      byte[] baseQuals,
      byte[] baseInsertionQuals,
      byte[] baseDeletionQuals,
      String readGroupId,
      boolean isNegativeStrand,
      boolean isReadPaired,
      boolean isSecondOfPair,
      int platformType);

  private native byte[] calculateBAQArrayNative(
      byte[] refBases,
      byte[] bases,
      byte[] baseQuals,
      byte[] cigarOps,
      int[] cigarLens,
      int refOffset);

  private native double[] calculateErrorsNative(
      byte[] bases,
      byte[] quals,
      byte[] refForBAQ,
      byte[] refBases,
      byte[] cigarOps,
      int[]  cigarLens,
      byte[] readBAQArray,
      boolean isExcludeFromBAQ,
      boolean isNegativeStrand,
      int refOffset);

  // This is the actual native impl for applications
  private native int updateTableNative(
      byte[] refForBAQ,
      byte[] refBases,
      byte[] bases,
      byte[] baseQuals,
      byte[] baseInsertionQuals,
      byte[] baseDeletionQuals,
      byte[] cigarOps,
      int[]  cigarLens,
      byte[] readBAQArray,
      String readGroupId,
      boolean isNegativeStrand,
      boolean isReadPaired,
      boolean isSecondOfPair,
      boolean isExcludeFromBAQ,
      int platformType,
      int refOffset,
      boolean[] skips);

  private native RecalDatumTable[] getTableNative();

  private native void updateReadGroupCovariatesNative(ReadGroupCovariate cov);

  private native void putRecalDatumNative();

  private native byte[][] recalibrateNative(
      byte[] bases,
      byte[] baseQuals,
      byte[] baseInsertionQuals,
      byte[] baseDeletionQuals,
      String readGroupId,
      boolean isNegativeStrand,
      boolean isReadPaired,
      boolean isSecondOfPair,
      int platformType);

  // for unit test
  public native double log10QempPriorNative(double qemp, double qreported);

  public native double log10QempLikelihoodNative(
      double qemp,
      long numObservations,
      long numMismatches);

  public native double bayesianEstimateOfEmpiricalQualityNative(
      long numObservations,
      long numMismatches,
      double qreported);

  private native void finalizeNative();
}
