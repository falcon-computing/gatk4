package com.falconcomputing.genomics.bqsr;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMFileHeader;

//import org.broadinstitute.hellbender.engine.datasources.reference.*;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
//import htsjdk.samtools.reference.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.*;
import org.broadinstitute.hellbender.tools.walkers.bqsr.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceContext.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.baq.BAQ;
//import org.broadinstitute.hellbender.utils.sam.*;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;

import java.io.File;
import java.nio.file.*;
import java.util.*;
import java.net.*;

//
import org.broadinstitute.hellbender.utils.recalibration.BaseRecalibrationEngine;


public class TestHelper extends BaseRecalibrationEngine {

  final static int BUFFER = 10000;
  final static byte[] Xs = new byte[BUFFER];

  private BAQ baq; // BAQ the reads on the fly to generate the alignment uncertainty vector

  private ReferenceDataSource refDataSource;
  //private ReferenceSequenceFile refReader;
  private GenomeLocParser glocParser;

  public TestHelper(Path referenceFile, final RecalibrationArgumentCollection recalArgs ,SAMFileHeader header) {
    super(recalArgs, header);

    //refDataSource = new ReferenceDataSource(referenceFile);
    refDataSource = ReferenceDataSource.of(referenceFile);
    //refReader = refDataSource.getReference();
    //glocParser = new GenomeLocParser(refReader);
    double BAQGOP = BAQ.DEFAULT_GOP;
    baq = new BAQ(BAQGOP); // setup the BAQ object with the provided gap open penalty
  }

  //public ReferenceSequenceFile getRefReader() {
  //  //return refDataSource.getReference();
  //  return refReader;
  //}

  public ReferenceDataSource getRefDataSource () {
    return refDataSource;
  }

/*
  //public int[] falconCalculateIsSNP(final GATKSAMRecord read, final ReferenceContext ref, final GATKSAMRecord org_read) {
  public int[] falconCalculateIsSNP(final SAMRecord read, final ReferenceContext ref, final SAMRecord org_read) {
    return calculateIsSNP(read, ref, org_read);
  }

  //public int[] falconCalculateIsIndel(final GATKSAMRecord read, EventType event) {
  public int[] falconCalculateIsIndel(final SAMRecord read, EventType event) {
    return calculateIsIndel(read, event);
  }

  public int falconNumEvents(int[] isSNP, int[] isInsertion, int[] isDeletion) {
    return nEvents(isSNP, isInsertion, isDeletion);
  }
*/
  public double[] falconCalculateFractionalErrorArray(final int[] errorArray, final byte[] baqArray) {
    return calculateFractionalErrorArray(errorArray, baqArray);
  }

  //public byte[] falconFlatBAQArray(final GATKSAMRecord read) {
  public byte[] falconFlatBAQArray(final GATKRead read) {
    return flatBAQArray(read);
  }

  //public byte[] falconCalculateBAQArray(final GATKSAMRecord read) {
  public byte[] falconCalculateBAQArray(final GATKRead read) {
      //baq.baqRead(read, refReader, BAQ.CalculationMode.RECALCULATE, BAQ.QualityMode.ADD_TAG);
      baq.baqRead(read, refDataSource, BAQ.CalculationMode.RECALCULATE, BAQ.QualityMode.ADD_TAG);
      return BAQ.getBAQTag(read);
  }


/*
  private byte[] getReferenceBases(GenomeLoc genomeLoc) {
    SAMSequenceRecord sequenceInfo = refReader.getSequenceDictionary().getSequence(genomeLoc.getContig());

    long start = genomeLoc.getStart();
    long stop = Math.min( genomeLoc.getStop(), sequenceInfo.getSequenceLength() );

    // Read with no aligned bases?  Return an empty array.
    if(stop - start + 1 == 0)
    return new byte[0];

    ReferenceSequence subsequence = refReader.getSubsequenceAt(genomeLoc.getContig(), start, stop);

    int overhang = (int)(genomeLoc.getStop() - stop);
    if ( overhang > 0 ) {
      if ( overhang > BUFFER ) // todo -- this is a bit dangerous
      throw new GATKException("Insufficient buffer size for Xs overhanging genome -- expand BUFFER");
      byte[] all = new byte[subsequence.getBases().length + overhang];
      System.arraycopy(subsequence.getBases(), 0, all, 0, subsequence.getBases().length);
      System.arraycopy(Xs, 0, all, subsequence.getBases().length, overhang);
      return all;
    } else {
      // fast path
      return subsequence.getBases();
    }
  }
*/
/*
  private class Provider implements ReferenceContext.ReferenceContextRefProvider {
    GenomeLoc loc;
    public Provider(GenomeLoc loc) {
      this.loc = loc;
    }
    public byte[] getBases() {
      return getReferenceBases(loc);
    }
  }

  //public ReferenceContext getRefContext(GATKSAMRecord read) {
  public ReferenceContext getRefContext(SAMRecord read) {
    final GenomeLoc loc = glocParser.createGenomeLoc(read);
    return new ReferenceContext(glocParser, loc, loc, new Provider(loc));
  }

*/
}

