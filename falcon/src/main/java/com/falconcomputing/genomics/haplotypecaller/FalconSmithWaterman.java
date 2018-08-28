package com.falconcomputing.genomics.haplotypecaller;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import com.falconcomputing.genomics.NativeLibraryLoader;
import com.falconcomputing.genomics.AccelerationException;

import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignerNativeBinding;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWNativeAlignerResult;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;

import java.io.File;
import java.lang.Object;
import java.nio.charset.Charset;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

/**
 * Provides a native SmithWaterman implementation accelerated by Falcon.
 */

public class FalconSmithWaterman implements SWAlignerNativeBinding {

    private final static Logger logger = LogManager.getLogger(FalconSmithWaterman.class);
    private static final String NATIVE_LIBRARY_NAME = "falcon_genomics";
    private String nativeLibraryName = "falcon_genomics";
    private static boolean initialized = false;
    private static boolean loaded = false;
    // private IntelGKLUtils gklUtils = new IntelGKLUtils();
    //
    private enum State {
        MATCH,
        INSERTION,
        DELETION,
        CLIP
    }

    private static CigarElement makeElement(final State state, final int length) {
        CigarOperator op = null;
        switch (state) {
            case MATCH: op = CigarOperator.M; break;
            case INSERTION: op = CigarOperator.I; break;
            case DELETION: op = CigarOperator.D; break;
            case CLIP: op = CigarOperator.S; break;
        }
        return new CigarElement(length, op);
    }

    void setNativeLibraryName(String nativeLibraryName) {
        this.nativeLibraryName = nativeLibraryName;
    }

    public FalconSmithWaterman() {

        setNativeLibraryName(NATIVE_LIBRARY_NAME);
    }


    /**
    * Loads the native library, if it is supported on this platform. <p>
    * Returns false if AVX is not supported. <br>
    * Returns false if the native library cannot be loaded for any reason. <br>
    *
    * @param tempDir  directory where the native library is extracted or null to use the system temp directoryt
    * @return  true if the native library is supported and loaded, false otherwise
    */

    @Override
    public synchronized boolean load(File tempDir) {
        if (!loaded) {
            if (!NativeLibraryLoader.load(tempDir, NATIVE_LIBRARY_NAME)) {
                return false;
            }
            loaded = true;
        }

        /*
         Initializes the function pointers to use machine specific optimized code
         */
        init_native();
        return true;
    }


    /**
     *Implements the native implementation of SmithWaterman, and returns the Cigar String and alignment_offset
     *
     * @param refArray array of reference data
     * @param altArray array of alternate data
     *
     */

    @Override
    public SWNativeAlignerResult align(byte[] refArray, byte[] altArray, SWParameters parameters, SWOverhangStrategy overhangStrategy)
    {
        if ( refArray == null || refArray.length == 0 || altArray == null || altArray.length == 0 )
            throw new IllegalArgumentException("Non-null, non-empty sequences are required for the Smith-Waterman calculation");
        
        int intStrategy = OverhangToInt(overhangStrategy);

        int[] state_list;
        int[] length_list;
        int cigar_list_size;
        int alignment_offset = 0;
        int ret_code;
        
        FalconSWResultsNative alignment_results_native = align_native(refArray, altArray, parameters.getMatchValue(), parameters.getMismatchPenalty(), parameters.getGapOpenPenalty(), parameters.getGapExtendPenalty(), intStrategy); /* remove cutoff*/

        state_list = alignment_results_native.state_list;
        length_list = alignment_results_native.length_list;
        cigar_list_size = alignment_results_native.cigar_list_size;
        alignment_offset = alignment_results_native.alignment_offset;
        ret_code = alignment_results_native.ret_code;

        if (ret_code < 0) {
            logger.debug("error happens during native call");
        }

        final List<CigarElement> lce = new ArrayList<CigarElement>(5);
        for(int i = 0; i < cigar_list_size; ++i) {
            lce.add(makeElement(IntToState(state_list[i]), length_list[i]));
        }

        return (new SWNativeAlignerResult(TextCigarCodec.encode(AlignmentUtils.consolidateCigar(new Cigar(lce))), alignment_offset));
    }

    /*
    public SmithWatermanAlignment align_top(byte[] refArray, byte[] altArray, SWParameters parameters, SWOverhangStrategy overhangStrategy)
    {
        if ( refArray == null || refArray.length == 0 || altArray == null || altArray.length == 0 )
            throw new IllegalArgumentException("Non-null, non-empty sequences are required for the Smith-Waterman calculation");
        
        int intStrategy = OverhangToInt(overhangStrategy);

        int[] state_list;
        int[] length_list;
        int cigar_list_size;
        int alignment_offset = 0;
        int ret_code;
        
        FalconSWResultsNative alignment_results_native = align_native(refArray, altArray, parameters.getMatchValue(), parameters.getMismatchPenalty(), parameters.getGapOpenPenalty(), parameters.getGapExtendPenalty(), intStrategy); 

        state_list = alignment_results_native.state_list;
        length_list = alignment_results_native.length_list;
        cigar_list_size = alignment_results_native.cigar_list_size;
        alignment_offset = alignment_results_native.alignment_offset;
        ret_code = alignment_results_native.ret_code;

        if (ret_code < 0) {
            logger.debug("error happens during native call");
        }

        final List<CigarElement> lce = new ArrayList<CigarElement>(5);
        for(int i = 0; i < cigar_list_size; ++i) {
            lce.add(makeElement(IntToState(state_list[i]), length_list[i]));
        }
        return (new SmithWatermanAlignment(AlignmentUtils.consolidateCigar(new Cigar(lce)), alignment_offset));
    }
    */

    public State IntToState(int x) {
        switch(x) {
        case 0: 
            return State.MATCH;
        case 1:
            return State.INSERTION;
        case 2:
            return State.DELETION;
        case 4:
            return State.CLIP;
        }
        return null;
    }

    public int OverhangToInt(SWOverhangStrategy x) {
        switch(x) {
        case SOFTCLIP:
            return 0;
        case INDEL:
            return 1;
        case LEADING_INDEL:
            return 2;
        case IGNORE:
            return 3;
        }
        return -1;
    }

    public void close()
    {
        done_native();
    }

    public class FalconSWResultsNative{
        int[] state_list;
        int[] length_list;
        int cigar_list_size;
        int alignment_offset;
        int ret_code;
        public FalconSWResultsNative(int[] state_list, int[] length_list, int cigar_list_size, int alignment_offset, int ret_code) {
            this.state_list = state_list;
            this.length_list = length_list;
            this.cigar_list_size = cigar_list_size;
            this.alignment_offset = alignment_offset;
            this.ret_code = ret_code;
        }
    }

    private native static void init_native();
    private native FalconSWResultsNative align_native(final byte[] reference, final byte[] alternate, final int w_match, final int w_mismatch, final int w_open, final int w_extend, final int overhang_strategy_native);
    private native static void done_native();
}
