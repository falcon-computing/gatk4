package org.broadinstitute.hellbender.utils.smithwaterman;


import com.falconcomputing.genomics.haplotypecaller.FalconSmithWaterman;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignerNativeBinding;
import org.broadinstitute.hellbender.exceptions.UserException;
/**
 * SmithWatermanFalconAligner class that converts instance of {@link SWAlignerNativeBinding} into a {@link SmithWatermanFalconAligner}
 * This is optimized for Intel Architectures and can fail if Machine does not support AVX and will throw {@link UserException}
 */

public final class SmithWatermanFalconAligner implements SmithWatermanAligner {



    private final SWAlignerNativeBinding smithWaterman = new FalconSmithWaterman();

    /*
    * Generate SWAlignerWrapper instance
    */
    private final SWNativeAlignerWrapper alignerWrapper = new SWNativeAlignerWrapper(smithWaterman);


    /**
     * Create a new SW pairwise aligner, which is implementation of smith waterman aligner
     */

    public SmithWatermanFalconAligner() throws UserException.HardwareFeatureException {
        final boolean isSupported = smithWaterman.load(null);
        if (!isSupported) {
            throw new UserException.HardwareFeatureException("Machine does not support Falcon SmithWaterman.");
        }
    }

    /**
     * Aligns the alternate sequence to the reference sequence
     *
     * @param reference  ref sequence
     * @param alternate  alt sequence
     */
    @Override
    public SmithWatermanAlignment align(final byte[] reference, final byte[] alternate, final SWParameters parameters, final SWOverhangStrategy overhangStrategy) {
        return alignerWrapper.align(reference, alternate, parameters, overhangStrategy);
    }


    /**
     * Close the aligner
     */
    @Override
    public void close() {
        alignerWrapper.close();
    }
}
