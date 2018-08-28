package org.broadinstitute.hellbender.utils.smithwaterman;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.SkipException;
import java.util.*;


public class SmithWatermanFalconAlignerUnitTest extends SmithWatermanAlignerAbstractUnitTest {

    /*
    *Test the Intel Aligner Native Interface
    */


    @Override
    protected SmithWatermanFalconAligner getAligner() {

        boolean loaded = true;
        SmithWatermanFalconAligner aligner = null;
        try {
            aligner = new SmithWatermanFalconAligner();
        }
        catch (UserException.HardwareFeatureException e ) {
            loaded = false;
        }
        if(!loaded) {
            //Skip test if correct version of AVX is not supported
            throw new SkipException("Falcon SmithWaterman is not supported on this system or the library is not available");
        }
        return aligner;
    }

}
