/**
 * @file CommonPeakFinderNew.java
 * This is the source file for edu.scripps.pms.util.spectrum.CommonPeakFinderNew
 * @author Tao Xu
 * @date $Date
 */

import java.util.ArrayList;
import java.util.Iterator;

public class ExtendedPeak {
    double xtractedmass;
    ArrayList<XtractedPeak> peaks = new ArrayList<XtractedPeak>(15);
    
    ExtendedPeak(double mass) {
        xtractedmass = mass;
    }
    public void addPeak(XtractedPeak p) {
        peaks.add(p);
    }
}
