/**
 * @file SpectrumFilter.java
 * This is the source file for edu.scripps.pms.util.spectrum.SpectrumFilter
 * @author Tao Xu
 * @date $Date: 2006/06/25 07:37:17 $
 */



import edu.scripps.pms.util.io.*;

//import edu.scripps.pms.util.io.DTASelectFilterReader; 
import edu.scripps.pms.util.dtaselect.Peptide; 
import edu.scripps.pms.util.dtaselect.Protein; 
import java.io.*;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import edu.scripps.pms.util.spectrum.*;

// this is the program to get ms2 spectrum that contain certain signature peak, 225 and 196 for Shascha.
// The spectrum has to have at least 20 peaks and the mass has to be >= 800   
public class SpectrumFilter {
    public static final String USAGE = "\n!!! USAGE: spectrumFilter ms2File massTolerance !!!\n";
    private static final double SIGNATUREMASS1 = 225;
    private static final double SIGNATUREMASS2 = 196;
    private static double maxSign1;
    private static double minSign1;
    private static double maxSign2;
    private static double minSign2;
    private static double maxSign;
    private static double minSign;
    private static final double MINMASS = 800;
    private static final int MINNUMPEAKS = 20;
//    private double [] xcorrs = new double[100000];
//    private double [] deltaMasses = new double[100000];
    
    private static void processMs1File(String ms2FileName) throws IOException {
        SpectrumReader sr = new SpectrumReader(ms2FileName, "ms2");
        boolean isFirstSpectrum = true;
        for(Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList pl = it.next();
            if(hasSignaturePeak(pl) && pl.numPeaks() >= MINNUMPEAKS && pl.getMaxPrecursorMass() >= MINMASS) {
                if(isFirstSpectrum) {
                    System.out.print(pl.getSpectrumWithHlines());
                    isFirstSpectrum = false;
                } else {
                    System.out.print(pl.getSpectrumWithoutHlines());
                }
            }
        }
    }

    public static boolean hasSignaturePeak(PeakList pl) {
        boolean hasSignaturePeak = false;
        for(Iterator<Peak> it = pl.getPeaks(); it.hasNext();) {
            Peak p = it.next();  
            double m2z = p.getM2z();
            if(m2z > maxSign) {
                break;
            }
            if(m2z >= minSign1 && m2z <= maxSign1) {
                return true;
            } 
            if(m2z >= minSign2 && m2z <= maxSign2) {
                return true;
            } 
        } 
        return hasSignaturePeak; 
    }
    public static void main(String args[]) throws Exception {
        try {
           double tolerance = Double.parseDouble(args[1]); 
            
           maxSign1 = SIGNATUREMASS1 + tolerance;
           minSign1 = SIGNATUREMASS1 - tolerance;
           maxSign2 = SIGNATUREMASS2 + tolerance;
           minSign2 = SIGNATUREMASS2 - tolerance;
           maxSign = maxSign1 > maxSign2? maxSign1 : maxSign2;
           minSign = minSign1 < minSign2? minSign1 : minSign2;
//System.out.println("maxSignature: " + maxSign + "\tminSignature: " + minSign);
           processMs1File(args[0]);
        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.out.println(USAGE);
            System.exit(1);
        }
       // System.out.println("Finished");
    }
}
