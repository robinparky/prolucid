/**
 * @file Ms1PeakFinder.java
 * This is the source file for edu.scripps.pms.util.spectrum.Ms1PeakFinder
 * @author Tao Xu
 * @date $Date: 2008/10/24 22:39:14 $
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

// this is the class for Christian's lock mass checking 
public class Ms1PeakFinder {
    public static final String USAGE = "java Ms1PeakFinder m2z tolerance_in_ppm  minRetentionTime maxRetentionTime ms1filename";
    public static final int MILLION = 1000000;
    public static double mz = 0;
//    private double [] xcorrs = new double[100000];
//    private double [] deltaMasses = new double[100000];
    
    private static void processMs1File(String ms1FileName, double minM2z, double maxM2z,  
                double minRT, double maxRT) throws IOException {
        SpectrumReader sr = new SpectrumReader(ms1FileName, "ms1");
        int numSpectra = 0;

        for(Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList pl = it.next();
            int scanno = pl.getLoscan();
           
            double retentionTime = pl.getRetentionTime();
            double injectiontime = pl.getIonInjectionTime();
            if(pl.numPeaks() > 0 && retentionTime >= minRT) {
                if(retentionTime > maxRT) {
                    break;
                }
                double deltamassinppm = 100000; 
                double intensity =pl.getMaxIntensity(minM2z, maxM2z); 
                Peak p = pl.getMaxIntensePeak(minM2z, maxM2z);
                if(p != null) {
                    deltamassinppm = ((p.getM2z()-mz)/mz)*MILLION;
                    System.out.println(scanno + "\t" + retentionTime + "\t" + injectiontime + "\t" + mz + "\t" + p.getM2z() + "\t" + deltamassinppm + "\t" + p.getIntensity()+ "\t" + p.getChargeState());
                } else {
   
                    System.out.println(scanno +"\t" + retentionTime + "\t" + injectiontime + "\t" +  mz + "\t0\t" + deltamassinppm + "\t0\t0");
                }
            }
        }
    }

    public static void main(String args[]) throws Exception {
        try {
           mz = Double.parseDouble(args[0]); 
           double toleranceinppm = Double.parseDouble(args[1]); 
            
           double minRetentionTime = Double.parseDouble(args[2]); 
           double maxRetentionTime = Double.parseDouble(args[3]); 
           String ms1file = args[4];
           double tolerance = mz*toleranceinppm/MILLION;
           double minM2z = mz - tolerance;
           double maxM2z = mz + tolerance;
           System.out.println("Find the " + mz + " peak in " +  ms1file + " with " +  toleranceinppm + " ppmtolerance and " + tolerance  + " amu tolerance");
           System.out.println("scan_number\tretentiontime\tioninjectiontime\ttargetmz\tmz\tdeltamassinpps\tintensity\tnumpeaksinrange");
           processMs1File(ms1file, minM2z, maxM2z, minRetentionTime, maxRetentionTime);

        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
        }
    }
}
