/**
 * @file PeakIntensityMean.java
 * This is the source file for edu.scripps.pms.util.spectrum.PeakIntensityMean
 * @author Tao Xu
 * @date $Date: 2014/08/12 19:29:19 $
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

// this is the class for Aaron's signal to noise computation
public class PeakIntensityMean {
    public static final String USAGE = "java PeakIntensityMean m2z tolerance percentageAsNoise minRetentionTime maxRetentionTime Spectrum_File_Type(ms1 or ms2)";
    public static final int MILLION = 1000000;
//    private double [] xcorrs = new double[100000];
//    private double [] deltaMasses = new double[100000];
    
    public static List<String> getFiles(String dir, String suffix) {
        
         ArrayList<String> files = new ArrayList<String>();
         File currentDir = new File(dir);
         for(String s : currentDir.list()) {
             if(s.endsWith(suffix)) {
                 files.add(s);
             }
         }
         return files;
    }
    public void output() throws IOException {
        /* 
        PrintWriter resultFile = new PrintWriter(new BufferedWriter(new FileWriter(percentage + "percent"+intensityThreshold+ms2FileName+".ratio")));
        PrintWriter problemFile = new PrintWriter(new BufferedWriter(new FileWriter(percentage + "percent"+intensityThreshold+ms2FileName+".problem")));
        resultFile.println("ScanNumber\tIntensity\tSignal2NoiseRatio\tzScore\tmean\tsigma\tdeltaMass(ppm)");
        resultFile.close();
        problemFile.close();
        */
    }
    private static void processMs1File(String ms1FileName, double minM2z, double maxM2z, int percentage, 
                double minRT, double maxRT, String fileformat) throws IOException {
        //SpectrumReader sr = new SpectrumReader(ms1FileName, "ms1");
        SpectrumReader sr = new SpectrumReader(ms1FileName, fileformat);
        int numSpectra = 0;
        double totalIntensity = 0;
        double totalSignal2Noise = 0;

        for(Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList pl = it.next();
            double retentionTime = pl.getRetentionTime();
            if(pl.numPeaks() > 0 && retentionTime >= minRT && retentionTime <= maxRT) {
                numSpectra++; 
                double intensity =pl.getMaxIntensity(minM2z, maxM2z); 
                double noiseIntensity = pl.getAvgIntensityOfLeastIntensePeaks(percentage);
//System.out.println(intensity + "\t" + noiseIntensity);
                totalIntensity += intensity;
                totalSignal2Noise += (intensity/noiseIntensity);
            }
        }
        System.out.println(ms1FileName + "\t" + totalIntensity/numSpectra + "\t" + totalSignal2Noise/numSpectra + "\t" + numSpectra);
    }

    public static void main(String args[]) throws Exception {
        try {
           double m2z = Double.parseDouble(args[0]); 
           double tolerance = Double.parseDouble(args[1]); 
           int percentage = Integer.parseInt(args[2]);
            
           double minRetentionTime = Double.parseDouble(args[3]); 
           double maxRetentionTime = Double.parseDouble(args[4]); 
           String filetype = args[5]; 
           List<String> files = getFiles(".", "ms1"); // get ms1 files
           Collections.sort(files);
           double minM2z = m2z - tolerance;
           double maxM2z = m2z + tolerance;
           System.out.println("SpectrumFile\tAvgIntensity\tAvgSignal2Noise\tNumSpectra");
           for(String file : files) {
               processMs1File(file, minM2z, maxM2z, percentage, minRetentionTime, maxRetentionTime, filetype);
           } 
        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("Finished");
    }
}
