import edu.scripps.pms.util.io.SQTParser;
//import java.io.BufferedReader;
import java.util.*;
import java.io.IOException;
import java.io.*;
import edu.scripps.pms.util.sqt.SQTPeptide;
import edu.scripps.pms.util.sqt.SQTPeptideSpComparator;
import edu.scripps.pms.util.sqt.MLine;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Tao Xu 
 * @version $Id: PeakFrequency.java,v 1.1 2006/06/25 07:37:17 taoxu Exp $
 */

import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.io.SpectrumReader;
// for Johaness's heatmap for frequency of peaks in give m2z bins 
public class PeakFrequency {
    public static final String USAGE = "java PeakFrequency [path] numBinsPerAMU";
    public static final int MAXPRCM2Z = 1400;
    public static final int MINPRCM2Z = 400;
    public static final int MAXFRAGM2Z = 2000;
    public static final int MINFRAGM2Z = 150;
    public static final String DELIMITER = "\t";
    private String scoreFile;
    private int numSpectra = 0; 
    private int numBinsPerAmu = 0;
    //private double [][] intensities; 
    private double [][] intensities; 
    private int minPrcIndex;
    private int maxPrcIndex;
    private int minFragIndex;
    private int maxFragIndex;
    private int numZeros = 0; 
    private int numNonZeros = 0; 
    private int numPeaks = 0; 
    private int [] prcFreq;   
    private int fragM2zMax = 0;    
    public PeakFrequency(int numBinsPerAmu) {
        this.numBinsPerAmu = numBinsPerAmu;
 
        maxPrcIndex = (MAXPRCM2Z-MINPRCM2Z)*numBinsPerAmu; 
        maxFragIndex = (MAXFRAGM2Z-MINFRAGM2Z)*numBinsPerAmu; 
        prcFreq = new int[maxPrcIndex]; 
        intensities = new double[maxPrcIndex][maxFragIndex];
        
        minPrcIndex = MINPRCM2Z*numBinsPerAmu;
        minFragIndex = MINFRAGM2Z*numBinsPerAmu;
        
    }

    public static void main(String args[]) throws IOException {
         ArrayList<String> fileNames = new ArrayList<String>();
         String path = ".";  
         int numBinsPerAmu = 0;
System.out.println("number of args: " + args.length);
         try {
             if(args.length == 2) {
                 path = args[0]; 
             
                 numBinsPerAmu = Integer.parseInt(args[1]);
             //System.err.println(USAGE);
             //System.exit(0);
             } else if(args.length == 1) {
                 numBinsPerAmu = Integer.parseInt(args[0]);
             } else {
                 System.out.println(USAGE);
                 System.exit(0);

             }
         } catch(Exception e) {

             e.printStackTrace();
             System.out.println(USAGE);
             System.exit(0);
         }
         
         PeakFrequency pf = new PeakFrequency(numBinsPerAmu);
         fileNames.addAll(getMs2Files("."));
         for(String file : fileNames) {
            System.out.println("Now processing " + file + "...");
            pf.addSpectra(file); 
         }
         
         pf.output();
    }
    public static List<String> getMs2Files(String dir) {
        
         ArrayList<String> ms2Files = new ArrayList<String>();
         File currentDir = new File(dir);
         for(String s : currentDir.list()) {
             if(s.endsWith(".ms2")) {
                 ms2Files.add(s);
             }
         }
         return ms2Files;
    }
    public void output() throws IOException {
        //System.out.println("\n\nReverseDeltaMass\tBothRightDeltaMAss\tprimaryRightDeltaMass\tspRightDeltaMass\n");
        PrintWriter outFile = new PrintWriter(new BufferedWriter(new FileWriter("m2zFrequency.txt")));
        //outFile.println("Minimum PrcM2z: " + MINPRCM2Z + "\tMAximum PrcM2z: " + MAXPRCM2Z + "\tMinimum FragM2z: " + MINFRAGM2Z + "\tMaximum FragM2z: " + MAXFRAGM2Z);
        outFile.print(minFragIndex/numBinsPerAmu);
        for(int i = 0; i < fragM2zMax; i++) {
            outFile.print(DELIMITER + (i+minFragIndex+0.0)/numBinsPerAmu);
        }

        outFile.println("");
        for(int i = 0; i < maxPrcIndex; i++) {
            outFile.print((i+minPrcIndex+0.0)/numBinsPerAmu + "(" + prcFreq[i] + ")"); 
            //for(int j = 0; j < maxFragIndex; j++) {
            for(int j = 0; j < fragM2zMax; j++) {
                outFile.print(DELIMITER + (int)(intensities[i][j]+0.5));
                if(intensities[i][j] == 0) {
                    numZeros++;
                } else { numNonZeros++; } 
            }
            outFile.print("\n");
        } 
        System.out.println("number of spectra: " + numSpectra + "\tnumber of Zeros: " + numZeros);
        System.out.println("number of peaks: " + numPeaks + "\tnumber of non-Zeros: " + numNonZeros);
        System.out.println("Minimum PrcM2z: " + MINPRCM2Z + "\tMAximum PrcM2z: " + MAXPRCM2Z + "\tMinimum FragM2z: " + MINFRAGM2Z + "\tMaximum FragM2z: " + MAXFRAGM2Z);
        System.out.println("Maximu fragement m2z: " + fragM2zMax);
    }
    public void addSpectra(String ms2FileName) throws IOException {

         SpectrumReader sr = new SpectrumReader(ms2FileName, "ms2");
         for (Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList pl = it.next();        
            numSpectra++;
            int prcM2zIndex = (int)((pl.getPrecursorMass()-MINPRCM2Z)*numBinsPerAmu+0.5);
            if(prcM2zIndex < maxPrcIndex && prcM2zIndex >= 0) {
                prcFreq[prcM2zIndex]++;
                for(Iterator<Peak> itp = pl.getPeaks(); itp.hasNext();) {
                    Peak p = itp.next();
                    int m2zIndex = (int)((p.getM2z()-MINFRAGM2Z)*numBinsPerAmu + 0.5) - minFragIndex;
                    if(m2zIndex < maxFragIndex && m2zIndex >= 0) {
                        //intensities[prcM2zIndex][m2zIndex] += p.getIntensity();   
                        fragM2zMax = fragM2zMax < m2zIndex? m2zIndex : fragM2zMax;
                        intensities[prcM2zIndex][m2zIndex]++;   
                        numPeaks++;
                    }
                
                }
            }
         }
    }
}
