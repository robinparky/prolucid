/**
 * @file AceMudpitFigure1.java
 * This is the source file for edu.scripps.pms.util.spectrum.AceMudpitFigure1
 * @author Tao Xu
 * @date $Date: 2008/10/24 22:38:25 $
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
import java.util.HashSet;
import edu.scripps.pms.util.spectrum.*;

// this is the class for Akira's aceMudpit paper, it should be used after running
// makeSubDirs and runCommandInSubDirectories DTASelect2.  It will compute number
// of additional unique peptides in each salt step 
public class AceMudpitFigure1 {
    public static final String USAGE = "java AceMudpitFigure1";
    private static int minStep = 1000;
    private static int maxStep = 0;
    private static ArrayList<HashSet<String>> peptideSets = new ArrayList<HashSet<String>>();
    private static boolean withRedundantPeptide = false;
//    private double [] xcorrs = new double[100000];
//    private double [] deltaMasses = new double[100000];
    
    public static List<String> getFiles(String dir) {
        
         ArrayList<String> files = new ArrayList<String>();
         File currentDir = new File(dir);
         for(String s : currentDir.list()) {
             if(s.startsWith("Step_")) {
                 files.add(s);
                 //System.out.println(s);
             }
         }
         
         Collections.sort(files);
         return files;
    }
    public static void output() throws IOException {
        HashSet<String> redundant = new HashSet<String>(1000000);
        int saltStep = 1;
        int sizeBefore = redundant.size();
        System.out.println("\n\nSaltStep\tNumberOfUniquePeptideAdded\tAccumulatedNumberOfPeptide\tCarryOverPercentage\tNumPeptidesInThisStep");
        for(HashSet<String> peptideSet : peptideSets) {
            redundant.addAll(peptideSet);
            int numPeptides = peptideSet.size(); // num of peptide identified in this step
            int sizeAfter = redundant.size();
            int numUniquePeptides = sizeAfter - sizeBefore; // numUniquePeptidesInThisStep
            int numCarryOver = numPeptides - numUniquePeptides;
            double carryOverRate = (numCarryOver/(numPeptides+0.0))*100; 
            System.out.println(saltStep++ + "\t" + numUniquePeptides + "\t" + sizeAfter + "\t" + carryOverRate + "\t" + numPeptides);
            sizeBefore = sizeAfter;
        }
  
        //System.out.println("Total number of unique peptides: " + sizeBefore);
       
        /* 
        PrintWriter resultFile = new PrintWriter(new BufferedWriter(new FileWriter(percentage + "percent"+intensityThreshold+ms2FileName+".ratio")));
        PrintWriter problemFile = new PrintWriter(new BufferedWriter(new FileWriter(percentage + "percent"+intensityThreshold+ms2FileName+".problem")));
        resultFile.println("ScanNumber\tIntensity\tSignal2NoiseRatio\tzScore\tmean\tsigma\tdeltaMass(ppm)");
        resultFile.close();
        problemFile.close();
        */
    }
    private static void processSaltStep(String subdir) throws IOException {
        int saltStep = Integer.parseInt(subdir.split("_")[1]);        
        HashSet<String> peptides = new HashSet<String>(50000); 
        peptideSets.add(peptides);
        minStep = minStep < saltStep? minStep : saltStep;
        maxStep = maxStep > saltStep? maxStep : saltStep;
        String dtaselectFilterFile = subdir + "/DTASelect-filter.txt";
       // System.out.println("saltStep: " + saltStep + "minStep: " + minStep + "\tmaxStep: " + maxStep);
        try {
            DTASelectFilterReader reader = new DTASelectFilterReader(dtaselectFilterFile);
            for (Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {
                Protein protein = itr.next();

                for(Iterator<Peptide> pepItr=protein.getPeptides(); pepItr.hasNext(); ) {
                    Peptide peptide = pepItr.next();
                    peptides.add(peptide.getSequence());
                    //Integer lowSan = new Integer(peptide.getLoScan());
        //System.out.println(peptide.getSequence());

                }
            }
            //peptideSets.add(peptides);
            //System.out.println(peptides.size() + " peptides found in " + dtaselectFilterFile);
        } catch(IOException e) {
           System.out.println(dtaselectFilterFile + " file does not exist");
       //e.printStackTrace(); 
        }

    }

    public static void main(String args[]) throws Exception {
        try {
           //System.out.println("Ms1File\tAvgIntensity\tAvgSignal2Noise\tNumSpectra");
           List<String> files = getFiles(".");  
           for(String file : files) {
               processSaltStep(file);
           } 
           output();
        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
            System.out.println(USAGE);
        }
    }
}
