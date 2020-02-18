/**
 * @file AceMudpitFigure2.java
 * This is the source file for edu.scripps.pms.util.spectrum.AceMudpitFigure2
 * @author Tao Xu
 * @date $Date: 2006/07/25 22:44:22 $
 */



import edu.scripps.pms.util.io.*;

//import edu.scripps.pms.util.io.DTASelectFilterReader; 
import edu.scripps.pms.util.dtaselect.Peptide; 
import edu.scripps.pms.util.dtaselect.Protein; 
import edu.scripps.pms.util.PmsUtil; 
import java.io.*;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import java.util.HashSet;
import edu.scripps.pms.util.spectrum.*;

// this is the class for Akira's orthogonality check, use the pepdites identified in both
// experiments
public class AceMudpitFigure2 {
    public static final String USAGE = "\n\n!!! USAGE: aceMudpitFigure2 folder1 folder2 yes(withRedundant, no for without redundant) yes (commonPeptides only, no for all peptides)!!!";
    private static int minStep = 1000;
    private static int maxStep = 0;
    private static ArrayList<HashSet<String>> peptideSets = new ArrayList<HashSet<String>>();
    private static ArrayList<String> resultOfDir1;
    private static ArrayList<String> resultOfDir2;
    private static List<String> dirs1; // subdirectories for experiment 1
    private static List<String> dirs2;  // subdiretories for experiment 2
    private static HashSet<String> peptidesInCommon = new HashSet<String>();
    private static boolean withRedundant = false;
    private static boolean commonPeptidesOnly = true;
//    private double [] xcorrs = new double[100000];
//    private double [] deltaMasses = new double[100000];
    
    public static List<String> getFiles(String dir) {
        
         ArrayList<String> files = new ArrayList<String>();
         File currentDir = new File(dir);
         for(String s : currentDir.list()) {
             if(s.startsWith("Step_")) {
                 files.add(dir+ "/" + s);
                 //System.out.println(s);
             }
         }
         System.out.println("number of file for " + dir + ": " + files.size());      
         return files;
    }

    public static void output() throws IOException {
        HashSet<String> redundant = new HashSet<String>(1000000);
        int saltStep = 1;
        int sizeBefore = redundant.size();
        System.out.println("\n\nSaltStep\tNumberOfUniquePeptideAdded\tAccumulatedNumberOfPeptide");
        for(HashSet<String> peptideSet : peptideSets) {
            redundant.addAll(peptideSet);
            int sizeAfter = redundant.size();
            System.out.println(saltStep++ + "\t" + (sizeAfter - sizeBefore) + "\t" + sizeAfter);
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
    private static HashSet<String> getPeptideSet(List<String> dirs) {

        HashSet<String> peptideSet = new HashSet<String>(100000);
        
        for(String subdir : dirs) {

            String dtaselectFilterFile = subdir + "/DTASelect-filter.txt";
            try {
                DTASelectFilterReader reader = new DTASelectFilterReader(dtaselectFilterFile);
                for (Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {
                    Protein protein = itr.next();
 
                    for(Iterator<Peptide> pepItr=protein.getPeptides(); pepItr.hasNext(); ) {
                        Peptide peptide = pepItr.next();
                        peptideSet.add(peptide.getSequence());
                        //Integer lowSan = new Integer(peptide.getLoScan());
            //System.out.println(peptide.getSequence());

                    }
                }
            } catch(IOException e) {
        
            }
        }
        return peptideSet;
    }
    private static void getPeptidesInCommon() {
        HashSet<String> peptideSet1 = getPeptideSet(dirs1); 
        HashSet<String> peptideSet2 = getPeptideSet(dirs2); 
        System.out.println("Number of peptide in experiment 1: " + peptideSet1.size());
        System.out.println("Number of peptide in experiment 2: " + peptideSet2.size());

        peptideSet1.retainAll(peptideSet2);    
        peptidesInCommon = peptideSet1;
        System.out.println("Number of peptide in common: " + peptidesInCommon.size());
    }
    private static ArrayList<String> getMs2FileName(String parentDir, String subdir) {
         ArrayList<String> ms2Files = new ArrayList<String>();
         File currentDir = new File(subdir);
         for(String s : currentDir.list()) {
 //System.out.println(s);
             if(s.endsWith(".sqt")) {
                 int index = s.indexOf(".sqt");
                 String ms2FileName = parentDir + "/" + s.substring(0, index) + ".ms2";
     //            System.out.println(ms2FileName);
                 ms2Files.add(ms2FileName); 
                 //System.out.println("ms2File in " + subdir + "is: " + ms2FileName);
             }
         }
         return ms2Files;
    }
    private static double [] getRetentionTime(String ms2FileName) throws IOException {
        double [] retentionTime = new double[100000];
         
        SpectrumReader sr = new SpectrumReader(ms2FileName, "ms2");
        for(Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList pl = it.next();
            retentionTime[pl.getLoscan()] = pl.getRetentionTime();
        //    System.out.println("Scan#: " + pl.getLoscan() + "\tretentionTime: " + pl.getRetentionTime());
        }
        return retentionTime;
    }
    private static ArrayList<String> outputExperiment(String dir, List<String> subdirs) throws IOException {
    //    System.out.println(dir);
        HashSet<String> peptideAdded = new HashSet<String>(100000);
        HashSet<String> result = new HashSet<String>(100000);
        
        for(String subdir : subdirs) {
            System.out.println("Processing subdirectory " + subdir);
            String [] content = subdir.split("_");
            int saltStep = Integer.parseInt(content[content.length-1]);
            ArrayList<String> ms2Files = getMs2FileName(dir, subdir); 
            double [] retentionTime = null; 
            int numMs2File = ms2Files.size();
//System.out.println("numMs2File: " + numMs2File);
            HashMap<String, double []> retentionTimes = new HashMap<String, double []>();
            for(Iterator<String> it = ms2Files.iterator(); it.hasNext();) {
                String ms2File = it.next();
                retentionTime = getRetentionTime(ms2File);
//System.out.println("ms2File: " + ms2File);
                if(numMs2File > 1) {
                    int slashIndex = ms2File.lastIndexOf("/") + 1;
                
                    slashIndex = slashIndex > 0? slashIndex : 0;
                    String fileNameWithoutMs2 = ms2File.substring(slashIndex, ms2File.lastIndexOf("."));
                    retentionTimes.put(fileNameWithoutMs2, getRetentionTime(ms2File)); 
                }
            }
            //double [] retentionTime = getRetentionTime(ms2File);
            String dtaselectFilterFile = subdir + "/DTASelect-filter.txt";
            try {
                DTASelectFilterReader reader = new DTASelectFilterReader(dtaselectFilterFile);
                for (Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {
                    Protein protein = itr.next();
 
                    for(Iterator<Peptide> pepItr=protein.getPeptides(); pepItr.hasNext(); ) {
                        Peptide peptide = pepItr.next();
                        String sequence = peptide.getSequence();
                        int lowScan = Integer.parseInt(peptide.getLoScan());
//System.out.println("ms2File for this peptide is: " + ms2File);
                        //String seq = sequence.substring(2, sequence.length() -2);
                        float pI = peptide.getPi();
                        if(numMs2File > 1) { 
                            String ms2File = peptide.getFileName();
                            retentionTime = retentionTimes.get(ms2File);
                        }
                        String output = sequence + "\t" +  peptide.getScanNum() + "\t" + saltStep + "\t" + retentionTime[lowScan] + "\t" + peptide.getChargeState() + "\t" + pI;
                        if(commonPeptidesOnly) {
                            if(peptidesInCommon.contains(sequence)) {
                               if(withRedundant) {
                                   result.add(output);
                               } else {// show only peptides that are not previously identified
                                   if(!peptideAdded.contains(sequence)) {  
                                       result.add(output);
                                       peptideAdded.add(sequence);
                                   }
                               }
                            }
                        } else {
                           if(withRedundant) {
                               result.add(output);
                           } else {// show only peptides that are not previously identified
                               if(!peptideAdded.contains(sequence)) {  
                                   result.add(output);
                                   peptideAdded.add(sequence);
                               }
                           }

                        }
                    }
                }
            } catch(IOException e) {
                // in case the DATSelect-filter.txt does not exist
                System.out.println("File " + dtaselectFilterFile + " does not exist");
                //e.printStackTrace();
                //System.exit(1); 
            }

        }
        ArrayList<String> results = new ArrayList<String>(100000);
        for(Iterator<String> it = result.iterator(); it.hasNext();) {
           results.add(it.next()); 
        }
        return results;

    }
    public static boolean checkBooleanOption(String option) {
       char t = option.charAt(0);
       if(t == 'T' || t == 'Y' || t == 't' | t == 'y') {
           return true;
       } else {
           return false;
       } 

    }
    public static void main(String args[]) throws Exception {
        try {
           //System.out.println("Ms1File\tAvgIntensity\tAvgSignal2Noise\tNumSpectra");
   
           dirs1 = getFiles(args[0]);  
           dirs2 = getFiles(args[1]);  
           withRedundant = checkBooleanOption(args[2]);
           commonPeptidesOnly = checkBooleanOption(args[3]); 
           if(commonPeptidesOnly) {
               getPeptidesInCommon();
           }
           resultOfDir1 = outputExperiment(args[0], dirs1);
           resultOfDir2 = outputExperiment(args[1], dirs2);
           ArrayList<String> moreList = null;
           ArrayList<String> lessList = null;
           int numMore = 0;
           int numLess = 0;
           String header = null;
       
           String columns =  "\nPepitde\tScanNumber\tSaltStep\tRetentionTime\tChargeState\tpI" +
                             "\tPeptide\tScanNumber\tSaltStep\tRetentionTime\tChargeState\tpI\n";
           if(resultOfDir1.size() > resultOfDir2.size()) {
               moreList = resultOfDir1;
               lessList = resultOfDir2;
               numMore = resultOfDir1.size();
               numLess = resultOfDir2.size();
               header = args[0] + "\t\t\t\t\t" + args[1] + columns; 

           } else {
               moreList = resultOfDir2;
               lessList = resultOfDir1;
               numMore = resultOfDir2.size();
               numLess = resultOfDir1.size();
               header = args[1] + "\t\t\t" + args[0] + columns; 

           }
           System.out.println(header);
           int i = 0;
           for(String s : moreList) {
               System.out.print(s);
               if(i < numLess) {
                   System.out.print("\t" + lessList.get(i));
               }
               System.out.println(); 
               i++;

           } 

        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
            System.out.println(USAGE);
        }
    }
}
