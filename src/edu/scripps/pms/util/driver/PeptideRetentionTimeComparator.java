/**
 * @file PeptideRetentionTimeComparator.java
 * This is the source file for edu.scripps.pms.util.spectrum.PeptideRetentionTimeComparator
 * @author Tao Xu
 * @date $Date
 */



import edu.scripps.pms.util.io.*;

//import edu.scripps.pms.util.io.DTASelectFilterReader; 
import edu.scripps.pms.util.dtaselect.Peptide; 
import edu.scripps.pms.util.dtaselect.Protein; 
import edu.scripps.pms.util.io.SpectrumReader; 
import java.io.*;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import java.util.HashSet;
import edu.scripps.pms.util.spectrum.*;

// comparing retention times of peptides identified in the sampe salt step in different MudPIT runs 
// to estimate the reproducibility of each salt step of LC/LC 
// this program expect a text file that contains all the directory names for comparison
// it also expect the number of salt steps the use wants to compare
public class PeptideRetentionTimeComparator {
    public static final String USAGE = "java PeptideRetentionTimeComparator directorynamefile numSaltStep";
    private static int minStep = 1000;
    private static int maxStep = 0;
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
    // assume only one ms2 file in the folder
    public static String getMs2File(String dir) {
         System.out.println("getting ms2 file for " + dir); 
         File currentDir = new File(dir);
         String [] files = currentDir.list();
         for(int i = 0; i < files.length; i++) {
             String s = files[i];
             if(s.endsWith("ms2")) {
                 return s;
             }
         }
         return null; 
    } 
    private static void processSaltStep(HashSet<String> dirs, int saltStep) throws IOException {

        System.out.println("Processing salt step " + saltStep);
        int numExp = dirs.size();

        ArrayList<HashMap<String, Double>> seq2Doubles = new ArrayList<HashMap<String, Double>>(); 
        ArrayList<HashMap<String, Peptide>> seq2Peptides = new ArrayList<HashMap<String, Peptide>>(); 

        for(String subdir : dirs) {
            //subdir = subdir + "/saltstep/";
            String steps = saltStep < 10? ("/Step_0" + saltStep) : ("/Step_" + saltStep);
            subdir = subdir + steps + "/";
            String dtaselectFilterFile = subdir + "/DTASelect-filter.txt";
           // System.out.println("saltStep: " + saltStep + "minStep: " + minStep + "\tmaxStep: " + maxStep);

            System.out.println("Processing " + dtaselectFilterFile);
            HashMap<String, Double> seq2RetentionTime = new HashMap<String, Double>();
            HashMap<String, Peptide> seq2Peptide = new HashMap<String, Peptide>();

            //SpectrumReader sr = new SpectrumReader(subdir+"/*.ms2", "ms2");
            String ms2file = subdir + getMs2File(subdir);
            SpectrumReader sr = new SpectrumReader(ms2file, "ms2");
            
            System.out.println("starting to get retention time array ");
            double [] retentiontimes = sr.getRetentionTimes(); 
            System.out.println("finished getting retention time array ");
            try {
                DTASelectFilterReader reader = new DTASelectFilterReader(dtaselectFilterFile);
                System.out.println("Processing " + dtaselectFilterFile);
                for (Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {
                    Protein protein = itr.next();

                    for(Iterator<Peptide> pepItr=protein.getPeptides(); pepItr.hasNext(); ) {
                        Peptide peptide = pepItr.next();
                        String seq = peptide.getSequence()+peptide.getChargeState();
                        double rettime = retentiontimes[peptide.getScanNumber()];
                        seq2RetentionTime.put(seq, new Double(rettime));                             
                        seq2Peptide.put(seq, peptide);                             
                    }
                }
                //peptideSets.add(peptides);
                //System.out.println(peptides.size() + " peptides found in " + dtaselectFilterFile);
            } catch(IOException e) {
               System.out.println(dtaselectFilterFile + " file does not exist");
           //e.printStackTrace(); 
            }
            seq2Doubles.add(seq2RetentionTime);
            seq2Peptides.add(seq2Peptide);
        }
        HashSet<String> seqs = new HashSet<String>(seq2Doubles.get(0).keySet());
        for(int i = 1; i < numExp; i++) {
            seqs.retainAll(seq2Doubles.get(i).keySet());
        }
        for(String s : seqs) {
            System.out.print(s + "\t");
            for(int i = 0; i < numExp; i++) {
                System.out.print(seq2Peptides.get(i).get(s).getRedundancy() + "\t");
            } 
            for(int i = 0; i < numExp; i++) {
                System.out.print(seq2Doubles.get(i).get(s).doubleValue() + "\t");
            } 
            System.out.println();
        }
    }

    public static HashSet<String> readDirs(String file) throws IOException {
        HashSet<String> dirs = new HashSet<String>();
        BufferedReader br = new BufferedReader(new FileReader(file), 4096);
        String line = null;
        while((line = br.readLine()) != null) {
            dirs.add(line.trim());    
        }
        br.close();
        return dirs; 
    }
    public static void main(String args[]) throws Exception {

        try {
           
           //System.out.println("Ms1File\tAvgIntensity\tAvgSignal2Noise\tNumSpectra");
           HashSet<String> dirs = readDirs(args[0]);
           System.out.println("Number of experiments: " + dirs.size());
           int numsaltstep = Integer.parseInt(args[1])+1;
           
           // ignore first salt step because it is often empty
           for(int i = 2; i < numsaltstep; i++) { // for each salt step
               processSaltStep(dirs, i);
           }
        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
            System.out.println(USAGE);
        }
    }
}
