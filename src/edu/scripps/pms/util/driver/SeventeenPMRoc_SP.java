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
 * @version $Id: SeventeenPMRoc_SP.java,v 1.2 2008/01/17 19:46:57 taoxu Exp $
 */

import edu.scripps.pms.mspid.MassCalculator;

public class SeventeenPMRoc_SP {
    ArrayList<SQTPeptide> sqts = new ArrayList<SQTPeptide>(1000000);
    public static final String USAGE = "java SevenTeenPMRoc 17PM_ chargeState interval, use 0 as chargeState for all chargeState";
    public static final String contaminant = "contaminant";
    public static final String REVERSE = "Reverse_";
    private int chargeState;
    private int interval;
    private String locus;
    private String reverseLocus;
    private String scoreFile;
    private double totalCpuTime = 0;
    private int totalNumSpectra = 0;
    private double totalNumCandidate = 0;
    private double numSpectraWithValidCpuTime = 0; 
    public static final int MINPEPLENGTH = (7+4)-1; // minimum peptide length 7 - 1 to avoid >=

    public static void main(String args[]) throws IOException {
         ArrayList<String> fileNames = new ArrayList<String>();
         if(args.length != 3) {
             System.err.println(USAGE);
             System.exit(0);
         }
         String locus = args[0];
         
         int charge = Integer.parseInt(args[1]);
         int interval = Integer.parseInt(args[2]);// determine how ofen to output fpf and tpf
         SeventeenPMRoc_SP roc = new SeventeenPMRoc_SP(locus, charge, interval);
         fileNames.addAll(getSqtFiles("."));
         for(String file : fileNames) {
            System.out.println("Now processing " + file + "...");
            roc.addSqts(file); 
         }
         roc.sortSqtPeptides();
         roc.outputRoc();
         
    }
    public SeventeenPMRoc_SP(String loc, int charge, int interval) {
        locus = loc;
        reverseLocus = "Reverse_" + locus;
        chargeState = charge;
        this.interval = interval;
        scoreFile = "xscores_" + charge + ".txt";
    }
    public void outputRoc() throws IOException {
        int numTruePositive = 0;
        int numFalsePositive = 0;
        int numReverse10PM = 0;
        int numTotal = 0;
        for(SQTPeptide s : sqts) {
            //if(s.topHitStartsWith(locus) || s.topHitStartsWith(contaminant)) {
            if(s.topHitStartsWith(locus)) {
                numTruePositive++;
            } else if (s.topHitStartsWith(REVERSE)){
                numFalsePositive++;
                if(s.topHitStartsWith(reverseLocus)) {
                    numReverse10PM++;
                }
            }
            numTotal++;
        }
        
        double [] scoreOfTrue = new double[numTruePositive];
        double [] deltaCnOfTrue = new double[numTruePositive];
        double [] deltaMassOfTrue = new double[numTruePositive];

        double [] scoreOfFalse = new double[numFalsePositive];
        double [] deltaCnOfFalse = new double[numFalsePositive];
        double [] deltaMassOfFalse = new double[numFalsePositive];
        int [] scanNumbers = new int[numFalsePositive];
        int numT = 0;
        int numF = 0;
        double fp = 0;
        double tp = 0;
        int numElements = numTotal/interval + 2;
//System.err.println("numElements: " + numElements);
        double [] tpfs = new double[numElements];
        double [] fpfs = new double[numElements];
        int index = 0;
        System.out.println("numTrueHits\tFPF\tTPF\tTPR\tFPR\tSpScore");
        int arrayIndex = 0;
        
        for(SQTPeptide s : sqts) {
                //if(s.topHitStartsWith(locus) || s.topHitStartsWith(contaminant)) {
                if(s.topHitStartsWith(locus)) {
                    //scoreOfTrue[numT] = s.getTopHit().getXcorrValue(); 
                    scoreOfTrue[numT] = s.getTopHit().getSpValue(); 
                    deltaCnOfTrue[numT] = s.getDeltaCn();
                    deltaMassOfTrue[numT] = s.getDeltaMass();
                    //deltaCnOfTrue[numT] = s.getTScore();
                    numT++;
                } else if (s.topHitStartsWith(REVERSE)){
                    //scoreOfFalse[numF] = s.getTopHit().getXcorrValue(); 
                    scoreOfFalse[numF] = s.getTopHit().getSpValue(); 
                    deltaCnOfFalse[numF] = s.getDeltaCn();
                    deltaMassOfFalse[numF] = s.getDeltaMass();
                    //deltaCnOfFalse[numF] = s.getTScore();
                    scanNumbers[numF] = Integer.parseInt(s.getLoScan());
                    numF++;
                }
                   
                if(index%interval == 0) {
                    fp = numF/(numF+numT+0.0);
                    tp = numT/(numT+numF+0.0);
                    tpfs[arrayIndex] = numT/(numTruePositive+0.0);
                    fpfs[arrayIndex] = numF/(numFalsePositive+0.0);
                    //System.out.println(fpfs[arrayIndex] + "\t" + tpfs[arrayIndex] + "\t" + tp + "\t" + fp+ "\t"+s.getTopHit().getXcorr());
                    System.out.println(numT + "\t" + fpfs[arrayIndex] + "\t" + tpfs[arrayIndex] + "\t" + numT + "\t" + tp + "\t" + fp+ "\t"+s.getTopHit().getXcorr());
                    arrayIndex++;

                }
                index++;
        }
        
        tpfs[arrayIndex] = numT/(numTruePositive+0.0);
        fpfs[arrayIndex] = numF/(numFalsePositive+0.0);
        fp = numF/(numF+numT+0.0);
        tp = numT/(numT+numF+0.0);
        System.out.println(numT + "\t" + fpfs[arrayIndex] + "\t" + tpfs[arrayIndex] + "\t" + tp + "\t" + fp);
        double areaUnderCurve = 0;
        System.out.println(numF/(numFalsePositive+0.0) + "\t" + numT/(numTruePositive+0.0));
        for(int i = 1; i < tpfs.length; i++) {
            double avg = (tpfs[i] + tpfs[i-1])/2;
            areaUnderCurve += (avg * (fpfs[i] - fpfs[i-1]));    
//System.out.println("A: " + areaUnderCurve + "\ttpfs: " + tpfs[i] + "\tavg: " + avg + "\tfpf[i]: " + fpfs[i] + "\tfpfs[i=1]: " + fpfs[i-1]);
        } 
        
        System.out.println("Area Under the ROC Curve: \t" + areaUnderCurve);
        System.out.println("num true positive: " + numT + "\tnum false positive: " + numF + "\tnum " + reverseLocus + ": " + numReverse10PM);
        System.out.println("Total number of spectra searched: " + totalNumSpectra + "\tCPU time in milliseconds: " + totalCpuTime + "\taverage time used per spectrum: " + totalCpuTime/numSpectraWithValidCpuTime);
        System.out.println("Avg number of candidate peptide per spectrum: " + totalNumCandidate/totalNumSpectra);
        // output true and false positive scores 
        PrintWriter outFile = new PrintWriter(new BufferedWriter(new FileWriter(scoreFile)));
        int maxSize = Math.max(numTruePositive, numFalsePositive);
       // System.out.println("output scores, maxSize: " + maxSize); 
        outFile.println("TrueScore\tTrueDeltaCn\tTrueDeltaMass\tFalseScore\tFalseDeltaCn\tFalseDeltaMass");
        for(int i = 0; i < maxSize; i++) {
            if(i < numTruePositive) {
                outFile.print(scoreOfTrue[i] + "\t" + deltaCnOfTrue[i] + "\t" + deltaMassOfTrue[i]);
            } else {

                outFile.print(" \t \t ");
            }
            outFile.print("\t");
            if(i < numFalsePositive) {
                outFile.print(scoreOfFalse[i] + "\t" + deltaCnOfFalse[i] + "\t" + deltaMassOfFalse[i] + "\t" + scanNumbers[i]);
            } else {

                outFile.print(" \t ");
            }
            outFile.println();
            outFile.print("");
            
        }
        outFile.close(); 
        
    }
    public void sortSqtPeptides() {

        System.out.println("number of sqts: " + sqts.size()); 
        int scoreType = SQTPeptideSpComparator.SORTBYSP;
        //int scoreType = SQTPeptideSpComparator.SORTBYXCORR;
        SQTPeptideSpComparator comparator = new SQTPeptideSpComparator(scoreType);
        Collections.sort(sqts, comparator);
    }
    public void addSqts(String sqtFile) throws IOException {

         SQTParser parser = new SQTParser(sqtFile);
         for(Iterator<SQTPeptide> itr = parser.getSQTPeptide(); itr.hasNext();) {
             SQTPeptide s = itr.next();
              
             if(s != null) { // for all
                 if(chargeState != 0 && s.getChargeStateInt() != chargeState) {
                     continue;
                 } 
    //         if(s != null && s.getChargeStateInt() == chargeState) { // for each charge state
                 totalNumSpectra++;
                 totalNumCandidate += Integer.parseInt(s.getNumSeq());
                 int cpuTime = Integer.parseInt(s.getTimeToProcess()); 
                 if(cpuTime > 0) { 
                     totalCpuTime += cpuTime; 
                     numSpectraWithValidCpuTime++;
                 }
                 //if(s.getNumMlines() > 0) {
                 if(s.getNumMlines() > 0 && s.getTopHit().getSequence().length() > MINPEPLENGTH) {
                     sqts.add(s);
                 }
             }
         }
    }
    public static List<String> getSqtFiles(String dir) {
        
         ArrayList<String> sqtFiles = new ArrayList<String>();
         File currentDir = new File(dir);
         for(String s : currentDir.list()) {
             if(s.endsWith(".sqt")) {
                 sqtFiles.add(s);
             }
         }
         return sqtFiles;
    }
}
