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
 * @version $Id: SeventeenPMCounter.java,v 1.4 2007/09/10 20:25:41 taoxu Exp $
 */

import edu.scripps.pms.mspid.MassCalculator;

public class SeventeenPMCounter {
    ArrayList<SQTPeptide> sqts = new ArrayList<SQTPeptide>(1000000);
    public static final String USAGE = "java SevenTeenPMRoc 17PM_ chargeState 0 for all chargeStates";
    public static final String contaminant = "contaminant";
    public static final String REVERSE = "Reverse_";
    private int chargeState;
    private String locus;
    private String reverseLocus;
    private String scoreFile;
    private double totalCpuTime = 0;
    private int totalNumSpectra = 0;
    private double totalNumCandidate = 0;
    private double numSpectraWithValidCpuTime = 0; 
 
    private int numSpectraWith17pmFound = 0;
    private int numBothRight = 0;
    private int numXcorrRight = 0;
    private int numSpRight = 0;

    // the following three variables are for hits that neither primary rank
    // nor the sp rank is one
    private int numPrimaryRankBetter = 0;
    private int numSpRankBetter = 0;
    private int numTie = 0;
    
    private int numReverse17PM = 0;
    private int numTopReverse17PM = 0;
    private int totalReverse = 0;
    private int totalForward = 0;
    private double[] bothRightDeltaMass = new double[200000];
    private double [] primaryRightDeltaMass = new double[200000];
    private double [] spRightDeltaMass = new double[200000];
    private double [] reverseDeltaMass = new double [200000];
    private String [] reverseAcc = new String[200000];
    private int numReverseDeltaMass = 0; 
    public static void main(String args[]) throws IOException {
         ArrayList<String> fileNames = new ArrayList<String>();
         if(args.length != 2) {
             System.err.println(USAGE);
             System.exit(0);
         }
         String locus = args[0];
         
         int charge = Integer.parseInt(args[1]);
         SeventeenPMCounter roc = new SeventeenPMCounter(locus, charge);
         fileNames.addAll(getSqtFiles("."));
         for(String file : fileNames) {
            System.out.println("Now processing " + file + "...");
            roc.addSqts(file); 
         }
         
         roc.output();
    }
    public void output() {
        System.out.println("numBothRight:\t" + numBothRight);
        System.out.println("numPrimaryRankRight:\t" + numXcorrRight);
        System.out.println("numSpRankRight:\t" + numSpRight);
        System.out.println("numPrimaryRankBettert:\t" + numPrimaryRankBetter);
        System.out.println("numSpRankBetter:\t" + numSpRankBetter);
        System.out.println("numTie:\t" + numTie);
        System.out.println("numTopReverse17PM topHits:\t" + numTopReverse17PM);
        System.out.println("numReverse17PM:\t" + numReverse17PM);
        System.out.println("totalNumSpectraSearched:\t" + totalNumSpectra);
        int max = numBothRight > numXcorrRight? numBothRight : numXcorrRight;
        max = max > numSpRight? max : numSpRight;
        max = max > numReverseDeltaMass? max : numReverseDeltaMass;
        System.out.println("max: " + max + "\t numReverseDeltaMass: " +  numReverseDeltaMass);
        System.out.println("\n\nReverseDeltaMass\tBothRightDeltaMAss\tprimaryRightDeltaMass\tspRightDeltaMass\n");
        for(int i = 0; i < max; i++) {
            if(i < numReverseDeltaMass) {
                System.out.print(reverseDeltaMass[i] + "\t"); 
                System.out.print(reverseAcc[i] + "\t");
            } else {

                System.out.print(" \t"); 
            }
            if(i < numBothRight) {
                System.out.print(bothRightDeltaMass[i] + "\t"); 
            } else {

                System.out.print(" \t"); 
            }
            if(i < numXcorrRight) {
                System.out.print(primaryRightDeltaMass[i] + "\t"); 
            } else {

                System.out.print(" \t"); 
            }
            if(i < numSpRight) {
                System.out.print(spRightDeltaMass[i] + "\t"); 
            } else {
                System.out.print(" \t"); 
            }
            System.out.print("\n");
        } 
    }
    public SeventeenPMCounter(String loc, int charge) {
        locus = loc;
        reverseLocus = "Reverse_" + locus;
        chargeState = charge;
    }
    public void addSqts(String sqtFile) throws IOException {

         SQTParser parser = new SQTParser(sqtFile);
         for(Iterator<SQTPeptide> itr = parser.getSQTPeptide(); itr.hasNext();) {
             SQTPeptide s = itr.next();
              
             //if(s != null && s.getChargeStateInt() == chargeState) {
             if(s != null) {
                 if(chargeState != 0 &&  s.getChargeStateInt() != chargeState) {
                     continue;
                 }
                 totalNumSpectra++;
                 totalNumCandidate += Integer.parseInt(s.getNumSeq());
                 int cpuTime = Integer.parseInt(s.getTimeToProcess()); 
                 if(cpuTime > 0) { 
                     totalCpuTime += cpuTime; 
                     numSpectraWithValidCpuTime++;
                 }
              
                 if(s.getNumMlines() > 0) {
                     //MLine mTop = s.getTopHit();
                     MLine mSpTop = s.getSpRankOneHit();
                     if(s.topHitStartsWith(REVERSE)) {
                         MLine mTop = s.getTopHit(REVERSE);
                         reverseAcc[numReverseDeltaMass] = mTop.getFirstLLine();
                         reverseDeltaMass[numReverseDeltaMass++] = s.getDeltaMass(mTop);
// System.out.println(mTop.getFirstLLine());
                     }
                     if(s.topHitStartsWith(reverseLocus)) {
                         numTopReverse17PM++;
                     }
                    
                     if(s.contains(reverseLocus)) { 
                         numReverse17PM++;
                     }
                     if(s.contains(locus)) {
                         MLine mTop = s.getTopHit(locus);
                         if(mTop != null && mTop.contains(locus)) {
//System.out.println(mTop.getMLine() + "SpRank: " + mTop.getSpRankInt());
                             if(mTop.getSpRankInt() == 1) {
                                 bothRightDeltaMass[numBothRight] = s.getDeltaMass(mTop);
                                 numBothRight++;
                             } else {
                                 primaryRightDeltaMass[numXcorrRight] = s.getDeltaMass(mTop);
                                 numXcorrRight++;
                             } 
                         } else {
                             if(mSpTop.contains(locus)){
                                 spRightDeltaMass[numSpRight] = s.getDeltaMass(mSpTop);
                                 numSpRight++;
                             } else {
                                 MLine m = s.getMLineStartsWith(locus);
                                 int pRank = m.getPrimaryRankInt();
                                 int spRank = m.getSpRankInt();
                                 if(pRank < spRank) {
                                     numPrimaryRankBetter++;
                                 } else if(pRank > spRank) {
                                     numSpRankBetter++;
                                 } else {
                                     numTie++;
                                 }
                             }
                         }
                     }
                     
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
