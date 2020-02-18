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
 * @version $Id: PeptideIdComparator.java,v 1.4 2008/10/24 22:38:25 taoxu Exp $
 */

import edu.scripps.pms.mspid.MassCalculator;

public class PeptideIdComparator {
    private static ArrayList<SQTPeptide> sqts1 = new ArrayList<SQTPeptide>(1000000);
    private static ArrayList<SQTPeptide> sqts2 = new ArrayList<SQTPeptide>(1000000);
    public static final String USAGE = "java PeptideIdComparator sqt1 sqt2";
    public static final double XCORRDIFF = 0.4;
    private static int [] firstBetterFreq = new int[256];
    private static int [] secondBetterFreq = new int[256];

    public static void main(String args[]) throws IOException {
         if(args.length != 2) {
             System.err.println(USAGE);
             System.exit(0);
         }
         
         //addSqts(args[0]); 
         //addSqts(args[1]); 

         HashMap<Integer, SQTPeptide> scan2Sqt1 = addSqts(args[0]);
         HashMap<Integer, SQTPeptide> scan2Sqt2 = addSqts(args[1]);

         HashSet<Integer> scans = new HashSet(10000); 
         
         for(Iterator<Integer> it = scan2Sqt1.keySet().iterator(); it.hasNext();) {
             scans.add(it.next());
         }
         for(Iterator<Integer> it = scan2Sqt2.keySet().iterator(); it.hasNext();) {
             scans.add(it.next());
         }
         StringBuffer sb = new StringBuffer(1000000); 
         sb.append(args[0] + "\t" + args[1] + "\n");
         sb.append("NumSqt in first sqt: " + scan2Sqt1.size() + "\n");
         sb.append("NumSqt in second sqt: " + scan2Sqt2.size() + "\n");
         sb.append("ScanNum\tSequence\tXcorr\tDeltaCN\tSequence\tXcorr\tDeltaCN\tDeltaXcorr\n");
         for(Iterator<Integer> it = scans.iterator(); it.hasNext();) {
             Integer scanNum = it.next();
             SQTPeptide s = scan2Sqt2.get(scanNum); 
             SQTPeptide s1 = scan2Sqt1.get(scanNum);
             sb.append(scanNum.intValue() + "\t");
             if(s1 != null) {
                 sb.append(s1.getTopHit().getSequence() + "\t");
                 sb.append(s1.getTopHit().getXcorr() + "\t");
                 sb.append(s1.getTopHit().getDeltCN() + "\t");
             } else {
                 sb.append("\t\t\t");
             }
             if(s != null) {
                 sb.append(s.getTopHit().getSequence() + "\t");
                 sb.append(s.getTopHit().getXcorr() + "\t");
                 sb.append(s.getTopHit().getDeltCN());
             } else {
                 sb.append("\t\t\t");
             }
             if(s1 != null && s != null) {
                 double xcorrdiff = s1.getTopHit().getXcorrValue() - s.getTopHit().getXcorrValue();
                 sb.append("\t" + xcorrdiff);
                 String seq1 = s1.getTopHit().getSequence();
                 String seq2 = s.getTopHit().getSequence();
                 if(seq1.equals(seq2)) {
                     if(xcorrdiff > XCORRDIFF) {
                         calcResidueFreq(seq1, firstBetterFreq);
                     } else if(xcorrdiff < -XCORRDIFF) {
                         calcResidueFreq(seq1, secondBetterFreq);
                     }
                 }
             }
             sb.append("\n");

         }
         System.out.println("First better");
         outputResidueFreq(firstBetterFreq);
         System.out.println("Second better");
         outputResidueFreq(secondBetterFreq);
         System.out.println(sb.toString()); 
    }
    private static void outputResidueFreq(int [] freq) {
        double total = 0;
        for(int count : freq) {
            total += count;
        }
        String residues = "";
        String freqs = "";
        for(int i = 65; i <=90; i++) {
            residues += "\t" + (char)i;
            freqs += "\t" + freq[i]/total;
        }
        System.out.println(residues);
        System.out.println(freqs);

    }
    private static void calcResidueFreq(String seq, int [] freq) {
        for(int i = 2; i < seq.length()-2; i++) {
            freq[seq.charAt(i)]++;
        }
    } 
    public static  HashMap<Integer, SQTPeptide> addSqts(String sqtFile) throws IOException {
         HashMap<Integer, SQTPeptide> scan2Sqt = new HashMap<Integer, SQTPeptide>();
         SQTParser parser = new SQTParser(sqtFile);
         for(Iterator<SQTPeptide> itr = parser.getSQTPeptide(); itr.hasNext();) {
             SQTPeptide s = itr.next();
              
             if(s != null) { // for all
                 //if(s.getNumMlines() > 0) {
                 if(s.getNumMlines() > 0 && s.getTopHit("17PM") != null) {
                     Integer scanNum = new Integer(Integer.parseInt(s.getLoScan()));
                     SQTPeptide oldsqt = scan2Sqt.get(scanNum);
                     
                     if(oldsqt == null || oldsqt.getTopHit().getXcorrValue() < s.getTopHit().getXcorrValue()) {
                         scan2Sqt.put(scanNum, s);
                     } 
                        
                 }
             }
        }
        return scan2Sqt;
    }
}
