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
 * @version $Id: CamsiFromSqt.java,v 1.1 2008/10/24 22:39:14 taoxu Exp $
 */

//import edu.scripps.pms.mspid.MassCalculator;

public class CamsiFromSqt {
    ArrayList<SQTPeptide> sqts = new ArrayList<SQTPeptide>(1000000);
    public static final String USAGE = "java CamsiFromSqt ProteinACCPrefix scoretype(sp or xc) deltamasscutoff";
    public static final String contaminant = "contaminant";
    public static final String REVERSE = "Reverse_";
    private int chargeState;
    private int interval;
    private String locus;
    private String reverseLocus;
    private double totalCpuTime = 0;
    private int totalNumSpectra = 0;
    private double totalNumCandidate = 0;
    private double numSpectraWithValidCpuTime = 0; 
    private static boolean getsp = false;
    private static double deltamasscutoff = 5; // 5 ppm
    public static void main(String args[]) throws IOException {
         ArrayList<String> fileNames = new ArrayList<String>();
         if(args.length < 2) {
             System.err.println(USAGE);
             System.exit(0);
         }
         String locus = args[0];
         if(args.length > 2) {
             deltamasscutoff = Double.parseDouble(args[2]);
         }
         getsp = "sp".equals(args[1]);
         CamsiFromSqt roc = new CamsiFromSqt(locus, 0); // get all chargestate
         fileNames.addAll(getSqtFiles("."));
         for(String file : fileNames) {
            //System.out.println("Now processing " + file + "...");
            roc.addSqts(file); 
         }
         roc.sortSqtPeptides();
         roc.outputRoc();
         
    }
    public CamsiFromSqt(String loc, int charge) {
        locus = loc;
        reverseLocus = "Reverse_" + locus;
        chargeState = charge;
        this.interval = interval;
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
        int numElements = 10000;
//System.err.println("numElements: " + numElements);
        double [] tpfs = new double[numElements];
        double [] fpfs = new double[numElements];
        System.out.println("FPF\tTPF\tTPR\tFPR\tPrimaryScore");
        int arrayIndex = 0;
        
        for(SQTPeptide s : sqts) {
                //if(s.topHitStartsWith(locus) || s.topHitStartsWith(contaminant)) {
                MLine mline = s.getTopHit();
                String peptideseq = mline.getSequence();
                Iterator<String> llines = mline.getLLine(); 
                //double score = s.getTopHit().getXcorrValue();
                double score = getScore(s); //.getTopHit().getXcorrValue();
                double deltamass = s.getDeltaMassInPpm();
                int scan = Integer.parseInt(s.getLoScan());
                boolean isReverseHit = false;
                if(deltamasscutoff > 0 && deltamass > deltamasscutoff) {
                    continue;
                }
                if(s.topHitStartsWith(locus)) {
                    scoreOfTrue[numT] = score; 
                    deltaCnOfTrue[numT] = s.getDeltaCn();
                    deltaMassOfTrue[numT] = deltamass; 
                    numT++;
                    
                } else if (s.topHitStartsWith(REVERSE)){
                    scoreOfFalse[numF] = score; 
                    deltaCnOfFalse[numF] = s.getDeltaCn();
                    deltaMassOfFalse[numF] = deltamass; 
                    numF++;
                    isReverseHit = true;
System.out.println(scan + "\t" + llines.next() + "\t" + score + "\t" + peptideseq + "\t" + fp + "\t" + deltamass); 
                }
                   
                fp = numF/(numF+numT+0.0);
                tp = numT/(numT+numF+0.0);
                tpfs[arrayIndex] = numT/(numTruePositive+0.0);
                fpfs[arrayIndex] = numF/(numFalsePositive+0.0);
                //System.out.println(fpfs[arrayIndex] + "\t" + tpfs[arrayIndex] + "\t" + tp + "\t" + fp+ "\t"+s.getTopHit().getXcorr());

                //System.out.println(fpfs[arrayIndex] + "\t" + tpfs[arrayIndex] + "\t" + numT + "\t" + tp + "\t" + fp+ "\t"+s.getTopHit().getXcorr());
                if(!isReverseHit) {
                    StringBuffer sb = new StringBuffer(500);
                    sb.append(scan);
                    sb.append("\t" + llines.next());

                    while(llines.hasNext()) {
                        sb.append(";" + llines.next());
                    }
                    sb.append("\t"); 
                    // may need to remove mass shift here
                    sb.append("\t" + peptideseq + "\t");
                    sb.append(score + "\t");
                    sb.append(fp);
                    sb.append("\t" + deltamass);
          
                    System.out.println(sb);
                }
                //System.out.println(fpfs[arrayIndex] + "\t" + tpfs[arrayIndex] + "\t" + numT + "\t" + tp + "\t" + fp+ "\t"+s.getTopHit().getXcorr());
                arrayIndex++;

        }
        
        tpfs[arrayIndex] = numT/(numTruePositive+0.0);
        fpfs[arrayIndex] = numF/(numFalsePositive+0.0);
        fp = numF/(numF+numT+0.0);
        tp = numT/(numT+numF+0.0);
        System.out.println(fpfs[arrayIndex] + "\t" + tpfs[arrayIndex] + "\t" + tp + "\t" + fp);
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
        
    }
    public void sortSqtPeptides() {

        System.out.println("number of sqts: " + sqts.size()); 
        int scoreType = SQTPeptideSpComparator.SORTBYXCORR;
        if(getsp) {
            scoreType = SQTPeptideSpComparator.SORTBYSP;
        }
         
//System.out.println("scoretype for sorting: " + scoreType); 
        SQTPeptideSpComparator comparator = new SQTPeptideSpComparator(scoreType);
        
        Collections.sort(sqts, comparator);

    }
    public double getScore(SQTPeptide s) {
        if(getsp) {
            return s.getTopHit().getSpValue();
        } else {
            return s.getTopHit().getXcorrValue();
        }
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
                 if(s.getNumMlines() > 0) {
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
