/**
 * @file CamsiReportGenerator.java
 * This is the source file for CamsiReportGenerator
 * @author Tao Xu
 * @date $Date
 */



import edu.scripps.pms.util.io.*;

//import edu.scripps.pms.util.io.DTASelectFilterReader; 
import edu.scripps.pms.util.dtaselect.Peptide; 
import edu.scripps.pms.util.dtaselect.Protein; 
import edu.scripps.pms.util.PmsUtil; 
import edu.scripps.pms.util.stat.StatCalc; 
import edu.scripps.pms.util.sqt.SQTPeptide;
import edu.scripps.pms.util.sqt.MLine;


import java.io.*;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import edu.scripps.pms.util.spectrum.*;

// command line processor
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

// this is the class for generate CAMSI repotring files. Modified from DeltaMassZScoreCalculator.java file
public class CamsiReportGenerator {
    public static final String USAGE = "\n\n!!! USAGE: camsireportgenerator -i input -o output [-p] [-d] !!!" 
                                       + "\nThe DTASelect should be run with -t 0 option!" + 
                                          "\n use -d option to retrieve low score hits from DTASelect.txt file" +
                                          "\n use -p together with -d to get all low score peptide ids from identified proteins";
    private static ArrayList<HashSet<String>> peptideSets = new ArrayList<HashSet<String>>();
//    private static HashSet<String> peptidesInCommon = new HashSet<String>();
    private static final int MAXDELTAMASS = 50;    
    private static String inputFile = null;
    private static String outputFile = null;
    private static double minPValue = 0.01;
    private static CamsiId []identifiedscans = new CamsiId[1000000]; // identified scans
    private static HashMap<String, ArrayList<CamsiId>> seq2CamsiIds = new HashMap<String, ArrayList<CamsiId>>(1000000);
    private static boolean isSorted = true;
    private static String sqtfile = null;

    private static void sortCamsiIds() {
        Set<String> keys = seq2CamsiIds.keySet();
        Iterator<String> peptides = keys.iterator();
        while(peptides.hasNext()) {
            ArrayList<CamsiId> id = seq2CamsiIds.get(peptides.next());
            if(id != null && id.size() > 1) {
                Collections.sort(id);
            }
        }
        isSorted = true;
    }
    public static void outputCamsiReport() {
        for(int i = 0; i < identifiedscans.length; i++) {
            if(identifiedscans[i] != null) {
                System.out.println(identifiedscans[i].output());
            } 
        }
        /*
        Set<String> seqs = seq2CamsiIds.keySet();
        Iterator<String> it = seqs.iterator();
        while(it.hasNext()) {
            //ArrayList<CamsiId> ids = seq2CamsiIds.get(it.next());
            Iterator<CamsiId> ids = seq2CamsiIds.get(it.next()).iterator();
            while(ids.hasNext()) {
                System.out.println(ids.next().output());
            }
        } 
        */
    }
    public static double getLowestScore(ArrayList<CamsiId> ids) {
        if(!isSorted) {
            sortCamsiIds();
        }
        return ids.get(0).getScore();
    }

    public static boolean isModifiedPeptide(String seq) {
        return seq.indexOf("(") != -1;
    }
    public static String removeModInfo(String seq) {
        if(seq.indexOf("(") == -1) {
            return seq;
        }
        String arr [] = seq.split("\\(");
        String [] arr1 = arr[1].split("\\)");
        return arr[0] + arr1[1];
    }
    public static void getDTASelectProteins(String dtaselectFilterFile) throws IOException {

        DTASelectFilterReader reader = new DTASelectFilterReader(dtaselectFilterFile);
        StatCalc calc = new StatCalc(); 
        for (Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {
            Protein p = itr.next();
            String proteinid = p.getAccession();
            if(proteinid.startsWith("Reverse")) {
                continue;
            }
            for(Iterator<Peptide> pepItr=p.getPeptides(); pepItr.hasNext(); ) {
                Peptide peptide = pepItr.next();
                float deltamass = peptide.getDeltaMass();
//System.out.println(deltamass);
                if(Math.abs(deltamass) < 50) {
                    calc.enter(deltamass);
                }
                int scan = peptide.getScanNumber(); 
                 if(sqtfile == null) {
                     sqtfile = peptide.getFileName() + ".sqt";
                 }
                CamsiId camsiid = null;
                if(identifiedscans[scan] != null) { 
                    camsiid = identifiedscans[scan];
                } else {
                    
                    camsiid = new CamsiId(scan);
                    identifiedscans[scan] = camsiid;
                    String seq = peptide.getSequence();
                    double confidence = peptide.getConfValue()/100.0;
//System.out.println("Confidence: " + confidence);
                    double zscore = peptide.getSpScoreValue();
                    camsiid.setScore(zscore);
                    camsiid.setFalsePositiveRate(1-confidence);
                
                    if(isModifiedPeptide(seq)) {
                        camsiid.isModified(true);
                        //camsiid.setPeptideSequence(removeModInfo(seq));
                        camsiid.setPeptideSequence(seq);
                    } else {
                        camsiid.setPeptideSequence(seq);
                    }
                }
                camsiid.addProteinId(proteinid);


                //String seq = peptide.getSequence() + peptide.getChargeState();
               /*
                ArrayList<CamsiId> ids = seq2CamsiIds.get(seq);

                if(ids == null) {
                    ids = new ArrayList<CamsiId>(5);
                    seq2CamsiIds.put(seq, ids);
                    ids.add(camsiid);
                    identifiedscans[scan] = camsiid; 
                } else {
                    if(identifiedscans[scan] == null) {
                        ids.add(camsiid);
                        identifiedscans[scan] = camsiid;
                    } 
                }            
                 
                */
            }
        }
        reader.close(); 
    }
    
    public static void readSqts(String sqtfile) throws IOException {
        SqtStatistics reader = new SqtStatistics(sqtfile);
         SQTPeptide peptide;
         MLine mLine;
         String temp;

         int numLowIdAdded = 0;
         for(Iterator<SQTPeptide> itr = reader.getSQTPeptide(); itr.hasNext();) {
             peptide = itr.next();
             int scan = peptide.getLoScanNumber();
             String charge = peptide.getChargeState();
             //if(identifiedscans[scan] != null) { // seqeunce passed dtaselct
                 MLine tophit = peptide.getTopHit();
                 numLowIdAdded += checkId(scan, charge, tophit);

                 //MLine topSpRankHit = peptide.getTopSpRankHit();
               //  numLowIdAdded += checkId(scan, charge, topSpRankHit);
              
                
             //}
        //     double prcMass = Double.parseDouble(peptide.getCalMZ().trim());


         }
         System.out.println("Number of low score hits added: " + numLowIdAdded);

    }
    private static int checkId(int scan, String charge, MLine m) {
        if(m == null) {
            return 0;
        }
        int numAdded = 1;
        if(identifiedscans[scan] != null) {
            
            numAdded = 0;
        }
        String seq = m.getSequence() + charge;
       
        CamsiId hit = null;

        ArrayList<CamsiId> hits = seq2CamsiIds.get(seq);
        if(hits != null) { // identified
            double score = getLowestScore(hits);  // dtaselect confidence score 
            if(identifiedscans[scan] == null) {
                hit = new CamsiId(scan);
                identifiedscans[scan] = hit; 
                hits.add(hit);
                boolean isModified = isModifiedPeptide(seq);
                if(isModifiedPeptide(seq)) { 
                    seq = removeModInfo(seq); 
                    hit.isModified(isModified);
                
                }
                hit.setPeptideSequence(seq);
            } else {
                hit = identifiedscans[scan];
            }

            hit.setFalsePositiveRate(1-score);
            score = m.getSpValue();
            hit.setScore(score);
            for(Iterator<String> accs = m.getLLine(); accs.hasNext();) {
                hit.addProteinId(accs.next());

            }
            return numAdded;
        }
        return 0;
    }
    public static boolean isAcceptable(String seq) {
        if(seq.indexOf("(") == -1) {
            return true;
        } else {

//System.out.println( seq );
        }
        double massShift = Double.parseDouble(seq.split("\\(")[1].split("\\)")[0]);
        char n = seq.charAt(0); 
        char c = seq.charAt(seq.length()-1); 
//System.out.println( seq + "\t" + massShift + "\tntermMassShift: " + ntermResidueMass );
//System.out.println(result + "\t" + seq + "\t" + massShift + "\tntermMassShift: " + ntermResidueMass + "\tlowLimit: " + lowLimit + "\thighlimit: " + highLimit);
        return true; 

    }    
    
    private static String getDeltaMassDistribution(int [] freq) {
        StringBuffer ppmLine = new StringBuffer(1000);
        StringBuffer freqLine = new StringBuffer(1000);
        StringBuffer relativeFreqLine = new StringBuffer(10000);
        int length = freq.length;
        int total = 0;
        for(int f : freq) {
            total += f;
        } 
        for(int i = 0; i < length; i++) {
            ppmLine.append(i-MAXDELTAMASS);
            ppmLine.append("\t");
            freqLine.append(freq[i]);
            freqLine.append("\t");
            relativeFreqLine.append(freq[i]/(0.0 + total) + "\t");
        }
        return "total\t" + total + "\n" + ppmLine + "\n" + freqLine + "\n" + relativeFreqLine; 
    }
    private static int getDeltaMassIndex(float deltaMass) {
        int index = 0;
        if(deltaMass <= -MAXDELTAMASS) {
            return 0;
        } 
        if(deltaMass >= MAXDELTAMASS) {
            return MAXDELTAMASS*2;
        }
              
        return (int) (deltaMass + MAXDELTAMASS + 0.5);
    }
    public static boolean isModifiedPeptide(Peptide p) {
        String seq = p.getSequence();
        boolean isModified = seq.indexOf("(") != -1 || seq.indexOf("*") != -1 || seq.indexOf("#") != -1 || seq.indexOf("@") != -1;
            //System.out.println(seq + "\t" + isModified);
        //System.out.println(seq + "\t" + isModified);
        return isModified;
    }
    public static void main(String args[]) throws Exception {
   

        // Command line processing
        Options opts = new Options();
        Option inputOpt = new Option
            ("i", "input", true, "Input file name - DTASelect-filer.txt run with --DM option");
        Option digpOpt = new Option
            ("p", "dig_protein", false,  "Dig all peptide of identified protein");
        Option outOpt = new Option
            ("o", "output", true,  "Output file name");
        Option digOpt = new Option
            ("d", "dig", false,  "dig deeper flag");
        //inputOpt.setRequired(true);
        //digpOpt.setRequired(true);
        //hprdOpt.setRequired(false);
        //opts.addOption(dbNameOpt);
        opts.addOption(inputOpt);
        opts.addOption(digpOpt);
        opts.addOption(outOpt);
        opts.addOption(digOpt);

        BasicParser cliParser = new BasicParser();
        CommandLine cli = null;

        try {
            cli = cliParser.parse(opts, args);
            
        } catch (Exception e) {
            HelpFormatter helpFormatter = new HelpFormatter();
            helpFormatter.printHelp(USAGE, opts);
            e.printStackTrace();
            System.out.println
                ("\nPlease use \"java -Xmx1500M ... \" when processing a large dataset.");
            System.exit(1);
        }
        inputFile = cli.getOptionValue("i");
        if(inputFile == null) {
            inputFile = "DTASelect-filter.txt";
        }
        outputFile = cli.getOptionValue("o");

        boolean digdeeper = cli.getOptionValue("d") != null;
        boolean digprotein = digdeeper && cli.getOptionValue("p") != null;
        try {
            //System.out.println("p value cutoff: " + minPValue);
            //output(getDTASelectProteins(inputFile));
            getDTASelectProteins(inputFile);
            //readSqts(sqtfile);
            outputCamsiReport();    

        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
        }
    }


}
