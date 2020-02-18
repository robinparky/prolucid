/**
 * @file CamsiResultEvaluator.java
 * This is the source file for CamsiResultEvaluator
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
public class CamsiResultEvaluator {
    public static final String USAGE = "\n\n!!! USAGE: camsireportgenerator -i input -o output -f falsePostiveRate -s scoreIndex -S !!!";
    private static ArrayList<HashSet<String>> peptideSets = new ArrayList<HashSet<String>>();
//    private static HashSet<String> peptidesInCommon = new HashSet<String>();
    private static final int MAXDELTAMASS = 50;    
    private static String inputFile = null;
    private static String proteinKeyFile = null;
    private static String outputFile = null;
    private static double fpr = 0.01; // false postive rate
    //private static CamsiId []identifiedscans = new CamsiId[1000000]; // identified scans
    private static ArrayList<CamsiId> ids = new ArrayList<CamsiId>(100000);
    private static HashSet<String> forwardproteins = new HashSet<String>(10000000);
    //private static HashSet<String> decopyproteins = new HashSet<String>(10000000);
    private static boolean isSorted = false;
    private static String sqtfile = null;
    private static HashSet<Integer> identifiedscans = new HashSet<Integer>(100000);
    private static HashSet<Integer> decoyscans = new HashSet<Integer>(100000);
    private static HashSet<String> identifiedpeptides = new HashSet<String>(1000000);
    private static HashSet<String> identifiedproteins = new HashSet<String>(1000000);
    private static int [] scan2freq = new int[100000];
    // number of decoyhits identified together with the foward hits
    private static HashSet<String> identifieddecoyproteins = new HashSet<String>(100000);

    private static void sortCamsiIds() {
        Collections.sort(ids);
        //isSorted = true;
        ArrayList sortedList = new ArrayList<CamsiId>(100000);
        for(int i = ids.size() -1; i > -1; i--) {
            sortedList.add(ids.get(i));
        }
        ids = sortedList;
    }
    public static void outputCamsiReport() throws IOException {
        int numforward = 0;
        int numdecoy = 0;
        int numboth = 0; // ids identified both decoy and forward
        int maxnumforward = 0;
        int maxuniqscan = 0;
        int maxnumdecoy = 0;
        int maxnumproteins = 0;
        int maxnumpeptides = 0;
        PrintStream ps = System.out;
        if(outputFile != null) {
            ps = new PrintStream(outputFile);
        }
        StringBuffer sb = new StringBuffer(1000000);
        for(int i = 0; i < ids.size(); i++) {
            CamsiId id = ids.get(i);
            if(isForwardProtein(id.getProteinIds())) {
                numforward++; 
                identifiedscans.add(new Integer(id.getScan()));
                identifiedpeptides.add(id.getPeptideSequence());
                for(Iterator<String> it = id.getProteinIds(); it.hasNext();) {
                    String protein = it.next();
                    if(isForwardProtein(protein)) {    
                        identifiedproteins.add(protein);
                    } else {
                        identifieddecoyproteins.add(protein);
                        numboth++;
//System.out.println("decoy: " + protein + "\t");
                    }
                }
            } else {
                numdecoy++;
                decoyscans.add(new Integer(id.getScan()));        
System.out.println("decoy: " + id.getProteinIds().next() + "\t" + id.getScan());
            }
            //if(numdecoy/(numforward+0.0) < fpr) {
            if(decoyscans.size()/(identifiedscans.size()+0.0) < fpr) {
                sb.append(id.getContentLine());
                sb.append("\n");
                maxnumforward = numforward;
                maxnumdecoy = decoyscans.size();
                maxuniqscan = identifiedscans.size();
                maxnumproteins = identifiedproteins.size();
                maxnumpeptides = identifiedpeptides.size();
            } else {
                ////break;
            }
        }
        int maxcount = 0;
        int scanwithmaxcount = 0;
        int numscanmorecount = 0;
        int numscan = 0;
        int totalm = 0;
        int morethan5 = 0;
        for(int i = 0; i < scan2freq.length; i++) {
            int freq = scan2freq[i];
            if(freq > 0) {
                numscan++;
                totalm += freq;
                if(freq > 1) {
                    numscanmorecount++;
                    if(freq > 4) {
                        morethan5++;
                    }
                    if(freq > maxcount) {
                        maxcount = freq;
                        scanwithmaxcount = i;
                    }
                }
            }
        }
        ps.println("TotalNumForward: " + numforward + "\tNumDecoy: " + numdecoy);
        ps.println("NumRowsPassed " + fpr + " filter: " + maxnumforward);
        ps.println("maxNumDecoyScan: " + maxnumdecoy + "\tmaxuniqscans: " + maxuniqscan);
        ps.println("NumUniqueForforwardPeptideSequences: " + maxnumpeptides);
        ps.println("NumUniqueForforwardProteinHits: " + maxnumproteins);
        ps.println("TotalNumUniqueScansForforwardHits: " + identifiedscans.size());
        //ps.println("NumUniqueForforwardProteinHits: " + identifiedproteins.size());
        ps.println("TotalNumUniqueDecoyProteindTogetherWithForforwardProteinHits: " + identifieddecoyproteins.size() + "\tNumBoth: " + numboth);
        //ps.println("NumUniqueForforwardPeptideSequences: " + identifiedpeptides.size());
        ps.println("SizeOfForwardproteins: " + forwardproteins.size());
        ps.println("scan with most ids: " + scanwithmaxcount + " " + maxcount + "\tnumscan: " + numscan + "\tnumscanwithmorethan1id: " + numscanmorecount + "\taveragenumidperscan: " + totalm/(numscan+0.0) + "\tscanwithmorethan5ids: " + morethan5); 
        ps.println(sb); 
    }
    public static double getLowestScore(ArrayList<CamsiId> ids) {
        if(!isSorted) {
            sortCamsiIds();
        }
        return ids.get(0).getScore();
    }

    public static boolean isForwardProtein(String ac) {
        if(forwardproteins.contains(ac)) {
            return true;
        } 
        return false;
    }
    public static boolean isForwardProtein(Iterator<String> proteins) {
        boolean isForward = false;
        while(proteins.hasNext()) {
             
            String id = proteins.next();    
            if(id.startsWith("REVERSED")) {
                return false;
            }
            //isForward = isForward && forwardproteins.contains(id);
            if(forwardproteins.contains(id)) {
//System.out.println("isforward: " + id);
                isForward = true;
            } else {
//System.out.println("notforward: " + id);
            }
        }
        return isForward;
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
    public static void readProteinKeys(String fileName) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(fileName));
        String line = reader.readLine();
        
//System.out.println("Line:\n" + line);
        while(line != null) {
            if("".equals(line)) {
                continue;
            }
            String [] arr = line.split("\\t");
            if(line != null) {
                if(arr[1].startsWith("Decoy")) {
                    //decoyproteins.add(arr[0]);    

//System.out.println("Decoy\n" + line);
                } else {
//System.out.println("not decoy\n" + line);
                    forwardproteins.add(arr[0]);
                }
            }
            line = reader.readLine();    
        }
        reader.close();

    }
    public static void readCamsiResult(String fileName, int scoreIndex) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(fileName));
        String line = reader.readLine();
         
        while(line != null) {
            if("".equals(line)) {
                line = reader.readLine();
                continue;
                
            }
            ////String [d] arr = line.split("\\s");
            //System.out.println("in readcamsiresult, line: " + line);
            if(line != null) {
                CamsiId id = new CamsiId(line, scoreIndex);
                scan2freq[id.getScan()]++;    
                ids.add(id);
            }
            line = reader.readLine();    
        }
        reader.close();
        if(!isSorted) {
            sortCamsiIds();
        } 
    }
    
    
    public static void main(String args[]) throws Exception {

        // Command line processing
        Options opts = new Options();
        Option inputOpt = new Option
            ("i", "input", true, "Input file name"); //submitted camsi result
        Option proteinKeyOpt = new Option
            ("p", "proteinkeyfile", true, "Protein Key file");
        Option pvalueOpt = new Option
            ("f", "falsepositive", true,  "false positive rate");
        Option outOpt = new Option
            ("o", "output", true,  "Output file name");
        Option scoreOpt = new Option
            ("s", "score", true,  "score index");

        scoreOpt.setRequired(false);
        pvalueOpt.setRequired(true);
        proteinKeyOpt.setRequired(true);
        inputOpt.setRequired(true);
        outOpt.setRequired(false);
        //hprdOpt.setRequired(false);
        opts.addOption(pvalueOpt);
        opts.addOption(proteinKeyOpt);
        opts.addOption(inputOpt);
        opts.addOption(outOpt);
        opts.addOption(scoreOpt);

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
        proteinKeyFile = cli.getOptionValue("p");
        outputFile = cli.getOptionValue("o");
        String sindex = cli.getOptionValue("s");
        int scoreIndex = sindex == null? 3 : Integer.parseInt(sindex);
        isSorted = sindex == null? false : true;
        fpr = Double.parseDouble(cli.getOptionValue("f"));
        //boolean digprotein = digdeeper && cli.getOptionValue("p") != null;
        try {
            //System.out.println("p value cutoff: " + minPValue);
            //output(readCamsiResult(inputFile));
            System.out.println("start reading results");
            readCamsiResult(inputFile, scoreIndex);
            //readSqts(sqtfile);
            System.out.println("finished reading results");
            readProteinKeys(proteinKeyFile);
            System.out.println("finished reading proteinkeys");
            outputCamsiReport();    
        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
        }
    }

}
