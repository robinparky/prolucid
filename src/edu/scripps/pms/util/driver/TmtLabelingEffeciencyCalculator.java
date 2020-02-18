package edu.scripps.pms.util.driver; /**
 * @file TmtLabelingEffeciencyCalculator.java
 * This is the source file for edu.scripps.pms.util.spectrum.TmtLabelingEffeciencyCalculator
 * @author Tao Xu
 * @date $Date: 2008/10/24 22:39:14 $
 */



import edu.scripps.pms.util.io.*;

//import edu.scripps.pms.util.io.DTASelectFilterReader; 
import edu.scripps.pms.util.dtaselect.Peptide; 
import edu.scripps.pms.util.dtaselect.Protein; 
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.PmsUtil; 
import edu.scripps.pms.util.stat.StatCalc; 
import java.io.*;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import java.util.HashSet;
import edu.scripps.pms.util.spectrum.*;

// command line processor
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

// this is the class for calculating the z-score for delta mass 
public class TmtLabelingEffeciencyCalculator {
    public static final String USAGE = "\n\n!!! USAGE: identifiedpeptidedb -i input -o output !!!";
    private static HashMap<String, Integer> peptide2count = new HashMap<String, Integer>(5000);
//    private static HashSet<String> peptidesInCommon = new HashSet<String>();
    private static String inputFile = null;
    private static String outputFile = null;
   
    private static int nLabeled = 0; 
    private static int nUnLabeled = 0; 
    private static int kLabeled = 0; 
    private static int kUnLabeled = 0; 
    private static int yLabeled = 0; 
    private static int yUnLabeled = 0; 
    
    public static ArrayList<Protein> getDTASelectProteins(String dtaselectFilterFile) throws IOException {

        ArrayList<Protein> result = new ArrayList(10000);
        DTASelectFilterReader reader = new DTASelectFilterReader(dtaselectFilterFile);
        for (Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {
            result.add(itr.next());
        }
        reader.close(); 
        return result;
    }
    public static void output(ArrayList<Protein> dtaselectProteins) throws IOException {
        PrintStream ps = null;
        if(outputFile == null) {
            ps = System.out;
        } else {
            ps = new PrintStream(outputFile);
        }
        
        for(Protein p : dtaselectProteins) {
            String acc = p.getAccession();
            HashSet<String> peptideadded = new HashSet<String>(100000); 
            for(Iterator<Peptide> pepItr=p.getPeptides(); pepItr.hasNext(); ) {
                Peptide peptide = pepItr.next();
                String seq = peptide.getMidSeq();
                 
                Integer count = new Integer(peptide.getRedundancy());
//System.out.println(seq + "\t" + count);
                if(!peptide2count.containsKey(seq)) {
                    peptide2count.put(seq, count);
                }
            }
        }
        for(Iterator<String> it = peptide2count.keySet().iterator(); it.hasNext();) {
            String seq = it.next();
            Integer count = peptide2count.get(seq);
//System.out.println(seq + "\t" + count);
            if(seq.startsWith("(")) {
                nLabeled += count;
            }else {
                nUnLabeled += count;
            }
            int lastindex = seq.length() - 1;
            byte [] bytes = seq.getBytes();
            for(int i = 0; i <= lastindex; i++) {
                if(bytes[i] == 'K') {
                    int nextindex = i + 1;
                    if(nextindex <= lastindex && bytes[nextindex] == '(') {
                        kLabeled += count;
                    } else {
                        kUnLabeled += count;
                    } 
                }
            }
            for(int i = 0; i <= lastindex; i++) {
                if(bytes[i] == 'Y') {
                    int nextindex = i + 1;
                    if(nextindex <= lastindex && bytes[nextindex] == '(') {
                        yLabeled += count;
                    } else {
                        yUnLabeled += count;
                    } 
                }
            }
                
        }
        System.out.println("Residue\tLabeled\tUnlabeled\tLabeling Effeciency"); 
        System.out.println("N-term\t" + nLabeled + "\t" + nUnLabeled + "\t" + nLabeled/(nLabeled+nUnLabeled+0.0)*100 + "%"); 
        System.out.println("K\t" + kLabeled + "\t" + kUnLabeled + "\t" + kLabeled/(kLabeled+kUnLabeled+0.0)*100 + "%"); 
        System.out.println("Y\t" + yLabeled + "\t" + yUnLabeled + "\t" + yLabeled/(yLabeled+yUnLabeled+0.0)*100 + "%"); 
    }
    public static void main(String args[]) throws Exception {
   

        // Command line processing
        Options opts = new Options();
        Option inputOpt = new Option
            ("i", "input", true, "Input file name - DTASelect-filer.txt run with -t 0 option");
        Option outOpt = new Option
            ("o", "output", true,  "Output file name");
        //inputOpt.setRequired(true);
        //pValueOpt.setRequired(true);
        //hprdOpt.setRequired(false);
        //opts.addOption(dbNameOpt);
        opts.addOption(inputOpt);
        opts.addOption(outOpt);

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


        try {
            //System.out.println("p value cutoff: " + minPValue);
            output(getDTASelectProteins(inputFile));

        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
        }
    }


}
