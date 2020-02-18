/**
 * @file MethionineCounter.java
 * This is the source file for edu.scripps.pms.util.spectrum.MethionineCounter
 * @author Tao Xu
 * @date $Date: 2009/01/29 19:07:02 $
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
public class MethionineCounter {
    public static final String USAGE = "\n\n!!! USAGE: identifiedpeptidedb -i input -o output !!!";
    private static ArrayList<HashSet<String>> peptideSets = new ArrayList<HashSet<String>>();
//    private static HashSet<String> peptidesInCommon = new HashSet<String>();
    private static String inputFile = null;
    private static String outputFile = null;
    
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
       
        int numProteins = 0;
        int numMContainingProteins = 0; // protein containing M
        int numOMContainingProteins = 0; // protein containing oxidated M
        
        HashSet<String> idedPeptides = new HashSet<String>(100000); 
        ps.println("Accession\tPeptideCount\tMPeptideCount\tOMPeptideCount\tDescription");
        for(Protein p : dtaselectProteins) {
            int peptidecount = 0;
            int mpeptidecount = 0;
            int ompeptidecount = 0; // oxidated m containing peptide count
             
            numProteins++;
            String acc = p.getAccession();
            HashSet<String> peptideadded = new HashSet<String>(100000); 
            boolean isMProtein = false;
            boolean isOMProtein = false;
            for(Iterator<Peptide> pepItr=p.getPeptides(); pepItr.hasNext(); ) {
                  
                Peptide peptide = pepItr.next();
                String seq = peptide.getMidSeq();
                idedPeptides.add(seq);
                peptidecount++;
                if(seq.indexOf("M") != -1) {
                    mpeptidecount++;
                    isMProtein = true;
                    if(seq.indexOf("M(") != -1) {
                        isOMProtein = true;
                        ompeptidecount++;
                    }
                } 
                
            }
            numMContainingProteins += isMProtein? 1 : 0;
            numOMContainingProteins += isOMProtein? 1 : 0;
            ps.println(acc + "\t" + peptidecount + "\t" + mpeptidecount + "\t" + ompeptidecount + "\t" + p.getDescription());
        }
        int numPeptides = idedPeptides.size();
        int numMContainingPeptides = 0;
        int numOMContainingPeptides = 0;
        for(Iterator<String> it = idedPeptides.iterator(); it.hasNext();) {
            String seq = it.next();
            if(seq.indexOf("M") != -1) {
                numMContainingPeptides++;
                if(seq.indexOf("M(") != -1) {
                    numOMContainingPeptides++;
                }
           }
        }
        ps.println("Number of identified proteins: \t" + numProteins); 
        ps.println("Number of identified proteins with identified M containg peptides: \t" + numMContainingProteins); 
        ps.println("Number of identified proteins with identified oxidated M containg peptides: \t" + numOMContainingProteins); 
        ps.println("Number of identified peptides: \t" + numPeptides); 
        ps.println("Number of identified M containg peptides: \t" + numMContainingPeptides); 
        ps.println("Number of identified oxidated M containg peptides: \t" + numOMContainingPeptides); 
    }
    public static void main(String args[]) throws Exception {
   

        // Command line processing
        Options opts = new Options();
        Option inputOpt = new Option
            ("i", "input", true, "Input file name - DTASelect-filer.txt run w/o -DM option");
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
