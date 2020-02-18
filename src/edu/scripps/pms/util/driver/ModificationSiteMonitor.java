/**
 * @file ModificationSiteMonitor.java
 * This is the source file for edu.scripps.pms.util.spectrum.ModificationSiteMonitor
 * @author Tao Xu
 * @date $Date: 2008/10/24 22:39:14 $
 */



import edu.scripps.pms.util.io.*;

//import edu.scripps.pms.util.io.DTASelectFilterReader; 
import edu.scripps.pms.util.dtaselect.Peptide; 
import edu.scripps.pms.util.dtaselect.Protein; 
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
public class ModificationSiteMonitor {
    public static final String USAGE = "\n\n!!! USAGE: java ModificationSiteMonitor -i input -o output !!!";
    private static int[][] frequency = new int[11][256]; 
    private static int[] modfreq = new int[256]; 
    private static int[] nonmodfreq = new int[256]; 
//    private static HashSet<String> peptidesInCommon = new HashSet<String>();
    private static String inputFile = null;
    private static boolean modifiedOnly = true;
    private static String outputFile = null;
    
        
    public static void output() throws IOException {
        PrintStream ps = null;
        if(outputFile == null) {
            ps = System.out;
        } else {
            ps = new PrintStream(outputFile);
        }
        char [] residues = {'G','A','S','P','V','T','C','L','I','N','D','Q','K','E','M','H','F','R','Y','W' };
        ps.println("residue\t0\t-5\t-4\t-3\t-2\t-1\t1\t2\t3\t4\t5\tmodfreq\tnonmodfreq");
        int totalmod = 0;
        int totalnonmod = 0;
        for(char c : residues) {
            totalmod += modfreq[c];
            totalnonmod += nonmodfreq[c];
        }
        for(char c : residues) {
            ps.print(c + "\t");

            for(int i = 0; i < 11; i++) {
                ps.print(frequency[i][c] + "\t"); 
            }
            ps.print(modfreq[c]/(totalmod+0.0) + "\t"); 
            ps.print(nonmodfreq[c]/(totalnonmod+0.0) + "\t"); 
            ps.println();
        }
    }    
    public static void calcFrequency(String dtaselectFilterFile) throws IOException {
        
        DTASelectFilterReader reader = new DTASelectFilterReader(dtaselectFilterFile);
        for(Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {
            Protein p = itr.next(); 
            String accession = p.getAccession();
            String description = p.getDescription(); 
            for(Iterator<Peptide> pepItr=p.getPeptides(); pepItr.hasNext(); ) {
                Peptide peptide = pepItr.next();
                if(isModifiedPeptide(peptide)) {
                    getModificationSiteInfo(peptide);
                    calcResidueFrequency(peptide, modfreq); 
                } else {
                    calcResidueFrequency(peptide, nonmodfreq); 
                }
            
            }
        }
    }
    public static void calcResidueFrequency(Peptide p, int [] freq) {
        String seq = p.getSequence();
        int maxindex = seq.length() - 2;
        for(int i = 2; i < maxindex; i++) {
            freq[seq.charAt(i)]++;
        }
    }
    public static void getModificationSiteInfo(Peptide p) {
        String seq = p.getSequence();
        int msite = seq.indexOf("~") - 1;
//System.out.println(seq + "\t" + msite);
        if(msite < 0) {
            return;
        }
        char c = seq.charAt(msite);
        frequency[0][c]++;
        for(int i = 1; i < 6; i++) {
            int left = msite - i; 
            int right = msite + i + 1; 
            if(left >= 0) {
                char cleft = seq.charAt(left);
                frequency[i][cleft]++;
            }
            if(right < seq.length()) {
                char cright = seq.charAt(right);
                frequency[5+i][cright]++;
            }
           
        } 
    }
    public static boolean isModifiedPeptide(Peptide p) {
        String seq = p.getSequence();
        boolean isModified = seq.indexOf("(") != -1 || seq.indexOf("~") != -1 || seq.indexOf("*") != -1 || seq.indexOf("#") != -1 || seq.indexOf("@") != -1;
            //System.out.println(seq + "\t" + isModified);
        //System.out.println(seq + "\t" + isModified);
        return isModified;
    }
    public static void main(String args[]) throws Exception {
   

        // Command line processing
        Options opts = new Options();
        Option inputOpt = new Option
            ("i", "input", true, "Input file name - DTASelect-filer.txt run with --DM option");
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
            calcFrequency(inputFile);
            output();

        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
        }
    }
}
