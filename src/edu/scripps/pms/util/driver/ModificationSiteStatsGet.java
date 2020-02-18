/**
 * @file ModificationSiteStatsGet.java
 * This is the source file for edu.scripps.pms.util.spectrum.ModificationSiteStatsGet
 * @author Tao Xu
 * @date $Date: 2009/01/29 19:07:02 $
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
public class ModificationSiteStatsGet {
    public static final String USAGE = "\n\n!!! USAGE: java ModificationSiteStatsGet -i input -o output -m minium_mass_shift -M maximum_mass_shift!!!";
    private static int minMassShift = 40; // negative mass shift extrem
    private static int maxMassShift = 200; // positive mass shift extrem
    private static int[][] frequency = new int[256][maxMassShift]; 
    private static int[][] negfrequency = new int[256][minMassShift]; // for negative mass shift
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
        ps.print("residues");
        for(int i = negfrequency[0].length -1; i > 0; i--) {
            ps.print("\t" + (0-i-1) );
        }
        for(int i = 0; i < frequency.length; i++) {
            ps.print("\t" + (0+i));
        }
        ps.println("\tmodfreq\tnonmodfreq");
        int totalmod = 0;
        int totalnonmod = 0;
        for(char c : residues) {
            totalmod += modfreq[c];
            totalnonmod += nonmodfreq[c];
        }
        for(char c : residues) {
            ps.print(c);
            int [] freq = negfrequency[c];
            for(int i = freq.length-1; i > 0; i--) {
                ps.print("\t" + freq[i]);
            }
            freq = frequency[c];
            for(int i = 0; i < freq.length; i++) {
                ps.print("\t" + freq[i]);
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

        String arr [] = seq.split("\\(");
        for(int i = 0; i < arr.length/2; i++) {
            char c = arr[i*2].charAt(arr[i*2].length()-1);
            float mass = Float.parseFloat(arr[i*2+1].split("\\)")[0]);
            int massindex = mass > 0 ? (int)(mass + 0.5) : (int)(Math.abs(mass - 0.5));
            if(mass > 0) {

                frequency[c][massindex]++;

            } else {
//System.out.println("mass: " + mass + "\tmassindex: " + massindex + "\tstring: " + seq);
                negfrequency[c][massindex]++;

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
        Option minOpt = new Option
            ("m", "min", true,  "min_mass_shift");
        Option maxOpt = new Option
            ("M", "max", true,  "max_mass_shift");
        //inputOpt.setRequired(true);
        //pValueOpt.setRequired(true);
        //hprdOpt.setRequired(false);
        //opts.addOption(dbNameOpt);
        opts.addOption(inputOpt);
        opts.addOption(outOpt);
        opts.addOption(minOpt);
        opts.addOption(maxOpt);

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

        String min = cli.getOptionValue("m");
        minMassShift = min == null? minMassShift : Integer.parseInt(min);
        minMassShift = Math.abs(minMassShift);
        String max = cli.getOptionValue("M");
        maxMassShift = max == null? maxMassShift : Integer.parseInt(max);


        frequency = new int[256][maxMassShift+2];
        negfrequency = new int[256][minMassShift+2];
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
