package edu.scripps.pms.util.driver;

/**
 * @file AdvancedDeltaMassFilter.java
 * This is the source file for edu.scripps.pms.util.spectrum.AdvancedDeltaMassFilter
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
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.mspid.ProteinDatabase;
import edu.scripps.pms.mspid.SearchParams;
import edu.scripps.pms.mspid.MassCalculator;

// command line processor
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

// this is the class for calculating the z-score for delta mass 
public class AdvancedDeltaMassFilter {
    public static final String USAGE = "\n\n!!! USAGE: AdvancedDeltaMassFilter -i input -o output -p pValueCutOff !!!";
    private static ArrayList<HashSet<String>> peptideSets = new ArrayList<HashSet<String>>();
//    private static HashSet<String> peptidesInCommon = new HashSet<String>();
    private static final int MAXDELTAMASS = 50;    
    private static String inputFile = null;
    private static double minPValue = 0.05;
    private static boolean modifiedOnly = true;
    private static String outputFile = null;
    private static double [] avgAaMasses = new double[256];
    private static double [] monoAaMasses = new double[256];
    private static String databaseName; 
    private static ProteinDatabase pdb;
    private static MassCalculator masscalc; 
 
    public static ArrayList<Protein> getDTASelectProteins(String dtaselectFilterFile) throws IOException {

        ArrayList<Protein> result = new ArrayList(10000);
        DTASelectFilterReader reader = new DTASelectFilterReader(dtaselectFilterFile);
        pdb = new ProteinDatabase(reader.getDbFilePathAndName());
        //databaseName = reader.getDbFileName();
        //System.out.println("databaseName: " + databaseName);
        for (Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {
            result.add(itr.next());
        }
        reader.close(); 
        return result;
    }
    public static boolean isAcceptable(String seq, Fasta f) {
        if(seq.indexOf("(") == -1) {
            return true;
        } 
        
//System.out.println( seq );
        return checkNTerm(seq, f) && checkCTerm(seq, f); 

    }    
    private static boolean checkNTerm(String seq, Fasta f) {
        boolean result = true;
        int index = seq.indexOf("(");
        if(index > 5 || index < 0) {
            return true;
        }
        // there is a mass shift located at the first, second or third N-term residue
        double massShift = Double.parseDouble(seq.split("\\(")[1].split("\\)")[0]);
        char n = seq.charAt(0); 
        double ntermResidueMass = monoAaMasses[n];

        double lowLimit = ntermResidueMass - 0.2;
        double highLimit = ntermResidueMass + 0.2;
        result &= !(massShift >= lowLimit && massShift <= highLimit);
//System.out.println(result + "\t" + seq + "\t" + massShift + "\tntermMassShift: " + ntermResidueMass + "\tlowLimit: " + lowLimit + "\thighlimit: " + highLimit);
        return result; 
        
    }
    private static boolean checkCTerm(String seq, Fasta f) {
        boolean result = true;
        int closeIndex = seq.lastIndexOf(")");
        
        if(seq.length() - closeIndex > 5 || closeIndex < 0) {
            return true;
        }
        int openIndex = seq.lastIndexOf("(") + 1;
        String shift = seq.substring(openIndex, closeIndex);
       //System.out.println("mass shiflt: " + shift + "\t" + seq); 
        
        // there is a mass shift located at the first, second or third C-term residue
        double massShift = Double.parseDouble(shift);

        String proteinseq = f.getSequence();
        int pos = proteinseq.indexOf(seq);
        int lastindex = pos - seq.length();

        char c = seq.charAt(seq.length()-1); 
         
        double ctermResidueMass = monoAaMasses[c];

        double lowLimit = ctermResidueMass - 0.02;
        double highLimit = ctermResidueMass + 0.02;
        result &= !(massShift >= lowLimit && massShift <= highLimit);
        return result; 
        
    }
    public static void output(ArrayList<Protein> dtaselectProteins) throws IOException {
        int length = MAXDELTAMASS*2+1;
        int [] normZscoreFreq = new int[length];
        int [] modZscoreFreq = new int[length];
        int [] normDeltaMassFreq = new int[length];
        int [] modDeltaMassFreq = new int[length];
        PrintStream ps = null;
        if(outputFile == null) {
            ps = System.out;
        } else {
            ps = new PrintStream(outputFile);
        }
        StringBuffer sb = new StringBuffer(1000000);
        StatCalc calc = new StatCalc(); 
        StatCalc mCalc = new StatCalc();
        for(Protein p : dtaselectProteins) {
            for(Iterator<Peptide> pepItr=p.getPeptides(); pepItr.hasNext(); ) {
                Peptide peptide = pepItr.next();
                float deltamass = peptide.getDeltaMass();
//System.out.println(deltamass);
                if(Math.abs(deltamass) < 10) {
                    if(!isModifiedPeptide(peptide)) {
                        calc.enter(deltamass);
                    } else {
                        mCalc.enter(deltamass);
                    }
                }
            }
        }
        double mean = calc.getMean();
        double sdv = calc.getStandardDeviation();
        String sumLine = mean + "\t" + sdv + "\t" + calc.getCount() + "\t";
        sumLine = sumLine + mCalc.getMean() + "\t" + mCalc.getStandardDeviation() + "\t" + mCalc.getCount();
        ps.println(sumLine);
        
        for(Protein p : dtaselectProteins) {
            String accession = p.getAccession();
            String description = p.getDescription(); 
            Fasta f = pdb.accession2Fasta(accession);

            for(Iterator<Peptide> pepItr=p.getPeptides(); pepItr.hasNext(); ) {
                Peptide peptide = pepItr.next();
                float deltamass = peptide.getDeltaMass();
                double zscore = Math.abs((deltamass - mean)/sdv);
                double pValue = StatCalc.zScore2PValue(zscore);
                 
                if(isModifiedPeptide(peptide) && pValue >= minPValue && peptide.getDeltCnValue() > 0.01) {
                    if(isAcceptable(peptide.getSequence(), f)) {
                        //sb.append(accession + "\t" + description + "\t");
                        //sb.append(peptide.getInfo() + "\t" + StatCalc.zScore2PValue(zscore) + "\n");
                        sb.append(accession + "\t" + peptide.getInfo() + "\t");
                        sb.append(StatCalc.zScore2PValue(zscore) + "\t" + description + "\n");
                    }
                }
                int [] deltaMassFreq = isModifiedPeptide(peptide)? modDeltaMassFreq : normDeltaMassFreq;
                deltaMassFreq[getDeltaMassIndex(deltamass)]++;
            
            }
        }
        ps.println("delta mass p value cutoff: " + minPValue);
        ps.println("norm\n" + getDeltaMassDistribution(normDeltaMassFreq));
        ps.println("mod\n" + getDeltaMassDistribution(modDeltaMassFreq));
        ps.println("\nAccession\t" + Peptide.getInfoHeader() + "\tp-value\tDescription");
        ps.println(sb.toString()); 
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

    public double getNtermDeltaMassInPpm(Fasta prot, String peptide, double peptidemass, double massshift) {
        double diff = 0;
        double mass = 0;
        String seq = prot.getSequence();
        int index = seq.indexOf(peptide);
        while(index > -1) {
            mass += masscalc.getPrecursorMass(seq.charAt(index--));
        }
        return diff/peptidemass*1000000;
    }


    public static void main(String args[]) throws Exception {
        masscalc = new MassCalculator(new SearchParams("search.xml")); 
        // Command line processing
        Options opts = new Options();
        Option inputOpt = new Option
            ("i", "input", true, "Input file name - DTASelect-filer.txt run with --DM option");
        Option pValueOpt = new Option
            ("p", "P-ValueCutOff", true,  "P-Value CutOff");
        Option outOpt = new Option
            ("o", "output", true,  "Output file name");
        //inputOpt.setRequired(true);
        //pValueOpt.setRequired(true);
        //hprdOpt.setRequired(false);
        //opts.addOption(dbNameOpt);
        opts.addOption(inputOpt);
        opts.addOption(pValueOpt);
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

        if(cli.getOptionValue("p") != null) {
            minPValue = Double.parseDouble(cli.getOptionValue("p"));        
        }

        try {
            //System.out.println("p value cutoff: " + minPValue);
            output(getDTASelectProteins(inputFile));

        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
        }
    }

    // initiatet the values of average and mono aa masses
    static {
        avgAaMasses['G'] =  57.05192f;   monoAaMasses['G'] =  57.0214636f;
        avgAaMasses['A'] =  71.07880f;   monoAaMasses['A'] =  71.0371136f;
        avgAaMasses['S'] =  87.07820f;   monoAaMasses['S'] =  87.0320282f;
        avgAaMasses['P'] =  97.11668f;   monoAaMasses['P'] =  97.0527636f;
        avgAaMasses['V'] =  99.13256f;   monoAaMasses['V'] =  99.0684136f;
        avgAaMasses['T'] = 101.10508f;   monoAaMasses['T'] = 101.0476782f;
        avgAaMasses['C'] = 103.13880f;   monoAaMasses['C'] = 103.0091854f;
        avgAaMasses['L'] = 113.15944f;   monoAaMasses['L'] = 113.0840636f;
        avgAaMasses['I'] = 113.15944f;   monoAaMasses['I'] = 113.0840636f;
        avgAaMasses['X'] = 113.15944f;   monoAaMasses['X'] = 113.0840636f;
        avgAaMasses['N'] = 114.10384f;   monoAaMasses['N'] = 114.0429272f;
        avgAaMasses['O'] = 114.14720f;   monoAaMasses['O'] = 114.0793126f;
        avgAaMasses['B'] = 114.59622f;   monoAaMasses['B'] = 114.5349350f;
        avgAaMasses['D'] = 115.08860f;   monoAaMasses['D'] = 115.0269428f;
        avgAaMasses['Q'] = 128.13072f;   monoAaMasses['Q'] = 128.0585772f;
        avgAaMasses['K'] = 128.17408f;   monoAaMasses['K'] = 128.0949626f;
        avgAaMasses['Z'] = 128.62310f;   monoAaMasses['Z'] = 128.5505850f;
        avgAaMasses['E'] = 129.11548f;   monoAaMasses['E'] = 129.0425928f;
        avgAaMasses['M'] = 131.19256f;   monoAaMasses['M'] = 131.0404854f;
        avgAaMasses['H'] = 137.14108f;   monoAaMasses['H'] = 137.0589116f;
        avgAaMasses['F'] = 147.17656f;   monoAaMasses['F'] = 147.0684136f;
        avgAaMasses['R'] = 156.18748f;   monoAaMasses['R'] = 156.1011106f;
        avgAaMasses['Y'] = 163.17596f;   monoAaMasses['Y'] = 163.0633282f;
        avgAaMasses['W'] = 186.21320f;   monoAaMasses['W'] = 186.0793126f;
    }


}
