/**
 * @file IdentifiedMs2DeltaMassGetter.java
 * This is the source file for edu.scripps.pms.util.spectrum.IdentifiedMs2DeltaMassGetter
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
public class IdentifiedMs2DeltaMassGetter {
    //public static final String USAGE = "\n\n!!! USAGE: identifiedMs2DeltaMassGetter -i input -o output -p for ppm !!!";
    public static final String USAGE = "\n\n!!! USAGE: identifiedMs2DeltaMassGetter -i input -o output !!!";
    private static ArrayList<HashSet<String>> peptideSets = new ArrayList<HashSet<String>>();
//    private static HashSet<String> peptidesInCommon = new HashSet<String>();
    private static final int MAXDELTAMASS = 50;    
    private static String inputFile = null;
    private static double minPValue = 0.05;
    private static boolean modifiedOnly = true;
    private static String outputFile = null;
    private static double [] avgAaMasses = new double[256];
    private static double [] monoAaMasses = new double[256];
    private static boolean outputppm = false;

    // use ms2 scannum as index to find the precursor scan number
    private static HashMap<String, int[]> file2precursorscans = new HashMap<String, int[]>(10000); 
    private static HashMap<String, StringBuffer[]> file2output = new HashMap<String, StringBuffer[]>(10000);
    private static final int MAXSCANNUM = 500000;    
    
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
        //StringBuffer sb = new StringBuffer(1000000);
        //sumLine = sumLine + mCalc.getMean() + "\t" + mCalc.getStandardDeviation() + "\t" + mCalc.getCount();
        //ps.println(sumLine);
        HashSet<Peptide> peptides = new HashSet<Peptide>(1000000);
 
        for(Protein p : dtaselectProteins) {
            //String accession = p.getAccession();
            //String description = p.getDescription(); 
            for(Iterator<Peptide> pepItr=p.getPeptides(); pepItr.hasNext(); ) {
                Peptide peptide = pepItr.next();
                if(!peptides.contains(peptide)) { 
                    peptides.add(peptide); 
                }

            }
        }
            //System.out.println("number of ids: " + peptides.size());
       
        for(Iterator<Peptide> it = peptides.iterator(); it.hasNext();) {
            Peptide peptide = it.next();
            float deltamass = peptide.getDeltaMass();
            String mplush = peptide.getMhPlus();
            String calcmplush = peptide.getCalcMHplus();
            String chargestate = peptide.getChargeState();
            double absdeltamass = Double.parseDouble(mplush) - Double.parseDouble(calcmplush);
            //double deltamass = peptide.getDeltaMass();
            String ms2file = peptide.getFileName();
            int ms2scannum = peptide.getScanNumber();
            int precursorscannum = getPrecursorScanNumber(ms2file, ms2scannum);
            StringBuffer sb = getOutput(ms2file, precursorscannum); 
            //sb.append("\t" + calcmplush + "\t" + mplush + "\t" + chargestate); 
            sb.append("\t" + calcmplush + "\t" + mplush + "\t" + deltamass + "\t" + chargestate); 
            
        } 
        for(Iterator<String> it = file2output.keySet().iterator(); it.hasNext();) {
            String file = it.next();
            StringBuffer [] sbs = file2output.get(file);
            ps.println("Identified in " + file);
            for(int i = 0 ; i < MAXSCANNUM; i++) {
                if(sbs[i] != null) {

                    ps.println(i + sbs[i].toString());
                }
            }
            ps.println("\n\n\n");
        }
        //ps.println(sb.toString()); 
    }
    private static StringBuffer getOutput(String ms2file, int precursorscannum) {
//System.out.println("precscan: " + precursorscannum);
        StringBuffer sbs [] = file2output.get(ms2file);
        if(sbs == null) {
            sbs = new StringBuffer[MAXSCANNUM];
            file2output.put(ms2file, sbs);
        }
        if(sbs[precursorscannum] == null) {
            sbs[precursorscannum] = new StringBuffer(100);
        }
        return sbs[precursorscannum];
    }
    private static int getPrecursorScanNumber(String ms2file, int ms2scannum) throws IOException {
        int[] precurscannums = file2precursorscans.get(ms2file);
        if(precurscannums == null) {
            precurscannums = getPrecursorScanNumbers(ms2file); 
            file2precursorscans.put(ms2file, precurscannums);
        }
        return precurscannums[ms2scannum]; 
    }

    private static int [] getPrecursorScanNumbers(String ms2file) throws IOException {
        int [] precursorscannums = new int[MAXSCANNUM];
        SpectrumReader sr = new SpectrumReader(ms2file+".ms2", "ms2"); 
        for(Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList spectrum = it.next();
            int scannum = spectrum.getLoscan();
            int precursorscan = spectrum.getPrecursorScanNumber();
            precursorscannums[scannum] = precursorscan;

        }
        return precursorscannums;
    }
    public static void main(String args[]) throws Exception {
   

        // Command line processing
        Options opts = new Options();
        Option inputOpt = new Option
            ("i", "input", true, "Input file name - DTASelect-filer.txt run with --DM option");
        Option pValueOpt = new Option
            ("p", "ppm", false,  "delta mass in ppm or absolute value");
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
            outputppm = true;        
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
