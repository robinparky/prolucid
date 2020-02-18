/**
 * @file IdentifiedSpectraSeparator.java
 * This is the source file for edu.scripps.pms.util.spectrum.IdentifiedSpectraSeparator
 * @author Tao Xu
 * @date $Date: 2009/10/09 21:02:03 $
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
public class IdentifiedSpectraSeparator {
    public static final String USAGE = "\n\n!!! USAGE: identifiedpeptidedb -i input -o output !!!";
    private static ArrayList<HashSet<String>> peptideSets = new ArrayList<HashSet<String>>();
//    private static HashSet<String> peptidesInCommon = new HashSet<String>();
    private static String inputFile = null;
    private static String outputFile = null;
    private static double [] avgAaMasses = new double[256];
    private static double [] monoAaMasses = new double[256];
    private static HashMap<String, HashSet<Integer>> filename2scans = new HashMap<String, HashSet<Integer>>(1000000);   
 
    public static ArrayList<Protein> getDTASelectProteins(String dtaselectFilterFile) throws IOException {

        ArrayList<Protein> result = new ArrayList(10000);
        DTASelectFilterReader reader = new DTASelectFilterReader(dtaselectFilterFile);
        for (Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {
            result.add(itr.next());
        }
        reader.close(); 
        return result;
    }
    // need to use --extra -DM -t 0 options in DTASelect2
    public static void getIdentifiedSpectra(ArrayList<Protein> dtaselectProteins) throws IOException {
        
        for(Protein p : dtaselectProteins) {
            int peptidecount = 1;
            String acc = p.getAccession();
            if(!acc.startsWith("Re")) {
                //HashSet<String> peptideadded = new HashSet<String>(100000); 
                for(Iterator<Peptide> pepItr=p.getPeptides(); pepItr.hasNext(); ) {
              
                    Peptide peptide = pepItr.next();
                    String ms2file = peptide.getFileName();
                    Integer scan = new Integer(peptide.getLoScan());
                    HashSet scans = filename2scans.get(ms2file);
                    if(scans == null) {
                        scans = new HashSet<Integer>(1000000);
                        filename2scans.put(ms2file, scans);
                    }
                    scans.add(scan);

                }
            }
        }
        
    }

    public static void outputSpectra(String ms2file) throws IOException {
        SpectrumReader sr = new SpectrumReader(ms2file+".ms2", "ms2");
        HashSet<Integer> identifiedscans = filename2scans.get(ms2file);
        PrintStream identified = new PrintStream(new File("identified/" + ms2file+"_identified.ms2"));
        PrintStream unidentified = new PrintStream(new File("identified/" + ms2file+"_unidentified.ms2"));
        boolean isfirstidentified = true;
        boolean isfirstunidentified = true; 
        for(Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList spectrum = it.next();
            int scan = spectrum.getLoscan();
            
            if(identifiedscans != null && identifiedscans.contains(new Integer(scan))) {
                if(isfirstidentified) {
                    identified.print(spectrum.getSpectrumWithHlines());
                    isfirstidentified = false;
                } else {
                    identified.print(spectrum.getSpectrumWithoutHlines());
                }
            } else {
                if(isfirstunidentified) {
                    unidentified.print(spectrum.getSpectrumWithHlines());
                    isfirstunidentified = false;
                } else {
                    isfirstunidentified = false;
                    unidentified.print(spectrum.getSpectrumWithoutHlines());
                }

            }
        }
        identified.close();
        unidentified.close();
    }
    public static ArrayList<String> getMs2Files(String dir) {

         ArrayList<String> ms2Files = new ArrayList<String>();
         File currentDir = new File(dir);
         for(String s : currentDir.list()) {
             if(s.endsWith(".ms2")) {
                 //System.out.println("got ms2 file: " + s);
                 ms2Files.add(s.split(".ms2")[0]);
             }
         }
         return ms2Files;
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
            getIdentifiedSpectra(getDTASelectProteins(inputFile));
            ArrayList<String> ms2files = getMs2Files(".");
            //Runtime.getRuntime().exec("mkdir identified");
            new File("identified").mkdir();
            for(Iterator<String> it = ms2files.iterator(); it.hasNext();) {
                String ms2file = it.next();
                System.out.print("Processing " + ms2file + ".ms2 file ...");
                outputSpectra(ms2file);
                System.out.println();
            }

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
