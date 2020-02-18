/**
 * @file Ms12Ms2.java
 * This is the source file for edu.scripps.pms.util.spectrum.Ms12Ms2
 * @author Tao Xu
 * @date $Date: 2009/06/23 23:27:48 $
 */



import edu.scripps.pms.util.io.*;

//import edu.scripps.pms.util.io.DTASelectFilterReader; 
import edu.scripps.pms.util.PmsUtil; 
import java.io.*;
import java.util.List;
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
public class Ms12Ms2 {
    public static final String USAGE = "\n\n!!! USAGE: ms12ms2 -i input -m m/z -z charge -f fragmentationMethod !!!";
//    private static HashSet<String> peptidesInCommon = new HashSet<String>();
    private static String inputFile = null;
    private static String fragmentationMethod = null;
    private static double mz = 0;
    private static int chargeState = 2;
 
    private static final double MASSH = 1.007276;;

    public static void outputSpectra(String ms1file) throws IOException {
        SpectrumReader sr = new SpectrumReader(ms1file+".ms1", "ms1");
        PrintStream ms2file = new PrintStream(new File(ms1file +"_from_ms1.ms2"));
        boolean isfirstspectrum = true;
        for(Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList spectrum = it.next();
            int scan = spectrum.getLoscan();
            double mplush = chargeState*mz - chargeState*MASSH + MASSH;           
            String frgmline = "I\tActivationType\t" + fragmentationMethod; 
 
            spectrum.addZline(new Zline(chargeState, mplush));
            spectrum.addIline(frgmline);
            spectrum.setPrecursorMass(mz);

            if(isfirstspectrum) {
                ms2file.print(spectrum.getSpectrumWithHlines());
                isfirstspectrum = false;
            } else {
                ms2file.print(spectrum.getSpectrumWithoutHlines());
            }

        }
        ms2file.close();
    }
    public static ArrayList<String> getMs1Files(String dir) {

         ArrayList<String> ms1files = new ArrayList<String>();
         File currentDir = new File(dir);
         for(String s : currentDir.list()) {
             if(s.endsWith(".ms1")) {
                 //System.out.println("got ms2 file: " + s);
                 ms1files.add(s.split(".ms1")[0]);
             }
         }
         return ms1files;
    }


    public static void main(String args[]) throws Exception {
   

        // Command line processing
        Options opts = new Options();
        Option inputOpt = new Option
            ("i", "input", true, "Input ms1 file name");
        Option frgmOpt = new Option
            ("f", "fm", true,  "Fragmentation method");
        Option mzOpt = new Option
            ("m", "mz", true,  "m/z");
        Option zOpt = new Option
            ("z", "chargestate", true,  "Charge state");

        zOpt.setRequired(true);
        mzOpt.setRequired(true);
        //inputOpt.setRequired(true);
        //pValueOpt.setRequired(true);
        //hprdOpt.setRequired(false);
        //opts.addOption(dbNameOpt);
        opts.addOption(inputOpt);
        opts.addOption(frgmOpt);
        opts.addOption(mzOpt);
        opts.addOption(zOpt);

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
        fragmentationMethod = cli.getOptionValue("f");
        fragmentationMethod = fragmentationMethod == null? "CID" : fragmentationMethod;
        mz = Double.parseDouble(cli.getOptionValue("m"));
        chargeState = Integer.parseInt(cli.getOptionValue("z"));
        
        inputFile = cli.getOptionValue("i");
        if(inputFile != null) {
            String file = inputFile.split(".ms1")[0];
            outputSpectra(file);
            return;
        }


        try {
            //System.out.println("p value cutoff: " + minPValue);
            ArrayList<String> ms1files = getMs1Files(".");
            for(Iterator<String> it = ms1files.iterator(); it.hasNext();) {
                String ms1file = it.next();
                System.out.print("Processing " + ms1file + ".ms1 file ...");
                outputSpectra(ms1file);
                System.out.println();
            }

        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
        }
    }



}
