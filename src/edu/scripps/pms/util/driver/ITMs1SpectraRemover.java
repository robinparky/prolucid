/**
 * @file ITMs1SpectraRemover.java
 * This is the source file for edu.scripps.pms.util.spectrum.ITMs1SpectraRemover
 * @author Tao Xu
 * @date $Date: 2011/01/29 00:13:58 $
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

// this is the class to remove mistaken itms ms1 scans in FTMS ms1 file
public class ITMs1SpectraRemover {
    public static final String USAGE = "\n\n!!! itms1spectraremover -i ms1file !!!";
//    private static HashSet<String> peptidesInCommon = new HashSet<String>();
    private static String instrumentType = "ITMS";
 

    public static void outputSpectra(String ms1file) throws IOException {
        
        String tempfile = ms1file +"_from_ms1";
        String inputfile = ms1file + ".ms1";
        System.out.print("Processing " + inputfile + " file ...");
        SpectrumReader sr = new SpectrumReader(inputfile, "ms1");
        PrintStream tempoutput = new PrintStream(new File(tempfile));
        boolean isfirstspectrum = true;
        for(Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList spectrum = it.next();
            if(instrumentType.equals(spectrum.getInstrumentType())) continue; 

            if(isfirstspectrum) {
                tempoutput.print(spectrum.getSpectrumWithHlines());
                isfirstspectrum = false;
            } else {
                tempoutput.print(spectrum.getSpectrumWithoutHlines());
            }

        }
        tempoutput.close();
        
        //Runtime.getRuntime().exec("rm -rf " + inputfile);
        //String rename = "mv " + tempfile + " " + inputfile;
        //System.out.println("now, " + rename);
        //Runtime.getRuntime().exec(rename);
        new File(inputfile).delete();
        new File(tempfile).renameTo(new File(inputfile));
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
        opts.addOption(inputOpt);

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
        
        String inputFile = cli.getOptionValue("i");
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
