/**
 * @file InclusionList.java
 * This is the source file for edu.scripps.pms.util.spectrum.InclusionList
 * @author Tao Xu
 * @date $Date
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

// for generating the inclusion list and exclusion list 
// list file contains a list of common separted values, each value on a line ended with a ","
// the list file can be either inclusion list or exclusion list
// for each value in the list file, the program output a line with two columns,
// the first is the mz value in the list, and the second will be 1 if the mz appeared in the given segment of the given ms2 file
// use -s 0 and -e 100000000 if there is just one segment
// use -m to specify the isolation mass tolerance
public class InclusionListChecker {
    public static final String USAGE = "\n\n!!! USAGE: inclusionlistchecker -l list_file -m ms2file -s start_time -e end_time -t mass_tolerance";
    public static final double MASSDIFFC12C13 = 1.003354826;
    public static final double TWOMASSDIFFC12C13 = MASSDIFFC12C13*2;
    public static final double THREEMASSDIFFC12C13 = MASSDIFFC12C13*3;
    

    
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
            ("l", "list_file", true, "inclusion or exclusion list file");
        Option ms2Opt = new Option
            ("m", "ms2_file", true, "ms2 file name");
        Option startOpt = new Option
            ("s", "start_time", true,  "segment start time in minutes");
        Option endOpt = new Option
            ("e", "end_time", true,  "segment end time in mimutes");
        Option toleranceOpt = new Option
            ("t", "mass_tolerance", true,  "mass tolerance in amu");
        opts.addOption(inputOpt);
        opts.addOption(ms2Opt);
        opts.addOption(startOpt);
        opts.addOption(endOpt);
        opts.addOption(toleranceOpt);

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

        try {
            String listfile = cli.getOptionValue("l");
            String ms2file = cli.getOptionValue("m");

            double starttime = Double.parseDouble(cli.getOptionValue("s"));
            double endtime = Double.parseDouble(cli.getOptionValue("e"));
            double tolerance = Double.parseDouble(cli.getOptionValue("t"));
                
            System.out.print("Processing " + ms2file + " file ...");

            SpectrumReader ms2sr = new SpectrumReader(ms2file, "ms2");

            boolean [] isselected = new boolean[4000000];

            BufferedReader br = new BufferedReader(new FileReader(listfile));
            String line = null; 
            while((line = br.readLine()) != null) {
                String mzl = line.trim().split(",")[0];
                if(mzl != null) {
                    double mz = Double.parseDouble(mzl);
                    //int mzint = (int)(mz*1000);
                    int lowlimit = (int)((mz - tolerance)*1000);
                    int highlimit = (int)((mz + tolerance)*1000);
                    for(int i = lowlimit; i <= highlimit; i++) {
                        if(i >= 0) {
                            isselected[i] = true;
                        }
                    }

                }
            }

            ArrayList<PeakList> ms2spectra = ms2sr.getSpectraList();
            for(Iterator<PeakList> it = ms2spectra.iterator(); it.hasNext();) {
                PeakList spectrum = it.next();
                double ret = spectrum.getRetentionTime();            
                if(ret < starttime) {
                    continue;
                }
                if(ret > endtime) {
                    break;
                }

                //int scan = spectrum.getLoscan();
                //int chargestate = spectrum.getFirstChargeState();
                double precmz = spectrum.getPrecursorMass();            
                int mzint = (int)(precmz*1000);
                if(isselected[mzint]) {
                    System.out.println("1\t" + precmz + "\ton list"); 
                } else {

                    System.out.println("0\t" + precmz + "\tnot on list");
                }
            }

            ms2sr.closeDataFile();
            br.close();

        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
        }
    }
    

}
