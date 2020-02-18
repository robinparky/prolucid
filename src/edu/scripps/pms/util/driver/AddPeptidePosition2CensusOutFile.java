/**
 * @file  AddPeptidePosition2CensusOutFile.java
 * This is the source file for edu.scripps.pms.util.spectrum. AddPeptidePosition2CensusOutFile
 * @author Tao Xu
 * @date $Date: 2011/08/30 19:56:11 $
 */



import edu.scripps.pms.util.io.*;

//import edu.scripps.pms.util.io.DTASelectFilterReader; 
import java.io.*;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;

// command line processor
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

// this is the class to add peptide postion to the census out file.
// requires two input files DB-Prot2Pep.txt census-out.txt
public class  AddPeptidePosition2CensusOutFile {
    public static final String USAGE = "\n\n!!! USAGE: addPeptidePosition2CensusOutFile -p DB-Prot2Pep.txt -c census-out.txt !!!";
//    private static HashSet<String> peptidesInCommon = new HashSet<String>();
    private static String prot2PepFile = null;
    private static String censusOutFile = null;
    private static String databasename;
    private static HashMap<String, String> spectrum2position = new HashMap(1000000);
    
    public static void getPeptidePositions(String prote2pepFile) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(prote2pepFile));
        String line = br.readLine(); // remove the header line in the DB-Prot2Pep.txt file
        while(line != null) {
            String [] arr = line.split("\t");
            if(arr.length > 2) {
                spectrum2position.put(arr[1], arr[2]);
            }         
            line = br.readLine();
        }
        br.close(); 
    }
    public static void output() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(censusOutFile));
        String line = br.readLine();
        while(line != null) {
            char c = line.charAt(0);
            if(c == '&' || c == 'S') {
                String [] arr = line.split("\t");
                String spectruminfo = "";
                if(c == 'S') {
                    spectruminfo = arr[13] + "." + arr[14] + "." + arr[14] + "."  + arr[15]; 
                } else {
                    spectruminfo = arr[12] + "." + arr[13] + "." + arr[13] + "."  + arr[14]; 

                }
//System.out.println(spectruminfo);
                System.out.println(line + "\t" + spectrum2position.get(spectruminfo));
                
            } else {

                System.out.println(line); 
            }
            line =  br.readLine(); 
        }
        br.close();
       
    }
    public static void main(String args[]) throws Exception {
   

        // Command line processing
        Options opts = new Options();
        Option pepOpt = new Option
            ("p", "DB-Prot2Pep.txt", true, "DB-Prot2Pep.txt from DTASelect --DB option");
        Option censusoutOpt = new Option
            ("c", "census-out.txt", true,  "Census Output file");
        //pepOpt.setRequired(true);
        //pValueOpt.setRequired(true);
        //hprdOpt.setRequired(false);
        //opts.addOption(dbNameOpt);
        opts.addOption(pepOpt);
        opts.addOption(censusoutOpt);

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
        prot2PepFile = cli.getOptionValue("p");
        if(prot2PepFile == null) {
            prot2PepFile = "DB-Prot2Pep.txt";
        }
        censusOutFile = cli.getOptionValue("c");
        if(censusOutFile == null) {
            censusOutFile = "census-out.txt";
        }

        try {
            //System.out.println("p value cutoff: " + minPValue);
            getPeptidePositions(prot2PepFile);
            output();

        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
        }
    }


}
