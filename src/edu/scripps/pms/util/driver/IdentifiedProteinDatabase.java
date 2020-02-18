/**
 * @file  IdentifiedProteinDatabase.java
 * This is the source file for edu.scripps.pms.util.spectrum. IdentifiedProteinDatabase
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
import edu.scripps.pms.mspid.ProteinDatabase;

// command line processor
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

// this is the class for calculating the z-score for delta mass 
public class  IdentifiedProteinDatabase {
    public static final String USAGE = "\n\n!!! USAGE: identifiedproteindb -i input -o output !!!";
//    private static HashSet<String> peptidesInCommon = new HashSet<String>();
    private static String inputFile = null;
    private static String outputFile = null;
    private static String databasename;
    private static HashSet<String> identifiedlocus = new HashSet<String>(100000);
    
    public static void getDTASelectProteins(String dtaselectFilterFile) throws IOException {

        DTASelectFilterReader reader = new DTASelectFilterReader(dtaselectFilterFile);
        databasename = reader.getDbFilePathAndName();
        for (Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {
            identifiedlocus.add(itr.next().getAccession());
        }
        reader.close(); 
    }
    public static void output() throws IOException {
        
        PrintStream ps = null;
        if(outputFile == null) {
            ps = System.out;
        } else {
            ps = new PrintStream(outputFile);
        }
       
        for(Iterator<Fasta> it = FastaReader.getFastas(databasename); it.hasNext();) {
             
            Fasta f = it.next(); 
            String acc = f.getAccession(); 
            String sequestacc = f.getSequestLikeAccession(); 
            //System.out.println("acc: " + acc + "\tseqeustacc: " + sequestacc);
            if(!f.isReversed() && (identifiedlocus.contains(acc) || identifiedlocus.contains(sequestacc))) {
                ps.println(">" + f.getDefline());
                ps.println(f.getSequence());
            }
        }
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
            getDTASelectProteins(inputFile);
            output();

        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
        }
    }


}
