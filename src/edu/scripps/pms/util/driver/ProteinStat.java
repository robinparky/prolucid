

/**
 * @file ProteinStat.java
 * Application to load PSI-MI files  to the database.
 * @author Tao Xu 
 * @date $Date
 */

// log4j
import org.apache.log4j.*;


import edu.scripps.pms.util.io.*;
import edu.scripps.pms.util.dtaselect.Peptide; 
import edu.scripps.pms.util.dtaselect.Protein; 
import edu.scripps.pms.util.dtaselect.ProteinItem; 
// command line processor
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

// java util
import java.io.*;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.HashSet;
import java.util.HashMap;


public class ProteinStat {

    private static Logger logger = Logger.getLogger(ProteinStat.class);

    public static final String USAGE = "java ProteinStat [-i inputFile/inputPath]";
            //" [-d database name] {-h | -p | -m}\n";

    
    public static void main (String[] args) throws IOException {

        // Set up a simple configuration that logs to STDERR.
        BasicConfigurator.configure
            (new ConsoleAppender
             (new PatternLayout(PatternLayout.TTCC_CONVERSION_PATTERN),
              ConsoleAppender.SYSTEM_ERR));
        // Set the default logging level and any exceptions.
         Logger.getRootLogger().setLevel(Level.FATAL);



        // Command line processing
        Options opts = new Options();
        Option inputOpt = new Option
            ("i", "input", true, "Input file that contains the directories");
        Option freqOpt = new Option
            ("f", "frequenceCutOff", true,  "Occurence in multiple experiments");
        Option subdirOpt = new Option
            ("s", "subdirecotry", true,  "Subdirectory name");
        Option seqCovOpt = new Option
            ("c", "seqCoverage", true,  "Sequence coverage cutoff");
        Option svmModelOpt = new Option
            ("m", "svmModel", true, "SVM model file");
        Option svmTrainOpt = new Option
            ("t", "svmTrain", false, "SVM training mode");
        //Option hprdOpt = new Option("h", "hprd", false, "Select this option to load HPRD PSI-MI files.");
        inputOpt.setRequired(true);
        freqOpt.setRequired(true);
        subdirOpt.setRequired(true); // the option is required 
        seqCovOpt.setRequired(true);
        svmModelOpt.setRequired(true);
        svmTrainOpt.setRequired(false);
        //hprdOpt.setRequired(false);
        //opts.addOption(dbNameOpt);
        opts.addOption(inputOpt);
        opts.addOption(freqOpt);
        opts.addOption(subdirOpt);
        opts.addOption(seqCovOpt);
        opts.addOption(svmModelOpt);
        opts.addOption(svmTrainOpt);

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

        //dbName = cli.getOptionValue("d");
        String inputPath = cli.getOptionValue("i");
        String subdir = cli.getOptionValue("s");

        int freqCutOff = Integer.parseInt(cli.getOptionValue("f"));        
        int seqCov = Integer.parseInt(cli.getOptionValue("c"));        
        int numPeptideCutOff = 1; 
        int compScore = 1; 
        // end command line processing
        HashMap<String, ProteinItem> map = getProteins(getDtaselectFilterFiles(inputPath, subdir));
        System.out.println("Number of proteins idenentified: " + map.size());
        double fpr = 1; // false positive rate
/*
        while (fpr > 0.05) {
            int numReverse = 0;
            int numPassed = 0;
            
            for(String s : map.keySet()) {
                ProteinItem pi = map.get(s);
                if(!pi.isProblematic()) {
                    //if(pi.getFrequency() > freqCutOff && pi.getNumPeptides() > numPeptideCutOff && pi.getSeqCoverage() >= seqCov) {
//System.out.println("compositeScore: " + pi.getCompositeScore());
                    if(pi.getCompositeScore() > compScore) {
                        if(pi.getAccession().startsWith("R")) {
                            numReverse++;
                        } else {
                            numPassed++;            
                        }
                    }
                }
            } 
        
            fpr = numReverse/(numPassed + numReverse + 0.0);
            System.out.println("NumPassed: " + numPassed + "\tNumReverse: " + numReverse + "\tfrequence: " + freqCutOff + "\tfp: " + fpr);
            //freqCutOff++;
            //seqCov++;
            //numPeptideCutOff++;
            //compScore += 10;
            compScore += 5;
            
        }
*/
        ArrayList<ProteinItem> forward = new ArrayList<ProteinItem>(10000);
        ArrayList<ProteinItem> reverse = new ArrayList<ProteinItem>(10000);

        int numTotal = 0;
        for(String s : map.keySet()) {
            ProteinItem p = map.get(s);
            if(!p.isProblematic()) {
                numTotal++;
            }
            //if((!p.isProblematic()) && p.getCompositeScore() > compScore) {
            if((!p.isProblematic())) {
             //   System.out.println(p.getAccession() + "\t" + p.getCompositeScore() + "\t" + p.getFrequence() + 
             //                      "\t" + p.getNumPeptides() + "\t" + p.getSeqCoverage()+ "\t" + p.getLength());
                if(p.getAccession().startsWith("Rev")) {
                    reverse.add(p);
                } else {
                    forward.add(p);
                }
            } 
        }  
        System.out.println("Number of proteins identified: " + numTotal);
        System.out.println(ProteinItem.getHeaderLine());
        int numReverse = 0;
        StringBuffer sb = new StringBuffer(1000000);

        for(ProteinItem f : forward) {
            sb.append(f.output());

            //if(numReverse < reverse.size()) {
            //    sb.append("\t");
            //    sb.append(reverse.get(numReverse++).output());
            //}
            sb.append("\n");
        }
       
        // output reverse hits
        for(ProteinItem r : reverse) {
            sb.append(r.output());
            sb.append("\n");

        }
        
        
        System.out.println(sb.toString());
    }
    private static Set<String> getDtaselectFilterFiles(String inputFile, String subdir) throws IOException {
         HashSet<String> files = new HashSet<String>();
         BufferedReader  br = new BufferedReader(new FileReader(inputFile));
         String s = "";
         while((s = br.readLine()) != null) {

             if(!s.equals("")) {
                 
                 //String file  = s.split("\t")[0] + "/" + subdir + "/notrypstat/DTASelect-filter.txt"; 
                 String file  = s.split("\t")[0] + "/" + subdir + "/DTASelect-filter.txt"; 
                 files.add(file);
             }
         }
         br.close(); 
         //Collections.sort(files);
         return files;
    }
    private static ArrayList<Protein> getProteins(String dtaselectFilterFile) throws IOException {

        DTASelectFilterReader reader = new DTASelectFilterReader(dtaselectFilterFile);
        ArrayList<Protein> proteins = new ArrayList<Protein>(2000);
        ArrayList<Protein> noPeptideProteins = null; // to keep track of the protein without peptides
        boolean problemProteinFound = false;
        for (Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {
            Protein protein = itr.next();
            if(protein.getPeptideSize() < 1) {
                if(!problemProteinFound) {
                    noPeptideProteins = new ArrayList<Protein>();
                    problemProteinFound = true; 
                }
                
                noPeptideProteins.add(protein);
            } else {
                if(problemProteinFound) {
                    for(Iterator<Peptide> itp = protein.getPeptides(); itp.hasNext();) {
                        Peptide p = itp.next();
                        for(Iterator<Protein> itpp = noPeptideProteins.iterator(); itpp.hasNext();) {
                            itpp.next().addPeptide(p);
                        }
                      
                    }
                    problemProteinFound = false; 
                } 
            }
            
            proteins.add(protein);
            String temps = "IPI00376164.2";
            if(protein.getAccession().startsWith(temps)) {
                System.out.println("in file " + dtaselectFilterFile + "found protein " + temps);
            }
        }
        return proteins; 
    }
    private static HashMap<String, ProteinItem> getProteins(Set<String> dtaselectFilterFiles) throws IOException {
        HashMap<String, ProteinItem> acc2ProteinItems = new HashMap<String, ProteinItem>(3000); 
        for(Iterator<String> it = dtaselectFilterFiles.iterator(); it.hasNext();) {
            String dtaselectFilterFile = it.next();
            System.out.println("Now processing " + dtaselectFilterFile);
            try {
                
                for (Iterator<Protein> itr = getProteins(dtaselectFilterFile).iterator(); itr.hasNext();) {
                    
                    Protein protein = null;
                    protein = itr.next();
                    String acc = protein.getAccession();
                    ProteinItem pi = acc2ProteinItems.get(acc);
                    if(pi == null) {
                        pi = new ProteinItem(protein);
                        acc2ProteinItems.put(acc, pi); 
                    } else {
                        pi.addProtein(protein);
                    }
                    for(Iterator<Peptide> pepItr=protein.getPeptides(); pepItr.hasNext(); ) {
                        Peptide peptide = pepItr.next();
                        //peptideSets.add(peptides);
                        //System.out.println(peptides.size() + " peptides found in " + dtaselectFilterFile);
                        //Integer lowSan = new Integer(peptide.getLoScan());
                        //System.out.println(peptide.getSequence());
                    }
                }
            } catch(IOException e) {
                System.out.println(dtaselectFilterFile + " file does not exist");
                //e.printStackTrace(); 
            }
        }
        return acc2ProteinItems; 
    }
   
}
