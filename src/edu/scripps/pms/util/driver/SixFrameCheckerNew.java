/**
 * @file  SixFrameCheckerNew.java
 * This is the source file for edu.scripps.pms.util.spectrum. SixFrameCheckerNew
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
public class  SixFrameCheckerNew {
    public static final String USAGE = "\n\n!!! USAGE: java SixFrameCheckerNew -i input -o output !!!";
//    private static HashSet<String> peptidesInCommon = new HashSet<String>();
    private static String inputFile = null;
    private static String outputFile = null;
    private static String databasename;
    private static HashSet<Protein> identifiedlocus = new HashSet<Protein>(100000);
    //private static HashMap<String, HashSet<String>> acc2frames = new HashMap<String, HashSet<String>>(10000);
    private static HashMap<String, HashSet<Protein>> locus2Proteins = new HashMap<String, HashSet<Protein>>(10000);
    private static HashMap<String, String> locus2Description = new HashMap<String, String>(10000);
    private static HashSet<String> identifiedaccs = new HashSet<String>(10000); // for accession with frame
    private static HashMap<String, Fasta> acc2Fastas = new HashMap<String, Fasta>(100000); 

    public static void getDTASelectProteins(String dtaselectFilterFile) throws IOException {

        ProteinDatabase pdb = new ProteinDatabase("/data/3/taoxu/db/rat_refeq_01192009/rat.rna.fna");
        DTASelectFilterReader reader = new DTASelectFilterReader(dtaselectFilterFile);
        databasename = reader.getDbFilePathAndName();
        for (Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {
            identifiedlocus.add(itr.next());
        }
        reader.close(); 

        for(Iterator<Protein> it = identifiedlocus.iterator(); it.hasNext();) {
            Protein p = it.next();
            String acc = p.getAccession();
            identifiedaccs.add(acc);
            char c = acc.charAt(0);
            String frame = acc.substring(0, 2);
            String locus = acc.substring(2, acc.length());
            
            //if(c == '-' || c == '+') {
            if(true) {
                //System.out.println(acc + "\t" + locus + "\t" + frame);
                HashSet<Protein> proteins = locus2Proteins.get(locus);
                if(proteins == null) {
                    proteins = new HashSet<Protein>(6);
                    locus2Proteins.put(locus, proteins);
//System.out.println(locus + "\t" + acc);
                    Fasta f = pdb.accession2Fasta(locus);
                    
                    String defline = f == null? "" : f.getDefline();
                    locus2Description.put(locus, defline);
                }
                proteins.add(p);

            }
        }
        pdb = null;
        System.gc();

       
        FileInputStream fis = new FileInputStream(new File(databasename));
        Iterator<Fasta> fastas = FastaReader.getFastas(fis);
        int i = 0;
        while(fastas.hasNext()) {
            //System.out.print("\r" + ++i);
            Fasta f = fastas.next();
            String acc = f.getAccession();
            if(identifiedaccs.contains(acc)) {
//System.out.println("found fasta\t" + acc + "\t" + f.getSequence());
               acc2Fastas.put(acc, f);
            }
            //System.out.println("ac: " + f.getAccession() + "\tdefline: " + f.getDefline());
        }
        fis.close();
    }
    public static void output() throws IOException {
        
        PrintStream ps = null;
        if(outputFile == null) {
            ps = System.out;
        } else {
            ps = new PrintStream(outputFile);
        }
        //ProteinDatabase pdb = new ProteinDatabase("/data/3/taoxu/db/rat_refeq_01192009/rat.rna.fna");
        int numMultipleFrames = 0; 
        int numGenes = 0; 
        ps.println("Accession\tnumFrames\tFrame\tPeptideCount\tUniquePeptideCount\tSpectrumCount\tminindex\tmaxindex\tFrame\tPeptideCount\tUniquePeptideCount\tSpectrumCount\tdescription");
        for(Iterator<String> it = locus2Proteins.keySet().iterator(); it.hasNext();) {
            String acc = it.next(); 
            numGenes++; 
            HashSet<Protein> proteins = locus2Proteins.get(acc);
            int numFrames = proteins.size();

            if(numFrames > 1) {
                numMultipleFrames++;
                String description = locus2Description.get(acc);
                ps.print(acc + "\t" + numFrames);
                for(Iterator<Protein> pit = proteins.iterator(); pit.hasNext();) {
                    Protein p = pit.next();
                    //ps.print("\t" + p.getAccession() + "\t" + p.getSeqCount() + "\t" + p.getSpectrumCount()); 
                    int numuniquepeptides = 0;
                    String protseq = acc2Fastas.get(p.getAccession()).getSequence();
// System.out.println(acc + "\t" + protseq);                    
                    int minindex = 100000000;
                    int maxindex = -1;
                    
                    for(Iterator<Peptide> pepit = p.getPeptides(); pepit.hasNext();) {
                        Peptide pp = pepit.next();
                        String pepseq = pp.getMidSeq();
//System.out.println(pepseq + "\t" + p.getAccession());
                        int index = protseq.indexOf(pepseq);
                        int lastindex = index + pepseq.length();                        
                        minindex = index != -1 && minindex > index? index : minindex;
                        maxindex = index != -1 && maxindex < lastindex? lastindex : maxindex; 
                        if(pp.isUnique()) {
                            numuniquepeptides++;
                        }
                    }
                    ps.print("\t" + p.getAccession() + "\t" + p.getSeqCount() + "\t" + numuniquepeptides + "\t" + p.getSpectrumCount()+ "\t" + minindex + "\t" + maxindex); 
                } 
                
               ps.println("\t" + description);
            }
        }
        ps.println("Number of genes identified:\t" + numGenes);
        ps.println("Number of genes that have multiple frames identified:\t" + numMultipleFrames);
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
