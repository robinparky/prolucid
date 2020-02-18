
import edu.scripps.pms.util.seq.Fasta;
//import edu.scripps.pms.mspid.ProteinDatabase;
import edu.scripps.pms.util.seq.UniProtProtein;
import edu.scripps.pms.util.io.UniprotDatReader;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.FileInputStream;
import java.util.Iterator;
//import java.util.LinkedList;
import java.io.IOException;
import java.io.PrintStream;


public class UniProtDat2Fasta {
  
    private static String usage = "!!! USAGE: uniprot2fasta  input_uniprot_dat_file_name output_fasta_file_name !!!";
    public static void main(String args[]) throws IOException {
     
        int numAdded = 0;
        String inputFile = null;
        String outputFile = null;
        if(args.length < 2) {
            System.err.println(usage);
            System.exit(1);
        }
        try {

            inputFile = args[0];
            outputFile = args[1];


            
            PrintStream ps = new PrintStream(outputFile);
            for (Iterator<UniProtProtein> itr = UniprotDatReader.getUniProtProteins(new FileInputStream(inputFile)); itr.hasNext(); ) { 
                UniProtProtein protein = itr.next();
                ps.println(protein.outputFasta());
                numAdded++;
            }
        } catch (Exception e) {
            System.err.println(usage);
            e.printStackTrace();
            System.err.println(usage);
            System.exit(1);
        }
        System.out.println("Number of proteins processed: " + numAdded);
    }

}
