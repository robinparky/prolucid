
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.mspid.ProteinDatabase;
import edu.scripps.pms.util.io.FastaReader;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.List;
import java.util.LinkedList;
import java.util.HashSet;
import java.io.IOException;
import java.io.FileReader;
import java.io.FileInputStream;
import java.io.PrintStream;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;


/**
 * @author  Tao Xu
 * @version $Id
 */
public class IPIId2RefseqAccession {
    public static String USAGE = "Usage: java -Xmx500M IPIId2RefseqAccession IPI_accession_file fasta_file";
    private static HashSet<String> ipiacs = new HashSet<String>(1000000); // swissprot of trembl id
    private static String fastadb = "/garibaldi/people-b/applications/yates/dbase/EBI-IPI_mouse_3.30_06-28-2007_con_reversed.fasta"; 
    private static String ipiFileName;
    private static HashMap<String, String> ipiac2refseqac = new HashMap<String, String>();

    public static void main(String args[]) throws IOException {
        try { 
            ipiFileName = args[0];
            if(args.length > 1) {
                fastadb = args[1];
            }
    
            getIdentifiedIds(); 
            int numrefseqac = 0;
            int numipiac = 0;
            for (Iterator<Fasta> itr = FastaReader.getFastas(new FileInputStream(fastadb)); itr.hasNext();) { 
                Fasta f = itr.next();
                String ipiac = f.getAccession();
                if(ipiacs.contains(ipiac)) {
                    //System.out.println("found " + ipiac + " refseqac: " + f.getRefSeqAccession());
                    String ac = f.getRefSeqAccession();
                    ac = ac == null? f.getSwissprotAccession() : ac;
                    ac = ac == null? f.getTremblAccession() : ac;
                    ac = ac == null? f.getEnsemblsAccession() : ac;
                    ac = ac == null? f.getVegaAccession() : ac;
                    ipiac2refseqac.put(ipiac, ac);
            
                }
            }

            output();        
       } catch(Exception e) {
           System.err.println(USAGE);
       }

      // System.out.println("NumAdded: " + numAdded);
    }

    public static void output() throws IOException {
       
        String prefix = ipiFileName.substring(0, ipiFileName.lastIndexOf("."));
        PrintStream ps = new PrintStream(prefix + "_with_refseq.txt");
        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(ipiFileName)));
        String line = null;
        while ((line = br.readLine()) != null) {
            String ipiac = line.trim();
            String refseqac = ipiac2refseqac.get(ipiac);
            if(refseqac != null) {
                ps.println(ipiac + "\t" + refseqac);
            } else {
                ps.println(ipiac + "\t" + "");
            } 
        }
        br.close();
        ps.close();
//System.out.println("Number of ids: " + ipiacs.size());
    }
    public static void getIdentifiedIds() throws IOException {


        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(ipiFileName)));
        //br.readLine(); // remove the first line
        String line = null;
        while ((line = br.readLine()) != null) {
            ipiacs.add(line.trim());
        }
        br.close();
//System.out.println("Number of ids: " + ipiacs.size());
    }
}
