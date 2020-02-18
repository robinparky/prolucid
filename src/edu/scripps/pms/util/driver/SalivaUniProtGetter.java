
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.seq.UniProtProtein;
import edu.scripps.pms.util.io.UniprotDatReader;
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


/**
 * @author  Tao Xu
 * @version $Id
 */
public class SalivaUniProtGetter {
    private static HashSet<String> identifiedProteinIds = new HashSet<String>(1000000);
    private static String uniprotdatfile = "/data/6/xmhan/database/uniprot_sprot_human_11-14-2007.dat";
    private static String humanFasta = "/data/6/xmhan/database/EBI-IPI_human_3.37_12-05-2007_original.fasta";
    private static String fastaoutputfile = "/data/2/taoxu/projects/xuemei/database/identifiedSalivaProteins.fasta";
    private static String datoutputfile = "/data/2/taoxu/projects/xuemei/database/identifiedSalivaUniProt.dat";
    private static String salivaIpis = "/data/2/taoxu/projects/xuemei/database/map.txt";

    public static void main(String args[]) throws IOException {
     
        getIdentifiedIpis(); 
        //for (Iterator<Fasta> itr = FastaReader.getFastas(new FileInputStream(humanFasta)); itr.hasNext(); ) { 
        //    Fasta f = itr.next();
       // }
        int numAdded = 0;

        PrintStream fastafile = new PrintStream(fastaoutputfile); 
        PrintStream datfile = new PrintStream(datoutputfile); 

        for (Iterator<UniProtProtein> itr = UniprotDatReader.getUniProtProteins(new FileInputStream(uniprotdatfile)); itr.hasNext(); ) { 
            UniProtProtein protein = itr.next();
            if(isIdentified(protein)) {
                fastafile.println(protein.outputFasta());
                datfile.println(protein.getOutput());
                numAdded++;
            }           
        }
        System.out.println("NumAdded: " + numAdded);
    }
    public static boolean isIdentified(UniProtProtein p) {
        for(Iterator<String> it = p.getAcs().iterator(); it.hasNext();) {
            String s = it.next();
            if(identifiedProteinIds.contains(s)) {
                //System.out.println(s);
                return true;
            }

        }
        return false;
    }
    public static void getIdentifiedIpis() throws IOException {
        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(salivaIpis)));
        br.readLine(); // remove the first line
        String line = null;
        while ((line = br.readLine()) != null) {
            String [] arr = line.split("\\|"); 
            //System.out.println(line + "\t");
            for(int i = 1; i < arr.length; i++) {
                String temp = arr[i];
                //System.out.println(temp);
                String [] ids = temp.split(":")[1].split(";");
                for(String s : ids) {
                    String myid = s.split("-")[0];
                    identifiedProteinIds.add(myid);
                    //System.out.println(myid);
                }
            }
        }
//System.out.println("Number of ids: " + identifiedProteinIds.size());
    }
}
