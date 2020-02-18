
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.mspid.ProteinDatabase;
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
import java.util.HashMap;


/**
 * @author  Tao Xu
 * @version $Id
 */
public class IPIId2UniProtCellularLocation {
    private static HashSet<String> identifiedProteinIds = new HashSet<String>(1000000); // swissprot of trembl id
    private static String uniprotdatfile = "/home/taoxu/taoxu_on_data/projects/lujian/db/uniprot_sprotplustrembl_rodents.dat";
    private static String fastadb = "/bluefish/people-b/applications/yates/dbase/EBI-IPI_rat_3.30_06-28-2007_con_reversed.fasta";
    private static String identifiedIpis = "/home/taoxu/taoxu_on_data/projects/lujian/db/mitIdentifiedIpi.txt";
    private static HashMap<String, String> id2Ipi = new HashMap<String, String>(10000);

    public static void main(String args[]) throws IOException {
     
        getIdentifiedIds(); 
        //for (Iterator<Fasta> itr = FastaReader.getFastas(new FileInputStream(fastadb)); itr.hasNext(); ) { 
        //    Fasta f = itr.next();
       // }
        int numAdded = 0;

        


        for (Iterator<UniProtProtein> itr = UniprotDatReader.getUniProtProteins(new FileInputStream(uniprotdatfile)); itr.hasNext(); ) { 
            UniProtProtein protein = itr.next();
            String ipi = isIdentified(protein);
            if(ipi != null) {
                String location = protein.getSubcellularLocation();
                if(location != null) {
                    System.out.println(ipi + "\t" + location);
                }
                numAdded++;
            }           
        }
        System.out.println("NumAdded: " + numAdded);
    }

    // return the IPI ac for this protein, null if not identified
    public static String isIdentified(UniProtProtein p) {
        for(Iterator<String> it = p.getAcs().iterator(); it.hasNext();) {
            String s = it.next();
            if(identifiedProteinIds.contains(s)) {
                //System.out.println(s);
                return id2Ipi.get(s);
            }

        }
        return null;
    }
    public static void getIdentifiedIds() throws IOException {

        ProteinDatabase pdb = new ProteinDatabase(fastadb);        
        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(identifiedIpis)));
        br.readLine(); // remove the first line
        String line = null;
        while ((line = br.readLine()) != null) {
            String ipiac = line.trim();
            String defline = pdb.accession2Fasta(ipiac).getDefline();
          
            String [] arr = defline.split("\\|"); 
            //System.out.println(line + "\t");
            for(int i = 1; i < arr.length; i++) {
                String temp = arr[i];
                //System.out.println(temp);
                String [] ids = temp.split(":")[1].split(";");
                for(String s : ids) {
                    String myid = s.split("-")[0];
                    identifiedProteinIds.add(myid);
                    id2Ipi.put(myid, ipiac);
                    //System.out.println(myid);
                }
            }
        }
//System.out.println("Number of ids: " + identifiedProteinIds.size());
    }
}
