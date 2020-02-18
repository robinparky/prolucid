
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
public class IPIId2GeneSymbol {
    public static String USAGE = "Usage: java -Xmx500M IPIId2GeneSymbol IPI_accession_file fasta_file";
    private static HashSet<String> ipiacs = new HashSet<String>(1000000); // swissprot of trembl id
    private static String fastadb = "/lustre/people/applications/yates/dbase/EBI-IPI_mouse_3.52_11-21-2008_reversed.fasta"; 
    private static String ipiFileName;
    private static HashMap<String, String> ipiac2genesymbol = new HashMap<String, String>();
    private static HashMap<String, String> ipiac2genedescription = new HashMap<String, String>();

    public static void main(String args[]) throws IOException {
        try { 
            ipiFileName = args[0];
            if(args.length > 1) {
                fastadb = args[1];
            }
    
            getIdentifiedIds(); 
            int numgenesymbol = 0;
            int numipiac = 0;
            for (Iterator<Fasta> itr = FastaReader.getFastas(new FileInputStream(fastadb)); itr.hasNext();) { 
                Fasta f = itr.next();
                String ipiac = f.getAccession();
                
                String ipiacnoversion = ipiac.split("\\.")[0];
//System.out.println("ipiac: " + ipiac + "\tipiacnoversion: " + ipiacnoversion);
                if(ipiacs.contains(ipiac) || ipiacs.contains(ipiacnoversion)) {
//System.out.println("found ipiac: " + ipiac + "\tipiacnoversion: " + ipiacnoversion);
                    String defline = f.getDefline();
                    String [] arr = defline.split("Gene_Symbol="); 
//System.out.println("arr.lenth: " + arr.length + "\tipiacnoversion: " + ipiacnoversion + "\tdefline: " + defline);
                    if(arr.length > 1) {
                        String [] arr2 = arr[1].split(" ");
                        String genesymbol = arr2[0];

                        //System.out.println("Gene_Symbol for " + ipiac + " is " + genesymbol);
                        ipiac2genesymbol.put(ipiac, genesymbol);
                        ipiac2genesymbol.put(ipiacnoversion, genesymbol);
                        StringBuffer sb = new StringBuffer();
                        for(int i = 1; i < arr2.length; i++) {
                            String s = arr2[i];
                            if(s.startsWith("similar") || s.startsWith("Similar")) {
                                i++; i++; // ignore similar to
                            }
                            sb.append(arr2[i]);
                            sb.append(" ");
                        }
                        ipiac2genedescription.put(ipiac, sb.toString());
                        ipiac2genedescription.put(ipiacnoversion, sb.toString());
                    }
            
                }
            }

            output();        
       } catch(Exception e) {
           System.err.println(USAGE);
           e.printStackTrace();
       }

      // System.out.println("NumAdded: " + numAdded);
    }

    public static void output() throws IOException {
       
        String prefix = ipiFileName.substring(0, ipiFileName.lastIndexOf("."));
        PrintStream ps = new PrintStream(prefix + "_with_gene_symbol.txt");
        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(ipiFileName)));
        String line = null;
        while ((line = br.readLine()) != null) {
            String ipiac = line.trim();
//System.out.println("line: " + ipiac + "\tipiac2genesymbol.szie: " + ipiac2genesymbol.size());
            String genesymbol = ipiac2genesymbol.get(ipiac);
            String description = ipiac2genedescription.get(ipiac); 
//System.out.println("genesymbol: " + genesymbol + "\tdesc: " + description);
            if(genesymbol != null) {
//System.out.println("Found genesymbol: " + genesymbol + "\tdesc: " + description);
                String [] arr = genesymbol.split(";");
                
                for(int i = 0; i < arr.length; i++) {
                    ps.println(ipiac + "\t" + arr[i] + "\t" + description);
                }
            } else {
                ps.println(ipiac + "\t" + "" + "\t" + description);
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
