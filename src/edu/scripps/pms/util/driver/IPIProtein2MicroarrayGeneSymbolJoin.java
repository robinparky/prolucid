
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
public class IPIProtein2MicroarrayGeneSymbolJoin {
    public static String USAGE = "Usage: java -Xmx500M IPIProtein2MicroarrayGeneSymbolJoin protein_file microarray_file  fasta_file";
    private static HashSet<String> ipiacs = new HashSet<String>(1000000); // swissprot of trembl id
    private static HashMap<String, String> ipiac2proteinentry = new HashMap<String, String>();
    private static String fastadb = "/lustre/people/applications/yates/dbase/EBI-IPI_mouse_3.52_11-21-2008_reversed.fasta"; 
    private static String ipiFileName = "/data/8/taoxu_on_data8/changchun/20110511/protein-expression.txt"; // protein files name (census output)
    private static String microarrayfile = "/data/8/taoxu_on_data8/changchun/20110511/mRNA-expression.txt";
    private static HashMap<String, HashSet<String>> ipiac2genesymbol = new HashMap();
    private static HashMap<String, HashSet<String>> genesymbol2ipiac = new HashMap();
    private static HashMap<String, String> ipiac2genedescription = new HashMap<String, String>();
    private static HashMap<String, String> genedescription2ipiac = new HashMap<String, String>(); // actually should be proteindescription


    public static void main(String args[]) throws IOException {
        try { 
            if(args.length == 3) {
                ipiFileName = args[0];
                microarrayfile = args[1];
                fastadb = args[2];
            }
    
            readProteinIds(); 
            int numgenesymbol = 0;
            int numipiac = 0;
            for (Iterator<Fasta> itr = FastaReader.getFastas(new FileInputStream(fastadb)); itr.hasNext();) { 
                Fasta f = itr.next();
                String ipiac = f.getAccession();
                if(ipiacs.contains(ipiac)) {
                    String defline = f.getDefline();
                    String [] arr = defline.split("Gene_Symbol="); 
                    if(arr.length > 1) {
                        String [] arr2 = arr[1].split(" ");
                        String genesymbol = arr2[0];

                         HashSet<String> symbols = ipiac2genesymbol.get(ipiac);
                         if(symbols == null) {
                             symbols = new HashSet();
                             ipiac2genesymbol.put(ipiac, symbols);
                         }

                        //System.out.println(genesymbol + " is Gene_Symbol for " + ipiac + " is " + genesymbol);
                        String [] arr3 = genesymbol.split(";");
                        for(int i = 0; i < arr3.length; i++) {
                            String sym = arr3[i].toUpperCase();
                            if(sym == null || "".equals(sym) || "-".equals(sym) || "_".equals(sym)) {

                            } else {
                                HashSet<String> acs = genesymbol2ipiac.get(sym);
                                if(acs == null) {
                                    acs = new HashSet();
                                    genesymbol2ipiac.put(sym, acs);
                                }
                                acs.add(ipiac);
//if("IK".equals(sym))  System.out.println(sym + "\t" + ipiac);
                                symbols.add(sym);
                            }
                        }

                        StringBuffer sb = new StringBuffer();
                        for(int i = 1; i < arr2.length; i++) {
                            String s = arr2[i];
                            if(s.startsWith("similar") || s.startsWith("Similar")) {
                                i++; i++; // ignore similar to
                            }
                            sb.append(arr2[i]);
                            sb.append(" ");
                        }
//System.out.println(ipiac + "\t" + sb.toString());
                        String pdesc = sb.toString().trim();
                        if("mRNA".equals(pdesc) || "Protein".equals(pdesc) || "Predicted".equals(pdesc) || "MRNA,".equals(pdesc) ||
                           pdesc.indexOf("kDa protein") != -1 || pdesc.indexOf("hypothetical protein") != -1 ||
                           pdesc.indexOf("Putative uncharacterized protein") != -1 || pdesc.indexOf("Novel protein") != -1) 
                        {
//System.out.println(pdesc + "\tfrom tostring is " + sb.toString());
                           // ignore these special cases
                        } else {
//if(pdesc.indexOf("mRNA") != -1 || pdesc.indexOf("MRNA") != -1) System.out.println(pdesc + "\tfrom tostring is " + sb.toString());
                           String desc1 = sb.toString().trim().toUpperCase();
                           ipiac2genedescription.put(ipiac, desc1);
                           genedescription2ipiac.put(desc1, ipiac);
                        }
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
        HashSet<String> positiveproteins = new HashSet(1000000);      
 
        String prefix = microarrayfile.substring(0, microarrayfile.lastIndexOf("."));
        PrintStream ps = new PrintStream(prefix + "_joined.txt");

        HashSet<String> proteinsymbols = new HashSet(1000000);
        HashSet<String> proteindescriptions = new HashSet(1000000);

        Iterator<HashSet<String>> symbolit = ipiac2genesymbol.values().iterator();
        while(symbolit.hasNext()) {
            Iterator<String> it = symbolit.next().iterator();
            while(it.hasNext()) proteinsymbols.add(it.next());
        }
        Iterator<String> it = ipiac2genedescription.values().iterator(); 
        while(it.hasNext()) proteindescriptions.add(it.next()); 
    

        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(microarrayfile)));
        String line = null;
        int num = 0;
        while ((line = br.readLine()) != null) {
            String [] arr = line.split("\t");
       //System.out.print("Now processing line " + ++num + "\r"); 
            String symbolfromarr = arr[8].toUpperCase().trim(); 
            String desc = arr[0];

            String symb = null;
            if(proteinsymbols.contains(symbolfromarr)) {
                symb = symbolfromarr;
            } else { 
                HashSet<String> prtdescs = geneSymbol2ProteinDescription(proteindescriptions, symbolfromarr); 
                if(prtdescs.size() > 0) {
                    // gene symbol is found in at least of one of the protein descriptions
                    // need to handle multiple
                    symb = prtdescs.iterator().next();
                } else {
                    // need to handle multiples here too
                    symb = getGeneSymbolFromDescription(proteinsymbols, proteindescriptions, desc);
                }
            }
            System.out.println(symb + "\t" + symbolfromarr);
            StringBuffer sb = new StringBuffer();
            
            sb.append(line);
            sb.append("\t" + symb);
            if(symb != null) {
                HashSet<String> acs = genesymbol2ipiac.get(symb); 
                if(acs != null) {
                    for(Iterator<String> ita = acs.iterator(); ita.hasNext();) {
                        String ac = ita.next();
                        sb.append("\t"+ ipiac2proteinentry.get(ac));    
                        positiveproteins.add(ac);
                    }
                } else {
                    String ac = genedescription2ipiac.get(symb);
                    sb.append("\t" + ipiac2proteinentry.get(ac));
                    positiveproteins.add(ac);
                }
            } else {
                sb.append("\tNA");    
            } 
            
            ps.println(sb.toString());
        }

        int numfound = 0;
        int numnotfound = 0;        
        for(Iterator<String> itno = ipiacs.iterator(); itno.hasNext();) {
            String ac = itno.next();
            if(positiveproteins.contains(ac)) {
                numfound++;
            } else {
                numnotfound++;
 System.out.println("Unmatched protein descriptin " + ipiac2genedescription.get(ac) + "\t" + ipiac2proteinentry.get(ac));               
                ps.println("Unmatched proteins\t" + ipiac2proteinentry.get(ac));
            }
        } 
        br.close();
        ps.close();
        System.out.println("Total number of proteins: " + ipiacs.size());
        System.out.println("Number of proteins matched: " + numfound);
        System.out.println("Number of proteins unmatched: " + numnotfound);
    }

    

    private static String getGeneSymbolFromProteinDescription() { 
        return null;
    }

    // return the protein descriptions that contains the genesymbol
    private static HashSet<String> geneSymbol2ProteinDescription(HashSet<String> proteindescs, String genesymbol) { 
        HashSet<String> descs = new HashSet();
        for(Iterator<String> it = proteindescs.iterator(); it.hasNext();) {
            String protdesc = it.next().toUpperCase();
            if(protdesc.indexOf(genesymbol) != -1 ) {
                descs.add(protdesc);
            }
        }
        return descs;
    }
    // symbol is gene_symbol from microarray file and desc is first column in microarray
    private static String getGeneSymbolFromDescription(HashSet<String> proteinsymbols, HashSet<String> proteindescs, String desc) {
        desc = desc.toUpperCase();
        for(Iterator<String> it = proteinsymbols.iterator(); it.hasNext();) {
            String sym = it.next();
            //desc.toUpperCase();
            //if(desc.equals(sym)) {
            if(desc.indexOf(" " + sym + " ") != -1 || desc.indexOf("(" + sym + ")") != -1)  {
                //if("IK".equals(sym)) System.out.println(sym + " form proteinsymbol");
                return sym;
            }
        }
        
        for(Iterator<String> it = proteindescs.iterator(); it.hasNext();) {
            String sym = it.next();
            if(desc.indexOf(sym) != -1) {
                //if("IK".equals(sym)) System.out.println(sym + " form proteindesc");
                return sym;
            }
        }

        System.out.println("null is returned for " + desc);
        return null;
    }

    public static void readProteinIds() throws IOException {


        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(ipiFileName)));
        //br.readLine(); // remove the first line
        String line = null;
        while ((line = br.readLine()) != null) {
            String [] arr = line.trim().split("\t");
            ipiacs.add(arr[0]);
            ipiac2proteinentry.put(arr[0], line);
        }
        br.close();
//System.out.println("Number of ids: " + ipiacs.size());
    }
}
