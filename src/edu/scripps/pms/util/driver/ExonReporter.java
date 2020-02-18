import java.util.*;
import java.io.*;
import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.dtaselect.Protein;
import edu.scripps.pms.util.dtaselect.Peptide;
import edu.scripps.pms.util.io.DTASelectFilterReader;
/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2011</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Tao xu 
 * @version $Id
 */
public class ExonReporter
{
    //private static string exonindexmapfile = "/data/2/evelyn/project/exon/NRXN1_X_exonindexmap.txt";
    //private static string fastafile = "/data/2/evelyn/project/exon/output_ExonCombinationsPlusRestricted.fasta";
    //private static string dtaselectresult = "/data/3/taoxu/projects/jsavas/exon/nrxn1x/DTASelect-filter.txt";
    // display two tables, one indicates which exon joined with the other, raw and column indicate the exons 
    // one display which exons is present
    // We should not allow any false positives in the DTASelect result, so better to use orbitrap data
    private String fastaFileName;
    private String exonMapFile; // format of exonMapFile is exonId\tEnsembleID\tLength
    private String dtaselectFile; // DTASelect-filter.txt file
    private int[][] joindExons;   
    private int [] exonFreq;   
    private int [] exonLengths = new int[1000]; // assume no more than 1000 exons per gene
    private String [] exonEnsembleIds = new String [1000];  // ensemble ID of each exon

    private HashMap<String, Peptide> seq2Peptide = new HashMap<String, Peptide>(1000000); // peptide seq to peptide objext
    private HashMap<String, Fasta> acc2Fasta = new HashMap<String, Fasta>(100000); // 
    int maxExonId = 0; 

 
    private static final String USAGE = "java ExonReporter DTASelect-filter_file exonidenxm_map_file fastadatabase_file";

    public static void main(String args[]) throws IOException
    {
        String currentdir = System.getProperty("user.dir");
        System.out.println("H\tWorking in " + currentdir);
        System.out.println("H\tDTASelect-filter.txt file: " + args[0]);
        System.out.println("H\tExon Index File: " + args[1]);
        System.out.println("H\tFasta Database File: " + args[2] + "\n");

        ExonReporter er = new ExonReporter(args[0], args[1], args[2]);
        er.outputResult();
    }

    public ExonReporter(String dtaselectfile, String exonmapfile, String fastafile) throws IOException
    {
        //fastaBr = new BufferedReader(new FileReader(plasmaFile));
        fastaFileName = fastafile;
        
        dtaselectFile = dtaselectfile ;
        exonMapFile = exonmapfile;        
        readExonMapFile(); 
        loadHashMap();
    }
    public void outputResult() throws IOException {
        //PrintStream plasma = new PrintStream("plasma.txt");

        System.out.println();
        System.out.println("Exon Length");
        StringBuffer header = new StringBuffer();
        StringBuffer content = new StringBuffer();
        for(int i = 1; i < exonFreq.length; i++) {
            header.append(i + "/" + exonEnsembleIds[i]);
            header.append("\t");
            content.append(exonLengths[i]);
            content.append("\t");
        }
        System.out.println(header);
        System.out.println(content);
        System.out.println();
        

        System.out.println("Exon Frequence");
        header = new StringBuffer();
        content = new StringBuffer();

        int numexonsidenfied = 0;
        StringBuffer exonsb = new StringBuffer();
        for(int i = 1; i < exonFreq.length; i++) {
            header.append(i + "/" + exonEnsembleIds[i]);
            header.append("\t");
            content.append(exonFreq[i]);

            content.append("\t");
            if(exonFreq[i] > 0) {
                numexonsidenfied++;
                exonsb.append("\t");
                exonsb.append(i);

            }
        }
        System.out.println(header);
        System.out.println(content);
        System.out.println();
        System.out.println("\nTotal Number of Exons Identified:\t" + numexonsidenfied + "\tidentified_exon_ids" + exonsb.toString() + "\n");
      
      
        System.out.println("Joined Exon Frequence");
        //header = new StringBuffer();
        content = new StringBuffer();
        for(int i = 1; i < exonFreq.length; i++) {
            content.append(i + "/" + exonEnsembleIds[i] + "\t");
            for(int j = 1; j < exonFreq.length; j++) {
                content.append(joindExons[i][j]);       
                content.append("\t");
                if(joindExons[i][j] > 0) {
                    System.out.println("Joined Exons Found " + i + "\t" + j + " with freq " + joindExons[i][j]);
                }
            }
            content.append("\n");
        }
         
        System.out.println();
        System.out.println("Exons\t" + header);
        System.out.println(content);
    } 
    private void readExonMapFile() throws IOException {
        BufferedReader br = new BufferedReader(new InputStreamReader
                              (new FileInputStream(exonMapFile)));
        String line = null;
        while ((line = br.readLine()) != null) {
            String [] arr = line.split("\t");
            if(arr.length > 2) {
                int exonid = Integer.parseInt(arr[0].trim()); 
                exonLengths[exonid] = Integer.parseInt(arr[2].trim()); 
                exonEnsembleIds[exonid] = arr[1]; 
                maxExonId = maxExonId > exonid? maxExonId : exonid; 
            }
        }
        br.close();
    }
    private void processPeptide(Fasta f, Peptide p) {
        String seq = p.getMidSeq();
        String acc = f.getAccession();
        int index = f.getSequence().indexOf(seq);
        String arr[] = acc.split("_");
        int frame = Integer.parseInt(arr[arr.length-1]);
        int lengths[] = new int[arr.length - 5 ]; //e.g. NRXN1_X_11_31_40_Frame_Plus_1
        int sumlength[] = new int[arr.length - 5 ]; //e.g. NRXN1_X_11_31_40_Frame_Plus_1
        int exonids [] = new int[arr.length - 5 ]; //e.g. NRXN1_X_11_31_40_Frame_Plus_1
        int sumlengthtemp = 0;

        for(int i = 2; i < arr.length -3; i++) {
            int exonid = Integer.parseInt(arr[i]); 
            exonids[i-2] = exonid;
            lengths[i-2] = exonLengths[exonid];
            sumlengthtemp += exonLengths[exonid]; 
            sumlength[i-2] = sumlengthtemp;
        }
        int start = index * 3 + frame;
        //int stop = (index + seq.length()) * 3 + frame;
        int stop = start + seq.length()*3;
        for(int i = 0; i < sumlength.length; i++) {
            if(stop <= sumlength[i]) {
                if(i > 0) {
                    if(start > sumlength[i-1]) {
                        int exonid = exonids[i];
                        exonFreq[exonid]++;
                    } else {
                        int exonid1 = exonids[i-1];
                        int exonid2 = exonids[i];
                        exonFreq[exonid1]++;
                        exonFreq[exonid2]++;
                        joindExons[exonid1][exonid2]++;
        System.out.println("two exon join: " + acc + "\t" + seq + "\t joined " + exonid1 + " and " + exonid2 + " index:" + index +  "\tand nucleotide position\t" + start + " stop: " + stop +  "\ti: " + i + "\tsumlength[i]: " + sumlength[i] + "\tsumlength[i-1]: " + sumlength[i-1]); 
                        int twobefore = i - 2;
                        if(twobefore > -1 && start <= sumlength[twobefore]) {
                            int exonid3 = exonids[twobefore];
                            exonFreq[exonid3]++;
                            joindExons[exonid3][exonid1]++;
        System.out.println("three exon join: " + acc + "\t" + seq + "\t joined " + exonid3 + " and " + exonid1 + " index:" + index +  "\tand nucleotide position\t" + start + " stop: " + stop +  "\ti: " + i + "\tsumlength[i]: " + sumlength[i] + "\tsumlength[i-1]: " + sumlength[i-1] + "\tsumlenth[i-2]: " + sumlength[i-2]); 
                        }
                    } 
                } else {
                    int exonid = exonids[i];  // first exon
                    exonFreq[exonid]++; 
                }
                break;
            }
        }
    }
    private void loadHashMap() throws IOException {
 
        HashSet<String> identifiedprotienaccs = new HashSet<String>(100000);
        DTASelectFilterReader dtareader = new DTASelectFilterReader(dtaselectFile);    
        ArrayList<Protein> proteins = dtareader.getProteinList();
        dtareader.close();
        for(Iterator<Protein> it = proteins.iterator(); it.hasNext();) {
            String acc = it.next().getAccession();
            identifiedprotienaccs.add(acc); 
//System.out.println("acc identified: " + acc);
        }

        FileInputStream fis = new FileInputStream(new File(fastaFileName));
        for (Iterator<Fasta> it = FastaReader.getFastas(fastaFileName); it.hasNext();) {
            Fasta f = it.next();
            String acc = f.getAccession();
            if(identifiedprotienaccs.contains(acc)) {
                acc2Fasta.put(acc, f);                
  
//System.out.println("acc to put: " + acc);
            } else {
//System.out.println("acc not identified: " + acc);
            }
        }
        fis.close();
        System.out.println("number of accs added: " + acc2Fasta.size() + "\tnumber of proteins identified: " +  proteins.size());
                
        HashSet<String> addedpeptides = new HashSet<String>(1000000);
        joindExons = new int[maxExonId + 1][maxExonId + 1];
        exonFreq = new int[maxExonId + 1];
        for(Iterator<Protein> it = proteins.iterator(); it.hasNext();) {
            Protein prot = it.next();
            String acc = prot.getAccession();
//System.out.println("acc to retrieve: " + acc);
            Fasta f = acc2Fasta.get(acc);
            if(f == null) {

                System.out.println("acc to retrieve but not found: " + acc); // this should never happen
            }
            for(Iterator<Peptide> pepit = prot.getPeptides(); pepit.hasNext();) {
                Peptide p = pepit.next();
                String seq = p.getMidSeq();
                if(!addedpeptides.contains(seq)) {
                    processPeptide(f, p); 
                    addedpeptides.add(seq); 
                }
            } 
        }
        

    } 

}
