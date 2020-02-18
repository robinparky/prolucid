import java.util.*;
import java.io.*;
import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;
/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Tao xu 
 * @version $Id: PlasmaTester.java,v 1.2 2005/10/17 16:45:07 taoxu Exp $
 */
public class PlasmaTester
{
    private String fastaFileName;
    private HashMap<String, String> refseq2Ipi = new HashMap<String, String>();
    private HashMap<String, String> ipi2Name = new HashMap<String, String>();
    private HashSet<String> plasmaAccs = new HashSet<String>();
    private HashSet<String> refseqs = new HashSet<String>();
    private static final String USAGE = "java PlasmaTester IPIFastaFile PlasmaAccFile RefseqAccFile";
    public static void main(String args[]) throws IOException
    {
        PlasmaTester pt = new PlasmaTester(args[0], args[1], args[2]);
        pt.outputResult();
    }

    public PlasmaTester(String ipiFastaFile, String plasmaFile, String refseqFile) throws IOException
    {
        //fastaBr = new BufferedReader(new FileReader(plasmaFile));
        fastaFileName = ipiFastaFile;
        loadHashMap();
        plasmaAccs = readAccs(plasmaFile);
        refseqs = readAccs(refseqFile);        
    }
    public void outputResult() throws IOException {
        PrintStream plasma = new PrintStream("plasma.txt");
        PrintStream nonPlasma = new PrintStream("nonPlasma.txt");
        PrintStream unknown = new PrintStream("unknown.txt");
        for(Iterator<String> it = refseqs.iterator(); it.hasNext();) {
            String refseqAcc = it.next();
            String ipi = refseq2Ipi.get(refseqAcc);
            if(ipi != null) {
                String name = ipi2Name.get(ipi);
                if(plasmaAccs.contains(ipi)) {
                    plasma.println(refseqAcc + "\t" + ipi + "\t" + name);
                } else {

                    nonPlasma.println(refseqAcc + "\t" + ipi + "\t" + name);
                }
            } else {
                unknown.println(refseqAcc);
            }
        }        
    } 
    private HashSet<String> readAccs(String accFileName) throws IOException {
        HashSet<String> accs = new HashSet<String>();
        BufferedReader br = new BufferedReader(new InputStreamReader
                              (new FileInputStream(accFileName)));
        String line = null;
        while ((line = br.readLine()) != null) {
//        System.out.println("before remove version: " + line);
             String acc = removeVersion(line.trim());
//        System.out.println("after remove version: " + acc);
             accs.add(acc);
        }
        br.close();
        return accs;
    }
    public static String removeVersion(String acc) {
        
//        System.out.println("before remove version: " + acc);
        if(acc.indexOf(".") > 0) {
            acc = acc.substring(0, acc.indexOf("."));
        }
//        System.out.println("after remove version: " + acc);
        return acc;
    }
    private void loadHashMap() throws IOException {
        FileInputStream fis = new FileInputStream(new File(fastaFileName));
        for (Iterator<Fasta> it = FastaReader.getFastas(fastaFileName); it.hasNext();) {
            Fasta f = it.next();
            parseDefline(f.getDefline());
        }
    } 
    private void parseDefline(String defline) {
        int index = defline.indexOf(" ");
        String accPart = defline.substring(0, index);
        String namePart = defline.substring(index+1);
        String [] def = accPart.substring(1).split("\\|");
        String ipiAcc = removeVersion(getAcc(def[0]));
        ipi2Name.put(ipiAcc, getName(namePart)); 
        for(int i = 1; i < def.length; i++) {
            if(def[i].startsWith("REFSEQ")) {
                String [] accs = getAccs(def[i]); 
                if (accs.length > 1) {
                    System.out.print("Number of REFSEQ accs: " + accs.length);
                    for(String s : accs) {
                        System.out.print("\t" + s);
                    }
                    System.out.println();
                }
                for(String acc : accs) {
                    refseq2Ipi.put(removeVersion(acc), ipiAcc);
                }
            }
        }
    }
    private String getName(String description) {
        int index = description.indexOf(" ") + 1;
//System.out.println(description + "\t" + index + "\t" + description.substring(index));
        return description.substring(index);
    }
    private String [] getAccs(String defPart) {
        int index = defPart.indexOf(":") + 1;
        return defPart.substring(index).split(";");
    }
    private String getAcc(String defPart) {
        int index = defPart.indexOf(":");
        return defPart.substring(index + 1);
    }

}
