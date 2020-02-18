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
 * @version $Id: FastaResidueFrequencyCalculator.java,v 1.1 2008/10/24 22:39:14 taoxu Exp $
 */
public class FastaResidueFrequencyCalculator
{
    private static  String fastaFileName;
    private static final String USAGE = "java FastaResidueFrequencyCalculator IPIFastaFile";
    public static void main(String args[]) throws IOException {
        int [] freq = new int[256];
        char [] residues = {'G','A','S','P','V','T','C','L','I','N','D','Q','K','E','M','H','F','R','Y','W' };
        try {
            fastaFileName = args[0];
            FileInputStream fis = new FileInputStream(new File(fastaFileName));
            for (Iterator<Fasta> it = FastaReader.getFastas(fastaFileName); it.hasNext();) {
                Fasta f = it.next();
                String seq = f.getSequence();
                int length = seq.length();
                for(int i = 0; i < length; i++) {
                    freq[seq.charAt(i)]++;
                }
            }
            int total = 0;
            for(char c : residues) {
                total += freq[c];
            }
            for(char c : residues) {
                System.out.println(c + "\t" + freq[c] + "\t" + (0.0+ freq[c])/total);
            }
        } catch (Exception e) {
            System.err.println(USAGE);
            System.exit(1);
        }
    }

}
