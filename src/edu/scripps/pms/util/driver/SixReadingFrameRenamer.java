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
 * @version $Id: SixReadingFrameRenamer.java,v 1.1 2009/01/29 19:07:02 taoxu Exp $
 */
public class SixReadingFrameRenamer
{
    private static  String fastaFileName;
    private static final String USAGE = "java SixReadingFrameRenamer FastaFile";
    public static void main(String args[]) throws IOException {
        try {
            fastaFileName = args[0];
            FileInputStream fis = new FileInputStream(new File(fastaFileName));
            for (Iterator<Fasta> it = FastaReader.getFastas(fastaFileName); it.hasNext();) {
                Fasta f = it.next();
                String seq = f.getSequence();
                String defline = f.getDefline();
                String acc = f.getAccession();
                int index = defline.length();
                String frame = defline.substring(index-2, index);
                //System.out.println(acc + "\t" + frame);
                String newdefline = ">" + frame + acc + " " + defline.replace("|Trans", "| Trans");
                System.out.println(newdefline); 
                System.out.println(seq); 
            }
            
        } catch (Exception e) {
            System.err.println(USAGE);
            System.exit(1);
        }
    }

}
