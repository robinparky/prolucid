import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import java.io.*;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.io.DTASelectFilterReader;
import edu.scripps.pms.util.dtaselect.Peptide;
import edu.scripps.pms.util.dtaselect.Protein;


/**
 * <p>Program to six-frame translate a nucleotide sequence</p>
 */

// It now ignores frames with more than 2 stop condons
public class JamieSeqCov {


    public static final String database = "/lustre/people/applications/yates/dbase/SGD_S-cerevisiae_na_12-16-2005_con_reversed.fasta";
    //public static final String interestedacsfile = "/data/3/taoxu/projects/catherine/jamie/spblist.txt";
    public static final String interestedacsfile = "/data/3/taoxu/projects/catherine/jamie/menlist.txt";
    //public static final String dtaselectfiles = "/data/2/rpark/ip2_data/catclw/Jamie/allfiles.txt";
    public static final String dtaselectfiles = "/data/2/rpark/ip2_data/catclw/Jamie/allfiles_no_tub4_invitro.txt";

    /**
     * Call this to get usage info, program terminates after call.
     */
    public static void help() {
	System.out.println(
		"usage: java JamieSeqCov ");
	System.exit( -1);
    }

    public static int getPeptideLocation(String peptide, Fasta protein) {
        return protein.getSequence().indexOf(peptide);
    }

    public static void main(String[] args) throws Exception{

	BufferedReader br = null;
        HashSet<String> interestedacs = new HashSet();
        HashMap<String, String> ac2fasta = new HashMap();
        HashMap<String, int[]> ac2seqcov = new HashMap();

	try {
	    br = new BufferedReader(new FileReader(interestedacsfile));
            String line = br.readLine();
             
            while(line != null) {
                line = line.trim();
                if(!"".equals(line)) {
                    interestedacs.add(line);
                } 
                line = br.readLine();
            } 
            br.close();

          
           
            FileInputStream fis = new FileInputStream(new File(database));
            Iterator<Fasta> fastas = FastaReader.getFastas(fis);
            while(fastas.hasNext()) {
                Fasta f = fastas.next();
                String myac = f.getAccession();
                if(interestedacs.contains(myac) ) {
                    String seq = f.getSequence();
                    ac2fasta.put(myac, seq);
                    ac2seqcov.put(myac, new int[seq.length()]);
                }
            }
 
	    br = new BufferedReader(new FileReader(dtaselectfiles));
            line = br.readLine();
            while(line != null) {
                line = line.trim();
                if(!"".equals(line)) {
                    System.out.println("Now processing " + line);
                    DTASelectFilterReader reader = new DTASelectFilterReader(line);
                    for (Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {

                        Protein p = itr.next();
                        String proteinid = p.getAccession();
                        
                        if(interestedacs.contains(proteinid)) {
                            int [] seqcov = ac2seqcov.get(proteinid);
                            String protseq = ac2fasta.get(proteinid);
                            for(Iterator<Peptide> pepItr=p.getPeptides(); pepItr.hasNext(); ) {
                                Peptide peptide = pepItr.next();
                                String pepseq = peptide.getSeqWithNoModification();
                                int length = pepseq.length();
                                int pos = protseq.indexOf(pepseq);
                                if(pos >= 0) {
                                    for(int i = pos; i < (pos + length); i++) {
                                        seqcov[i]++;
                                    }
                                } else {
                                    
                                    System.out.println("Something wrong with " + pepseq);
                                }

                            }
                        }


                    }
                } 
                line = br.readLine();
            } 
            br.close();

            System.out.println("\n\n\n\nProtein Accession \tTotal Sequence Coverage");
            for(Iterator<String> it = interestedacs.iterator(); it.hasNext();) { 
                String ac = it.next();
                int [] seqcov = ac2seqcov.get(ac);
                int numzero = 0;
                for(int i = 0; i < seqcov.length; i++) {
                    if(seqcov[i] == 0) {
                        numzero++;
                    } 
                }
                int length = seqcov.length;
                //System.out.println(ac + "\t" + length + "\t" + numzero + "\t" + (length - numzero)/(length+0.0)*100 + "%" );
                System.out.println(ac + "\t" + (length - numzero)/(length+0.0)*100 + "%" );
            }
            

	} finally {
	    //tidy up
	    if(br != null){
		br.close();
	    }
	}
    }
}

