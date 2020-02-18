import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import java.io.*;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.io.DTASelectFilterReader;
import edu.scripps.pms.util.dtaselect.Peptide;
import edu.scripps.pms.util.dtaselect.Protein;
import edu.scripps.pms.util.dtaselect.ProteinGroup;


/**
 * <p>Program to six-frame translate a nucleotide sequence</p>
 */

// It now ignores frames with more than 2 stop condons
public class TomMod {


    public static final String database = "/data/2/rpark/ip2_data/taoxu/database/DrosMelano7227_reversed.fasta";
    //public static final String interestedacsfile = "/data/3/taoxu/projects/catherine/jamie/spblist.txt";
    public static final String interestedacsfile = "/data/3/taoxu/projects/james/tom/20120310/dynein_subunits.txt";
    //public static final String dtaselectfiles = "/data/2/rpark/ip2_data/catclw/Jamie/allfiles.txt";
    public static final String dtaselectfiles = "/data/3/taoxu/projects/james/tom/20120310/3modfolders.txt";

    /**
     * Call this to get usage info, program terminates after call.
     */
    public static void help() {
	System.out.println(
		"usage: java TomMod ");
	System.exit( -1);
    }

    public static int getPeptideLocation(String peptide, Fasta protein) {
        return protein.getSequence().indexOf(peptide);
    }

    public static void main(String[] args) throws Exception{

        String modstring = "79";
	BufferedReader br = null;
        HashSet<String> interestedacs = new HashSet();
        HashMap<String, Fasta> ac2fasta = new HashMap();
        HashMap<String, int[]> ac2seqcov = new HashMap();
        HashMap<String, IdentifiedProtein> ac2IdentifiedProtein = new HashMap();

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
                    //ac2fasta.put(myac, seq);
                    ac2fasta.put(myac, f);
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
                    for (Iterator<ProteinGroup> itr = reader.getProteinGroupList().iterator(); itr.hasNext();) {

                        ProteinGroup pg = itr.next();
//System.out.println("Number of proteins is " + pg.getNumProteins() + " and number of peptides is : " + pg.getNumPeptides());
                        for(Iterator<Protein> itp = pg.getProteins(); itp.hasNext();) {
                            Protein p = itp.next();
                            String proteinid = p.getAccession();
                            
                            if(interestedacs.contains(proteinid)) {
//System.out.println("Protein " + proteinid + " is included in the interest list" );
                                int [] seqcov = ac2seqcov.get(proteinid);
                                String protseq = ac2fasta.get(proteinid).getSequence();
                                IdentifiedProtein idp = ac2IdentifiedProtein.get(proteinid);
                                if(idp == null) {
                                    idp = new IdentifiedProtein(ac2fasta.get(proteinid));
                                    ac2IdentifiedProtein.put(proteinid, idp);
                                }

                                for(Iterator<Peptide> pepItr=pg.getPeptides(); pepItr.hasNext(); ) {
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

                                    idp.addPeptide(peptide);

                                }
                            }
                        }


                    }
                } 
                line = br.readLine();
            } 
            br.close();

            System.out.println("\n\n\n\n\nProtein Accession\tTotal Sequence Coverage\tProtein Length\tNumPeptides\tNumModifiedPeptides\tDescription");
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
                Fasta f =  ac2fasta.get(ac);
                IdentifiedProtein idp = ac2IdentifiedProtein.get(ac);
                System.out.println(ac + "\t" + (length - numzero)/(length+0.0)*100 + "%" + "\t" + f.getLength() +  "\t" + idp.getNumPeptides() + "\t" + idp.getNumModifiedPeptides() + "\t" + f.getDescription());
                //System.out.println(ac + "\t" + (length - numzero)/(length+0.0)*100 + "%" + "\t" + f.getLength() +  "\t" + f.getDescription());
            }
           
            //System.out.println("P_for_Protein_p_for_Peptide\tProtein Accession \tTotal Sequence Coverage\tProtein Length\tNumPeptides\tNumModifiedPeptides\tDescription");
            for(Iterator<String> it = interestedacs.iterator(); it.hasNext();) { 

                String ac = it.next();

                System.out.println("\n\n\nPTMs identified on Protein " + ac);

                System.out.println("Sequence\tOccurence\tXCorr\tDeltaCN\tZScore");
                int [] seqcov = ac2seqcov.get(ac);
                int numzero = 0;
                for(int i = 0; i < seqcov.length; i++) {
                    if(seqcov[i] == 0) {
                        numzero++;
                    } 
                }
                int length = seqcov.length;
                //System.out.println(ac + "\t" + length + "\t" + numzero + "\t" + (length - numzero)/(length+0.0)*100 + "%" );
                Fasta f =  ac2fasta.get(ac);
                IdentifiedProtein idp = ac2IdentifiedProtein.get(ac);
                //System.out.println("\n\n\nP\t" + ac + "\t" + (length - numzero)/(length+0.0)*100 + "%" + "\t" + f.getLength() +  "\t" + idp.getNumPeptides() + "\t" + idp.getNumModifiedPeptides() + "\t" + f.getDescription());
                for(Iterator<IdentifiedPeptide> itp = idp.getIdentifiedPeptides(); itp.hasNext();) {
                    IdentifiedPeptide idpep = itp.next();
                    if(idpep.isModifiedPeptide()) {
                        Peptide bestpeptide = idpep.getBestPeptide();
                        System.out.println(idpep.getSequence() + "\t" + idpep.getNumPeptides() + "\t" + bestpeptide.getXCorr() + "\t" + bestpeptide.getDeltCN() + "\t" + bestpeptide.getSpScore());
                    }
                }
            }

	} finally {
	    //tidy up
	    if(br != null){
		br.close();
	    }
	}
    }
}

class IdentifiedProtein {

    Fasta fasta;
    HashMap<String, IdentifiedPeptide> seq2IdentifiedPetide = new HashMap(); 
    IdentifiedProtein(Fasta f) {
        fasta = f;
    }


    int getNumPeptides() {
        return seq2IdentifiedPetide.values().size();
    }

    int getNumModifiedPeptides() {
        int num = 0;
        for(Iterator<IdentifiedPeptide> it = seq2IdentifiedPetide.values().iterator(); it.hasNext();) {
            if(it.next().isModifiedPeptide()) {
                num++;
            }
        }
        return num;
    }
    void addPeptide(Peptide p) {
        String seq = p.getSequence();
        IdentifiedPeptide ip = seq2IdentifiedPetide.get(seq);
        if(ip == null) {
            ip = new IdentifiedPeptide(seq);
            seq2IdentifiedPetide.put(seq, ip);
        } 
        ip.addPeptide(p);
    } 
    Iterator<IdentifiedPeptide> getIdentifiedPeptides() {
        return seq2IdentifiedPetide.values().iterator();
    }
}

class IdentifiedPeptide {
    ArrayList<Peptide> peptides = new ArrayList();
    String sequence;

    Peptide getBestPeptide() {
        Peptide best = null; 
        for(Iterator<Peptide> it = peptides.iterator(); it.hasNext();) {
            Peptide p = it.next();
            if(best == null || p.getXCorrValue() > best.getXCorrValue()) {
                best = p;
            }
        }
        return best;
    }

    IdentifiedPeptide(String seq) {
        sequence = seq;
    }

    void addPeptide(Peptide p) {
        peptides.add(p);
    }
    boolean isModifiedPeptide() {
        return sequence.indexOf("(") != -1;
    }
    int getNumPeptides() {
        return peptides.size();
    }
    String getSequence() {
        return sequence;
    }
}
