//#package edu.scripps.pms.util.driver;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import java.io.*;

import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.mspid.MassCalculator;
import edu.scripps.pms.mspid.MassSpecConstants;
import edu.scripps.pms.util.enzyme.Protease;

/**
 * <p>Program to six-frame translate a nucleotide sequence</p>
 */

// It now ignores frames with more than 2 stop condons
public class TrypticPeptideLibrary {

static int totalpep = 0;
static int uniquepep = 0;
static HashSet<String> pepseqs = new HashSet(1000000);
    /**
     * Call this to get usage info, program terminates after call.
     */
    public static void help() {
	System.out.println(
		"usage: proteintrypticpeptides <fastafile>");
	System.exit( -1);
    }

    //public static int processFasta(Fasta f, Protease p, MassCalculator mc, double highLimit, double lowLimit, int minlength) {
    public static int processFasta(Fasta f, MassCalculator mc, double highLimit, double lowLimit, int minlength, HashMap<String, ArrayList<String>> pepseq2acs) {
        double [] precMasses = mc.getPrecMasses();
        int counts = 0;
        int proteinlength =  f.getLength();
        int lastIndex = f.getLength() - 1;
        byte [] seq = f.getSequenceAsBytes();
            // mass of H2O need to be added
            // because the mass of the residues does not count H2O
        double mass = MassSpecConstants.MASSH3O;
        double tempmass = 0;
        int templeft = 0;
        //int rightEnd = -1;
        int leftEnd = 0;

        while(leftEnd < proteinlength) {

            int rightEnd = leftEnd;
            
            while(rightEnd < lastIndex) {
                rightEnd++;
                if(rightEnd == lastIndex || seq[rightEnd-1] == 'K' ||  seq[rightEnd-1] == 'R') { // tryptic end
                    if(rightEnd - leftEnd > minlength-1) {
                        double precmass = mc.getPrecursorMass(seq, leftEnd, rightEnd);

                        if(precmass >= lowLimit && precmass <= highLimit) {
                            counts++;
                            String pepseq = f.getSequence().substring(leftEnd, rightEnd); 
                            ArrayList<String> acs = pepseq2acs.get(pepseq);
                            totalpep++;
pepseqs.add(pepseq);
                            if(acs == null) {
                                acs = new ArrayList();
 //                               pepseq2acs.put(pepseq, acs);
                                uniquepep++;
                            }
   //                         acs.add(f.getAccession());
                        }
                    }
                    break;
                }
            }
            leftEnd = rightEnd+1;

        }

        return counts;

    }

    public static void main(String[] args) throws Exception{

        
        HashSet<String> acset = new HashSet();
        ArrayList<String> aclist = new ArrayList();
        HashMap<String, ArrayList<String>> pepseq2acs = new HashMap();

        //Protease p = new Protease("Trypsin");
        //p.setType(true);
        //p.addCleavageSite('R');
        //p.addCleavageSite('K');

        MassCalculator mc = new MassCalculator();


	try {

            String fastaFileName = args[0];
int numfastaprocessed = 0;
           
            HashMap<String, ArrayList<String>> pepseq2ac2 = new HashMap(); 
            FileInputStream fis = new FileInputStream(new File(fastaFileName));
            for (Iterator<Fasta> it = FastaReader.getFastas(fastaFileName); it.hasNext();) {
                Fasta f = it.next();
                String acc = f.getAccession();
//                aclist.add(acc);
//                acset.add(acc);
                processFasta(f, mc, 4800, 600, 6, pepseq2ac2);
 System.out.println(++numfastaprocessed + "\t" + totalpep + "\t" + uniquepep + "\t" + pepseqs.size()); 
            }
            fis.close();
System.out.println("Number of ac2: " + aclist.size());
System.out.println("Number of unique ac2: " + acset.size());
System.out.println("Number of peptides: " + pepseq2ac2.size());

	} catch(Exception e){
            help();
            e.printStackTrace(); 
	} finally {
	    //tidy up
	}
    }
}

