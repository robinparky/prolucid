package edu.scripps.pms.util.seq;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import java.io.*;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.io.DTASelectFilterReader;
import edu.scripps.pms.util.dtaselect.Peptide;
import edu.scripps.pms.util.dtaselect.Protein;
import edu.scripps.pms.util.spectrum.PeakList;
import edu.scripps.pms.util.io.SpectrumReader;


/**
 * <p>Program to six-frame translate a nucleotide sequence</p>
 */

// It now ignores frames with more than 2 stop condons
public class PrefixDb {


   
    private ArrayList<Fasta>[][][][] prefixmap = null;

    /**
     * Call this to get usage info, program terminates after call.
     */
    public static void help() {
	System.out.println("usage: java PrefixDb fasta_file");
	System.exit(-1);
    }


    public PrefixDb(String fastafile) throws IOException {
        prefixmap = sortFastaByPrefex(fastafile);
    }

    // the length of pep has to be at least 4 
    // remember to remove the special char, such as ptm symbols in the peptide sequence before calling this method
    // return the FASTA items in this database that contains the first 4 char s and the last 4 chars of the given pep
    // seq. It does NOT guarantee the returned FASTA items to contain the pep string
    public HashSet<Fasta> peptideseq2Fastas(String pep) {
        HashSet<Fasta> result = new HashSet();
        int len = pep.length();
        try{
            
            ArrayList<Fasta> prefixarr = prefixmap[pep.charAt(0)][pep.charAt(1)][pep.charAt(2)][pep.charAt(3)];
            if(prefixarr != null) {
                //result.addAll(prefixmap[pep.charAt(0)][pep.charAt(1)][pep.charAt(2)][pep.charAt(3)]);
                result.addAll(prefixarr);
            } else {

                System.err.println(pep + "\t prefix not found");
            }

            HashSet<Fasta> suffixprots = new HashSet();
            ArrayList<Fasta> sufixarr = prefixmap[pep.charAt(len-4)][pep.charAt(len-3)][pep.charAt(len-2)][pep.charAt(len-1)];
            if(sufixarr != null) {
                //suffixprots.addAll(prefixmap[pep.charAt(len-4)][pep.charAt(len-3)][pep.charAt(len-2)][pep.charAt(len-1)]);
                suffixprots.addAll(sufixarr);
            } else {
                System.err.println(pep + "\t sufix not found");
            }
            //System.out.println("Number of suffix proteins: " + suffixprots.size());
            //System.out.println("Number of prefix proteins: " + result.size());
            result.retainAll(suffixprots);
            //System.out.println("Number of proteins after suffixe: " + result.size());
        } catch (Exception e) { 
            e.printStackTrace();
            System.out.println("Illegal Char found in " + pep);
        }
        return result;
    }
    private ArrayList<Fasta>[][][][] sortFastaByPrefex(String database) throws IOException {
        ArrayList<Fasta>[][][][] prefixtree = new ArrayList[91][91][91][91];
        long starttime = System.nanoTime();
        System.out.println("Start sortFastaByPrefex");
        FileInputStream fis = new FileInputStream(new File(database.trim()));
        Iterator<Fasta> fastas = FastaReader.getFastas(fis);
        int j = 0;
        while(fastas.hasNext()) {
            Fasta f = fastas.next();

            String myac = f.getAccession();
   //System.out.println("protien ac " + myac + "\tdescription " + f.getDescription());
            //if(myac != null && !myac.startsWith("Reverse")) {
            byte [] seq = f.getSequenceAsBytes();
            for(int i = 0; i < seq.length-3; i++) {
                ArrayList al = prefixtree[seq[i]][seq[i+1]][seq[i+2]][seq[i+3]];
                if(al == null) {
                    al = new ArrayList();
                    prefixtree[seq[i]][seq[i+1]][seq[i+2]][seq[i+3]] = al;
                }
                al.add(f); 
            }
//System.out.println(++j + "\r");
        }
        long finishtime = System.nanoTime();
        fis.close(); 
        System.out.println("finished sortFastaByPrefex. time used: " + (finishtime - starttime)/1000000);
        return prefixtree;
    }

    public static void main(String[] args) throws Exception{

	BufferedReader br = null;
        String database = "/lustre/people/applications/yates/dbase/UniProt_mouse_05-17-2011_reversed.fasta";
	try {
            if(args.length > 0) {
                database = args[0];
            }
            PrefixDb pdb = new PrefixDb(database); 
            String peptideseq = "VDAASSNGPFQPV";
            HashSet<Fasta> prots = pdb.peptideseq2Fastas(peptideseq);
            if(prots != null) {
                System.out.println("Number of proteins in the database contains the peptide " + peptideseq + " is "  + prots.size());
                for(Iterator<Fasta> it = prots.iterator(); it.hasNext();) {
                    Fasta f = it.next();
                    String protseq = f.getSequence();
                    int start = protseq.indexOf(peptideseq);
                    if(start == -1) {
                        System.out.println(f.getAccession() + " does not contain " + peptideseq);
                    } else {
                        System.out.println(f.getAccession() + " containis " + peptideseq + " at " + start);

                    }
                }
            } else {
                System.out.println("No protein in the database contains the peptide " + peptideseq);
            } 

            

	} finally {
	    //tidy up
	    if(br != null){
		br.close();
	    }
	}
    }

}

