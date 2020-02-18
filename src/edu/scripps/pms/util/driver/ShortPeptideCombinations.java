/**
 * @file Accession2Fastas.java
 * This is the source file for Accession2Fastas class
 * @author Tao Xu
 * @date $Date
 */

import edu.scripps.pms.mspid.ProteinDatabase;

import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.io.FastaReader;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.HashMap;
import java.util.Set;
import java.io.FileInputStream;
import java.io.File;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.IOException;
import java.io.*;
import java.util.*;


public class ShortPeptideCombinations {

    public static final String USAGE = "!!! Usage: java ShortPeptideCombinations databaseFile coverage (e.g., 0.95) !!!";
    // databaseName - the path and name of the protein database

    //public static final char [] residues = {'G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'I', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W'};
    public static final char [] residues = {'A', 'C', 'E', 'D', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
    // assuming the mass tolerance is smaller than the mass of any AA residue
    private static ArrayList<String> shortPeptides = null;
    private static String proteindb = "";
    private static int numResidues = 3;
    private static double coverage = 0.90;


    private static ArrayList<String> getShortPeptides() {
        if(shortPeptides != null) {
            return shortPeptides;
        }
        ArrayList<String> peptides = new ArrayList();
        ArrayList<String> tempold = new ArrayList();

        tempold.add(""); 
        for(int i = 0; i < numResidues; i++) {
            ArrayList<String> temp = new ArrayList();
            
            for(Iterator<String> it = tempold.iterator(); it.hasNext();) {
                String s = it.next();
                if(s.length() == i) {
                    for(int j = 0; j < residues.length; j++) {
                        temp.add(s+residues[j]);
                    }
                    
                }
            }
            tempold = temp;
        }
        for(Iterator<String> it = tempold.iterator(); it.hasNext();) {
            String s = it.next(); 
            if(s != null && s.length() == numResidues) {
                 peptides.add(s);

                 //System.out.println("Peptide: " + s);
            }
        }
        System.out.println("Number of Peptides added: " + peptides.size());
        return peptides; 
    }

    public static void main(String [] args) throws Exception {
        try {
            String databaseName = args[0];
            if(args.length > 1) {
                coverage = Double.parseDouble(args[1]);
            }
            //numResidues = Integer.parseInt(args[1]);
            ArrayList<String> peptides = getShortPeptides();
         
            HashMap<String, PeptideAccs> peptide2accs = new HashMap(1000000); 
            HashMap<String, AccPeptideAccs> acc2peptides = new HashMap(1000000); 
            FileInputStream fis = new FileInputStream(new File(databaseName));
            Iterator<Fasta> fastas = FastaReader.getFastas(fis);
            int i = 0;
         
            for(Iterator<String> it = peptides.iterator(); it.hasNext();) {
                String pepseq = it.next();
                peptide2accs.put(pepseq, new PeptideAccs(pepseq));
            } 
            HashSet<String> proteinseqs = new HashSet<String>(10000000);
            int numfasta = 0;
            int numuniqueprotein = 0;
            StringBuffer duplicates = new StringBuffer(1000000);
            while(fastas.hasNext()) {

                numfasta++; 
                Fasta f = fastas.next();
                String myac = f.getAccession();
                //String myacnoversion = f.getAccessionWithNoVersion();
                //System.out.print("Checking " + myac + "\r"); 
                String seq = f.getSequence();
                if(proteinseqs.add(seq)) { // cannot add means dupicate sequence
                    AccPeptideAccs apa = new AccPeptideAccs(myac);

                    numuniqueprotein++;
                     
                    for(Iterator<String> it = peptides.iterator(); it.hasNext();) { 
                        String peptide = it.next();
                        //if(seq.indexOf(pepitde) > -1) {
                        if(seq.contains(peptide)) {
                            PeptideAccs pa = peptide2accs.get(peptide);
                            pa.addAcc(myac);
                            //peptide2accs.get(peptide).addAcc(myac);
                            apa.addPeptide(pa);
                         
                        }        

                    } 
                    acc2peptides.put(myac, apa);
                //System.out.println(">" + f.getDefline());
                } else {
                    //System.out.println("Duplicate protein sequence\t" + myac);
                    duplicates.append("Duplicate protein sequence\t" + myac + "\n");
                }

            }
            fis.close();


//            for(Iterator<String> it = peptides.iterator(); it.hasNext();) {
//                String peptide = it.next();
//                System.out.println(peptide + "\t" + peptide2accs.get(peptide).getAccs().size());
//            }

        
            ArrayList<PeptideAccs> paset = new ArrayList<PeptideAccs>();
            paset.addAll(peptide2accs.values());   
            
            Collections.sort(paset); 
            HashSet<String> addedaccs = new HashSet<String>(100000);
            int numaccsadded = 0;
            int additionalpeptides = 25;
           
            ArrayList<PeptideAccs> peptideused = new ArrayList<PeptideAccs>();

            System.out.println("Sequence\tSize before adding\tsize after adding\tnumber of new accs added\tCoverage");
            for(Iterator<PeptideAccs> it = paset.iterator(); it.hasNext();) {
                PeptideAccs pa = it.next();
                //pa = it.next();
                //pa = it.next();

                peptideused.add(pa);
 
                numaccsadded++;
                int sizebefore = addedaccs.size();
                addedaccs.addAll(pa.getAccs());
                int sizeafter = addedaccs.size();
                int newlyadded = sizeafter - sizebefore;
                System.out.println(pa.getSeq() + "\t" + sizebefore + "\t" + sizeafter + "\t" + newlyadded + "\t" + (sizeafter + 0.0)/numfasta );
                if(sizeafter >= numuniqueprotein*coverage) additionalpeptides--;
                if(additionalpeptides == 0) break;
 
                //System.out.println("After sorting: " + pa.getSeq() + "\t" + pa.getAccs().size());
            }
            
            HashSet<String> uniquestatus = new HashSet<String>(100000); 
            int numnotunique = 0;
            for(Iterator<String> it = acc2peptides.keySet().iterator(); it.hasNext();) {
                String acc = it.next();
                AccStatus as = new AccStatus(acc);
                AccPeptideAccs apa = acc2peptides.get(acc);
                for(Iterator<PeptideAccs> pit = peptideused.iterator(); pit.hasNext();) {
                    //PeptideAccs pa = peptide2accs.get(it.next());
                    PeptideAccs pa = pit.next();
                    if(apa.contains(pa)) {
                        as.setStatus(true);                        
                    } else {
                        
                        as.setStatus(false);                        
                    }
                }
                if(!uniquestatus.add(as.getStatus())) {
                    numnotunique++;
                    System.out.println(as.getAcc() + " does not have uqnitue status\t" + as.getStatus());
                }
            }

            System.out.println("Number of protein with unique status: " + uniquestatus.size() + " number of repeated status: " + numnotunique);
            System.out.println("Used " + numaccsadded + " Peptides to cover " + coverage + " of proteome in " + args[0]);
            System.out.println("\n\nNumFasta\t" + numfasta + "\tNumUniqueProteinSeq\t" + numuniqueprotein);
            System.out.println(duplicates);


        } catch (Exception e) {
            System.err.println(USAGE);
            e.printStackTrace(); 
        }
    }
}

class AccStatus {
    private String acc;
    private StringBuffer status = new StringBuffer(50);
    private String statusString = null;    
 
    public AccStatus(String s) {
        acc = s;
    }
    public String getAcc() {
        return acc;
    }
    public void setStatus(boolean b) {
        if(b) {
            status.append('1');
        } else {

            status.append('0');
        }
    }
    public String getStatus() {
       if(statusString == null) {
           statusString = status.toString();
       }
       return statusString;
    }
}

class AccPeptideAccs implements Comparable {
    private String acc;
    private HashSet<PeptideAccs> pas = new HashSet(500);

    public AccPeptideAccs(String a) {
        acc = a;
    }

    public void addPeptide(PeptideAccs pa) {
        pas.add(pa);
    }

    public void addAccs(Set<PeptideAccs> accset) {
        for(Iterator<PeptideAccs> it = accset.iterator(); it.hasNext();) {
            pas.add(it.next());
        }
    }
    public String getAcc() {
        return acc;
    }
    public int getNumPeptides() {
        return pas.size();
    }
    public HashSet<PeptideAccs> getPeptides() {
        return pas;
    }
    public boolean contains(PeptideAccs pa) {
        return pas.contains(pa);
    }
    
    public int compareTo(Object o) {
        AccPeptideAccs p2a = (AccPeptideAccs) o;
        int diff = this.getNumPeptides() - p2a.getNumPeptides();
        if(diff == 0) {
            return 0;
        } else if(diff > 0) {
            return -1;
        } else {
            return 1;
        }
    }
    
}

class PeptideAccs implements Comparable {

    private String seq;
    private HashSet<String> accs = new HashSet(500);

    public PeptideAccs(String s) {
        seq = s;
    }

    public void addAcc(String acc) {
        accs.add(acc);
    }

    public void addAccs(Set<String> accset) {
        for(Iterator<String> it = accset.iterator(); it.hasNext();) {
            accs.add(it.next());
        }
    }
    public String getSeq() {
        return seq;
    }
    public int getNumAccs() {
        return accs.size();
    }
    public HashSet<String> getAccs() {
        return accs;
    }

    public int compareTo(Object o) {
        PeptideAccs p2a = (PeptideAccs) o;
        int diff = this.getNumAccs() - p2a.getNumAccs();
        if(diff == 0) {
            return 0;
        } else if(diff > 0) {
            return -1;
        } else {
            return 1;
        }
    }
      

 
}
