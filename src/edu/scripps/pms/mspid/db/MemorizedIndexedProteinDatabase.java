/**
 * @file MemorizedIndexedProteinDatabase.java
 * This is the source file for edu.scripps.pms.mspid.MemorizedIndexedProteinDatabase
 * @author Tao Xu
 * @date $Date: 2007/09/24 21:30:27 $
 */
package edu.scripps.pms.mspid.db;

import blazmass.dbindex.DBIndexer;
import blazmass.dbindex.IndexedProtein;
import blazmass.dbindex.IndexedSequence;
import edu.scripps.pms.mspid.*;
import edu.scripps.pms.mspid.Peptide;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.TimeUtils;
import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.ByteArrayConverter;
//import edu.scripps.pms.util.enzyme.Protease;
import java.util.*;
import java.io.*;
import java.lang.reflect.Array;

// assume all peptide index info kept in RAM 
public class MemorizedIndexedProteinDatabase extends ProteinDatabase {
    /**
     * The protein index start from 0 to number of proteins - 1
     * The .pin the same number of records as the number of peptides,
     * Each record has 4 bytes for mass (equals to (int)Math.round(realMass*ACCURACYFACTOR),
     * 4 bytes for protein index, 2 bytes for start (unsighed short, same as char), 
     * 2 bytes for end (unsighed short, same as char)
     * The .min file contains offset of masses in the .pin file, each record is a long
     * 
     */
    private TimeUtils timer = new TimeUtils();
    private static final String USAGE = "java MemorizedIndexedProteinDatabase dbfile";
    private String peptideIndexFile;
    private String massIndexFile;
    //private RandomAccessFile raf;
    private static final int ACCURACYFACTOR = 1000;
    // keep 1pp accuracy for peptides with mass 1000 dalton
    // for the frequence of peptide length
    private static final int RECORDSIZE = 8;
    //private static final int DEFAULTNUMRECORDS = 1000000; 
    private long [] index; // mass index for offsets of each mass
    private int MAXINDEX;

    private int [] proteinIndex;
    private char [] peptideStartIndex;
    private char [] peptideEndIndex;

    private DBIndexer dbIndexer;

    long totalNumPeptides = 0;
    // the following variable will be used by the heap sort for swapping
    //private byte [] b = null; // to read the records

    public MemorizedIndexedProteinDatabase(DBIndexer dbIndexer) {
        this.dbIndexer = dbIndexer;
    }

    public MemorizedIndexedProteinDatabase() {
    }


    public MemorizedIndexedProteinDatabase(String fasta) throws IOException {

        System.out.println("This version of prolucid does not support old pin min based indexing. Try on different version");
        System.exit(0);
        /*super(fasta);
        // .pin for peptide index, .min for mass index
        peptideIndexFile = fasta + ".pin";
        massIndexFile = fasta + ".min";

        timer.startTiming();
        loadMassIndex();
        System.out.println("Time used to load mass index: " + timer.getTimeUsed());

        MAXINDEX = index.length - 1;
       // b = new byte[DEFAULTNUMRECORDS*RECORDSIZE];

        timer.startTiming();
        readPeptideIndex();
        System.out.println("Time used to load peptide index: " + timer.getTimeUsed());
        */

        //raf = new RandomAccessFile(peptideIndexFile, "rws");
        //raf.seek(0);
    }
    private void readPeptideIndex() throws IOException {

        File pif = new File(peptideIndexFile); 
        int numpeptides = (int)(pif.length()/RECORDSIZE);
        proteinIndex = new int[numpeptides];    
        peptideStartIndex = new char[numpeptides];    
        peptideEndIndex = new char[numpeptides];    
        System.out.println("Number of Peptide from pin file: " + numpeptides);
        FileInputStream fis = new FileInputStream(pif);
        DataInputStream dis = new DataInputStream(new BufferedInputStream(fis));
        int i = 0;

        System.out.println("dis.available: " + dis.available());
        //while(dis.available() > 0) { // may return negative number if greater than 2G
        while(dis.available() != 0) {
            proteinIndex[i] = dis.readInt();
            peptideStartIndex[i] = dis.readChar();
            peptideEndIndex[i] = dis.readChar();

            i++;
            //System.out.println("dis.available: " + dis.available());
        }

        System.out.println("dis.available: " + dis.available());
        System.out.println("Number of Peptide read from pin file: " + i);
        dis.close();
    } 

    private void diffModSearch(ProcessedPeakList ppl, SearchResult sr, Modifications m) throws IOException {

        // temp solution, will use m.getMassShift()
        //double massShift = m.getDiffNTermMod();
        double massShift = m.getMassShift();

        SearchParams params = ppl.getSearchParams();
        int numIsotopes = params.getNumIsotopicPeaks();
        float precMassAccuracy = params.getPrecursorToleranceFloat();

//        double precMassAccuracy = params.getPrecursorTolerance();

        timer.startTiming();


        float prcMass = (float)ppl.getPrecursorMass();
        prcMass -= massShift;



       // precMassAccuracy =prcMass*precMassAccuracy/1000f;
        //prcMass -= ppl.getSearchParams().getStaticTerminalMods();

        //double acc = precMassAccuracy/1000000.0f;
       
        // all static mods except N and C terminal static mods have been 
        // considered while the database are processed, so N and C terminal
        // static mods and all diff mods need special attention here



        int i = 0;

      //  System.out.println("===============");
        do{

            //System.out.println("1111\t" + dPrecursorMass +"\t====================\t" + precMassAccuracy + "\t"+ prcMass + "\t" + i  + "\t" + MassSpecConstants.MASSDIFFC12C13);
            float dPrecursorMass = prcMass - i*MassSpecConstants.MASSDIFFC12C13;
           // System.out.println(i + "\t" + numIsotopes + "\t" + dPrecursorMass +"\t====================\t" + precMassAccuracy + "\t"+ prcMass + "\t" + ppl.getSearchParams().getHighPrecursorTolerance());


            List<IndexedSequence> sequenceList = dbIndexer.getSequences(dPrecursorMass, precMassAccuracy);

          //  if(numIsotopes == 0) {
          //      sequenceList = dbIndexer.getSequences(dPrecursorMass, precMassAccuracy+6);
         //   } else {
           //     sequenceList = dbIndexer.getSequences(dPrecursorMass, precMassAccuracy);
         //   }


            Fasta fasta = null;

            if(sequenceList.size()>0) {
                for(Iterator<IndexedSequence> itr=sequenceList.iterator(); itr.hasNext(); ) {
                    IndexedSequence seq = itr.next();

                    List<IndexedProtein> plist = this.dbIndexer.getProteins(seq);

                    fasta = new Fasta(seq.getWholeCleanSequence());

                  //  if(seq.getWholeCleanSequence().contains("WHLKTEI"))
                  //      System.out.println("@@@@@============\t" + seq.getMass() + "++\t" + seq.getWholeCleanSequence() + "\t" + seq.getSequence() + "\t" + seq.getResLeft() + "\t" + seq.getResLeft());
                    //else System.out.print(".");

                    for(Iterator<IndexedProtein> iptr=plist.iterator(); iptr.hasNext(); ) {
                        IndexedProtein eachProtein = iptr.next();
                        fasta.addDefList(eachProtein.getAccession());
                    }
                    //     System.out.println("=========" + seq.getWholeCleanSequence() + " "+ seq.getResRight() + " " + seq.getProteinIds() + "\t" + seq.getProteinDescArray() + " " + plist);
                    getPeptides(sr, seq, fasta, m);

               //     getPeptides(sr, seq, fasta);
                }
            }


            //prs[i] = new PeptidesReader(sr, startOffset, len); 
            //prs[i].start();
        } while (++i < numIsotopes);
    }

    public synchronized     SearchResult search(ProcessedPeakList ppl) throws IOException {
        SearchParams params = ppl.getSearchParams();
        int numIsotopes = params.getNumIsotopicPeaks();
        float precMassAccuracy = params.getPrecursorToleranceFloat();

        timer.startTiming();
        SearchResult sr = new SearchResult(ppl);
        float prcMass = (float)ppl.getPrecursorMass();
        precMassAccuracy =prcMass*precMassAccuracy/1000f;
        prcMass -= ppl.getSearchParams().getStaticTerminalMods();


        double minPrecMass = ppl.getSearchParams().getMinPrecursorMass();
        double maxPrecMass = ppl.getSearchParams().getMaxPrecursorMass();

        int i = 0;
        //for(int i = 0; i < numIsotopes; i++) {
        do{

            float dPrecursorMass = prcMass - i*MassSpecConstants.MASSDIFFC12C13;

            List<IndexedSequence> sequenceList = dbIndexer.getSequences(dPrecursorMass, precMassAccuracy);
            Fasta fasta = null;

            if(sequenceList.size()>0) {
                for(Iterator<IndexedSequence> itr=sequenceList.iterator(); itr.hasNext(); ) {
                    IndexedSequence seq = itr.next();

                    List<IndexedProtein> plist = this.dbIndexer.getProteins(seq);

                    fasta = new Fasta(seq.getWholeCleanSequence());
           //         System.out.println("============\t" + seq.getMass() + "++\t" + seq.getWholeCleanSequence() + "\t" + seq.getSequence() + "\t" + seq.getResLeft() + "\t" + seq.getResLeft());

                    for(Iterator<IndexedProtein> iptr=plist.iterator(); iptr.hasNext(); ) {
                        IndexedProtein eachProtein = iptr.next();
                        fasta.addDefList(eachProtein.getAccession());
                    }
               //     System.out.println("=========" + seq.getWholeCleanSequence() + " "+ seq.getResRight() + " " + seq.getProteinIds() + "\t" + seq.getProteinDescArray() + " " + plist);
                    getPeptides(sr, seq, fasta);
                }
            }
        } while (++i < numIsotopes);

        for(Iterator<Modifications> it = ppl.getModifications(); it.hasNext();) {
            Modifications m = it.next();
            if(m != null && m.getDiffModsShift() != 0 ) {
//System.out.println("modification searches, massShfit: " + m.getDiffModsShift());
                diffModSearch(ppl, sr, m);
            }

        }
//System.out.println("diffmodshift: " + m.getDiffModsShift());
        return sr;
    }


    private void getPeptides(SearchResult sr, IndexedSequence seq, Fasta fasta, Modifications m) throws IOException {

        int size = seq.getSequence().length();
        int start = 3;
        int end = start+size-1;

        if(sr.getProcessedPeakList().isDeCharged()) {
            //start += 3;
            sr.addPeptideHit(new DeChargedModifiedPeptideHit(fasta, start, end, m));
        } else {
            sr.addPeptideHit(new ModifiedPeptideHit(fasta, start, end, m));
        }
    }

    private void getPeptides(SearchResult sr, IndexedSequence seq, Fasta fasta) throws IOException {
       //    System.out.println("processing " + proteinIndex[i] + " for " + sr);
            //sr.addPeptideHit(new PeptideHit(getFasta(proteinIndex[i]), peptideStartIndex[i], peptideEndIndex[i]));

        int size = seq.getSequence().length();
        int start = 3;
        int end = start+size-1;

        if(sr.getProcessedPeakList().isDeCharged()) {
            sr.addPeptideHit(new DeChargedPeptideHit(fasta, start, end));

        } else {
            sr.addPeptideHit(new PeptideHit(fasta, start, end));
        }
    }



    public SearchResult search_orig(ProcessedPeakList ppl) throws IOException {
        SearchParams params = ppl.getSearchParams();
        int numIsotopes = params.getNumIsotopicPeaks();
        double precMassAccuracy = params.getPrecursorTolerance();

        timer.startTiming();
        SearchResult sr = new SearchResult(ppl);
        double acc = precMassAccuracy/1000000.0f;

        // all static mods except N and C terminal static mods have been
        // considered while the database are processed, so N and C terminal
        // static mods and all diff mods need special attention here

        double prcMass = ppl.getPrecursorMass();

//System.out.println("staticTerminalMods: " + ppl.getSearchParams().getStaticTerminalMods());
        prcMass -= ppl.getSearchParams().getStaticTerminalMods();


        //MassCalculator mc = ppl.getMassCalculator();
        // for multi threading
        //PeptidesReader [] prs = new PeptidesReader[numIsotopes];
        int i = 0;
        //for(int i = 0; i < numIsotopes; i++) {
        do{
            double mass = 1000*(prcMass - i*MassSpecConstants.MASSDIFFC12C13);
            double diff = mass * acc;
            diff = diff < 500? diff : 500;
            int highLimit = (int)Math.round(mass + diff);
            int lowLimit = (int)Math.round(mass - diff);

            if(numIsotopes == 0) { // traditional sequest like
                highLimit = (int)Math.round(mass + ppl.getSearchParams().getHighPrecursorTolerance());
                lowLimit = (int)Math.round(mass - ppl.getSearchParams().getLowPrecursorTolerance());
                //System.out.println("highLimit: " + highLimit + "\tlowLimit: " + lowLimit);
            }

            System.out.println("==========" + highLimit + " " + lowLimit);

            if(highLimit > MAXINDEX) { highLimit = MAXINDEX; }
            if(lowLimit > MAXINDEX) { lowLimit = MAXINDEX; }

            System.out.println("==========" + highLimit + " " + lowLimit);

//            System.out.print("highLimit: " + highLimit + "\tlowLimit: " + lowLimit);
            long startOffset = index[lowLimit-1];
            long endOffset = index[highLimit];
            int len = (int)(endOffset - startOffset);
//            System.out.println("\tstartOffSet: " + startOffset + "\tendOffset: " + endOffset + "\tlen: " + len);

            if(len > 0) {
                //getPeptides(sr, readPeptides(startOffset, len));
                getPeptides(sr, startOffset/RECORDSIZE, endOffset/RECORDSIZE);
            }
            //prs[i] = new PeptidesReader(sr, startOffset, len);
            //prs[i].start();
        } while (++i < numIsotopes);
        for(Iterator<Modifications> it = ppl.getModifications(); it.hasNext();) {
            Modifications m = it.next();
            if(m != null && m.getDiffModsShift() != 0 ) {
//System.out.println("modification searches, massShfit: " + m.getDiffModsShift());
                diffModSearch(ppl, sr, m);
            }

        }
//System.out.println("diffmodshift: " + m.getDiffModsShift());
        return sr;
    }

    // for diff mods search
    private void getPeptides(SearchResult sr, long start, long end, Modifications m) throws IOException {
     
        for(int i = (int)start; i <= (int)end; i++) {
            //sr.addPeptideHit(new ModifiedPeptideHit(getFasta(proteinIndex[i]), peptideStartIndex[i], peptideEndIndex[i], m));
            addPeptideHit(sr, getFasta(proteinIndex[i]), peptideStartIndex[i], peptideEndIndex[i], m);
        }       
    }
    private void getPeptides(SearchResult sr, long start, long end) throws IOException {
        for(int i = (int)start; i <= (int)end; i++) {

            //    System.out.println("processing " + proteinIndex[i] + " for " + sr);
            //sr.addPeptideHit(new PeptideHit(getFasta(proteinIndex[i]), peptideStartIndex[i], peptideEndIndex[i]));
            addPeptideHit(sr, getFasta(proteinIndex[i]), peptideStartIndex[i], peptideEndIndex[i]);
        }  
    }
    private void loadMassIndex() throws IOException {
        DataInputStream dis = new DataInputStream(new BufferedInputStream
                                     (new FileInputStream(massIndexFile)));
        int numIndex = dis.available()/8; // 64 is number of byte for a long
        //System.out.println("numIndex: " + numIndex);
        //timer.startTiming();
        index = new long[numIndex];
        for(int i = 0; i < numIndex; i++) {
            index[i] = dis.readLong();
            //if(i > 1587000)
           // System.out.println(i +"\t" + index[i]);
        }
        dis.close();
        //System.out.println("Time used to load mass index: " + timer.getTimeUsed());        
    }
    public static void main(String args[]) throws Exception {
        try { 
            TimeUtils timer = new TimeUtils();
            String fasta = args[0];
            timer.startTiming();
            MemorizedIndexedProteinDatabase se = new MemorizedIndexedProteinDatabase(fasta);
            // peptide index file
            //se.outputPeptides(); 
            //se.sortPeptides(); 
            for(int i = 100000; i < se.index.length; i += 100000) { 
                long size = se.index[i] - se.index[i-1000];
                System.out.println("mass: " + i/1000 + "\tsize: " + size);
            }
            timer.stopTiming();
            long timeUsed = timer.getTimeUsedMillis(); 
            System.out.println("Time used for process database: " + timeUsed);
            System.out.println();
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
        }
    }

    public DBIndexer getDbIndexer() {
        return dbIndexer;
    }

    public void setDbIndexer(DBIndexer dbIndexer) {
        this.dbIndexer = dbIndexer;
    }
}

