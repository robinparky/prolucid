/**
 * @file ProteinDatabase.java
 * This is the source file for edu.scripps.pms.util.spectrum.ProteinDatabase
 * @author Tao Xu
 * @date $Date: 2013/01/08 00:36:00 $
 */



package edu.scripps.pms.mspid;

import blazmass.dbindex.DBIndexer;
import blazmass.dbindex.DBIndexerNoSQL;
import blazmass.dbindex.IndexedProtein;
import blazmass.dbindex.IndexedSequence;
import edu.scripps.pms.protinf.ProteinData;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.enzyme.Protease;
import gnu.trove.*;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.set.TIntSet;

import java.util.*;
import java.io.FileInputStream;
import java.io.File;
import java.io.IOException;
import java.io.*;

public class ProteinDatabase implements Searchable {

    private int maxMisCleavage = -1; // -1 means unlimited miscleavage allowed
    private Protease protease;
    private int enzymeSpecificity = 0;
    protected ArrayList<Fasta> sequences = new ArrayList<Fasta>(100000);
    //protected HashMap<String, Fasta> ac2Fasta = new HashMap<String, Fasta>(100000);
    protected HashMap<String, Fasta> ac2Fasta = null;
    
    private int[][] freq = new int[5][41];   
    //private boolean isDeCharged = false; 
    private double [] precMasses;
    //private DBIndexer dbIndexer;
    private DBIndexer dbIndexer;
    private HashMap<Integer, HashMap<String, IndexedSequence>> dbMap;
    private int startRange=600;
    private int endRange=6000;
    HashMap<Integer,String> proteinMap;

    public ProteinDatabase() {}

    public  ProteinDatabase(DBIndexer dbIndexer) {
        this.dbIndexer = dbIndexer;
    }

    public ProteinDatabase(String databaseName) throws IOException {
        FileInputStream fis = new FileInputStream(new File(databaseName));
        Iterator<Fasta> fastas = FastaReader.getFastas(fis); 
        int i = 0;

        while(fastas.hasNext()) {
            //System.out.print("\r" + ++i);
            Fasta f = fastas.next();
            sequences.add(i++, f);
            //System.out.println("ac: " + f.getAccession() + "\tdefline: " + f.getDefline());
        }
      //  //System.out.println("Number of proteins in the database: " + i);
        fis.close(); 
        //fastas = null; // let gc remove the Iterator object
    }

    // databaseName - the path and name of the protein database
    public ProteinDatabase(String databaseName, HashMap<Integer, HashMap<String, IndexedSequence>> dbMap, HashMap<Integer,String> proteinMap, int startRange, int endRange) throws IOException {
        FileInputStream fis = new FileInputStream(new File(databaseName));
        Iterator<Fasta> fastas = FastaReader.getFastas(fis); 
        int i = 0;
        this.dbMap = dbMap;
        this.proteinMap = proteinMap;
        this.startRange = startRange;
        this.endRange = endRange;
        while(fastas.hasNext()) {
            //System.out.print("\r" + ++i);
            Fasta f = fastas.next();
            sequences.add(i++, f);
            //System.out.println("ac: " + f.getAccession() + "\tdefline: " + f.getDefline());
        }
      //  //System.out.println("Number of proteins in the database: " + i);
        fis.close(); 
        //fastas = null; // let gc remove the Iterator object
    }


    public void setMaxMisCleavage(int mc) {
        maxMisCleavage = mc;
    }
    public void setProtease(Protease p) {
        protease = p;
    }
    public void setEnzymeSpecificity(int es) {
        enzymeSpecificity = es;
    }
    protected void addPeptideHit(SearchResult sr, Fasta f, int start, int end) {
        if(enzymeSpecificity == 0 || (protease.checkEnzymeSpecificityStrict(f, start, end) >= enzymeSpecificity)) {
            if(maxMisCleavage == -1 || (protease.getNumInternalMissCleavage(f, start, end) <= maxMisCleavage)) {
                sr.addPeptideHit(f, start, end);
            }
        }
    

    }
    protected void addPeptideHit(SearchResult sr, Fasta f, int start, int end, Modifications m) {

        if(enzymeSpecificity == 0 || (protease.checkEnzymeSpecificityStrict(f, start, end) >= enzymeSpecificity)) {
            if(maxMisCleavage == -1 || (protease.getNumInternalMissCleavage(f, start, end) <= maxMisCleavage)) {
                sr.addPeptideHit(f, start, end, m);
            }
        }
    

    }
    private void populateAc2Fasta() {
        ac2Fasta = new HashMap<String, Fasta>(1000000);
        for(Iterator<Fasta> it = sequences.iterator(); it.hasNext();) {
            Fasta f = it.next();
            ac2Fasta.put(f.getAccession(), f);
        }
    }
    public Fasta accession2Fasta(String ac) {
        if(ac2Fasta == null) {
            populateAc2Fasta();
        }
        return ac2Fasta.get(ac);        
    }
    public int getNumSequences() {
        return sequences.size();
    }
    public Iterator<Fasta> getFastas() {
        return sequences.iterator();
    }
    public Fasta getFasta(int index) {
        return sequences.get(index);
    }
   
    /* 
    protected PeptideHit createPeptideHit(Peptide p) {
        if(isDeCharged) {
            return new DeChargedPeptideHit(p);
        } else {
           return new PeptideHit(p);
        }
    }
    protected ModifiedPeptideHit createModifiedPeptideHit(Peptide p, Modifications m) {
        if(isDeCharged) {
            return new DeChargedModifiedPeptideHit(p, m);
        } else {
           return new ModifiedPeptideHit(p, m);
        }
    }
    */
    protected void getPeptideHits(SearchResult result,  
                             double highLimit, double lowLimit) {
        for(Fasta f : sequences) {
            int lastIndex = f.getLength() - 1;
            byte [] seq = f.getSequenceAsBytes();
             // mass of H2O need to be added 
             // because the mass of the residues does not count H2O   
            double mass = MassSpecConstants.MASSH3O;
            double tempmass = 0;
            int templeft = 0; 
            int rightEnd = -1;
            int leftEnd = 0; 
            while(leftEnd < lastIndex && rightEnd < lastIndex) { 
                while (mass < lowLimit && rightEnd < lastIndex) {
                    mass += precMasses[seq[++rightEnd]]; 
                }
                if (mass < lowLimit) {
                    break;
                }
                tempmass = mass;
                templeft = leftEnd;
                while(mass == tempmass) { // ignore non-AA characters
                    mass -= precMasses[seq[leftEnd++]]; 
                }
               
                // now the mass should be >= lowLimit
                if (tempmass <= highLimit) {  
                    //result.addPeptideHit(new PeptideHit(f, leftEnd, rightEnd));
                    //result.addPeptideHit(createPeptideHit(new Peptide(f, templeft, rightEnd)));
                    //result.addPeptideHit(f, templeft, rightEnd);
                    addPeptideHit(result, f, templeft, rightEnd);
                } 

                // check the last residue in case removed the left the mass is within tolerance
                // because this is the last peptide, so no need to adjust the leftEnd and mass
                if(rightEnd == lastIndex && mass <= highLimit && mass >= lowLimit) {
    
                    addPeptideHit(result, f, leftEnd, rightEnd);
                }

            }
        }
    }
    protected void getPeptideHits(SearchResult result,  
                             double highLimit, double lowLimit, Modifications m) {
        for(Fasta f : sequences) {
//System.out.println("in diff");
            int lastIndex = f.getLength() - 1;
            byte [] seq = f.getSequenceAsBytes();
             // mass of H2O need to be added 
             // because the mass of the residues does not count H2O   
            
            double mass = MassSpecConstants.MASSH3O;
            double tempmass = 0;
            int templeft = 0;
            int rightEnd = -1;
            int leftEnd = 0; 
            while(leftEnd < lastIndex && rightEnd < lastIndex) { 
                while (mass < lowLimit && rightEnd < lastIndex) {
                    mass += precMasses[seq[++rightEnd]]; 
                }
                if (mass < lowLimit) {
                    break;
                }
                tempmass = mass;
                templeft = leftEnd;
                
                while(mass == tempmass) {
                    mass -= precMasses[seq[leftEnd++]]; 
                }
                // now the mass should be >= lowLimit
                //if (tempmass <= highLimit && rightEnd - templeft > params.getMinimumPeptideLength()-2) {
                if (tempmass <= highLimit) {  
                    //result.addPeptideHit(new ModifiedPeptideHit(f, templeft, rightEnd, m));
                    //result.addPeptideHit(createModifiedPeptideHit(new Peptide(f, templeft, rightEnd), m));
                    //result.addPeptideHit(f, templeft, rightEnd, m);
                    addPeptideHit(result, f, templeft, rightEnd, m);
                }

                // check the last residue in case removed the left the mass is within tolerance
                // because this is the last peptide, so no need to adjust the leftEnd and mass
                if(rightEnd == lastIndex && mass <= highLimit && mass >= lowLimit) {
    
                    addPeptideHit(result, f, leftEnd, rightEnd, m);
                }
            }
        }
    }
    protected void getPeptideHits(SearchResult result, double [] highLimits, 
                        double [] lowLimits, int numPeaks, Modifications m) {
        int num = 0; 
        double lowLimit = lowLimits[numPeaks-1];
        double highLimit = highLimits[0];
//System.out.println("low: " + lowLimit + "\thigh: " + highLimit);
        for(Fasta f : sequences) {
           //num++;
        //System.out.println("numproteins: " + num);
            int lastIndex = f.getLength() - 1;
            byte [] seq = f.getSequenceAsBytes();
         // mass of H2O need to be added 
         // because the mass of the residues does not count H2O   
            double mass = MassSpecConstants.MASSH3O;
            double tempmass = 0;
            int templeft = 0; 
            int rightEnd = -1;
            int leftEnd = 0; 
            while(leftEnd < lastIndex && rightEnd < lastIndex) { 
                while (mass < lowLimit && rightEnd < lastIndex) {
                    mass += precMasses[seq[++rightEnd]]; 
                }
                if (mass < lowLimit) {
                    //return;
                    break;
                }
                tempmass = mass;
                templeft = leftEnd;
                while(mass == tempmass) {
                    mass -= precMasses[seq[leftEnd++]]; 
                }

                // now the mass should be >= lowLimit
                if (tempmass <= highLimit) {  
                    for(int i = 0; i < numPeaks; i++) {
                        if(tempmass >= lowLimits[i] && tempmass <= highLimits[i]) {
                            //result.addPeptideHit(new ModifiedPeptideHit(f, templeft, rightEnd, m));
                            //result.addPeptideHit(createModifiedPeptideHit(new Peptide(f, templeft, rightEnd), m));
                            //result.addPeptideHit(f, templeft, rightEnd, m);
                            addPeptideHit(result, f, templeft, rightEnd, m);
                            break;
                        }
                    }
                } 

                // check the last residue in case removed the left the mass is within tolerance
                // because this is the last peptide, so no need to adjust the leftEnd and mass
                if(rightEnd == lastIndex && mass <= highLimit && mass >= lowLimit) {
                    for(int i = 0; i < numPeaks; i++) {
                        //if(tempmass >= lowLimits[i] && tempmass <= highLimits[i]) {
                        if(mass >= lowLimits[i] && mass <= highLimits[i]) {
                            addPeptideHit(result, f, leftEnd, rightEnd, m);
                            break;
                        }
                    }
                }
            }
        }
    }

    private synchronized void diffModSearch(ProcessedPeakList ppl, SearchResult sr, Modifications m) throws IOException {


        
         double massShift = m.getMassShift();

        SearchParams params = ppl.getSearchParams();
        int numIsotopes = params.getNumIsotopicPeaks();
        double precMassAccuracy = params.getPrecursorTolerance();
        double acc = precMassAccuracy/1000000.0f;
        double prcMass = ppl.getPrecursorMass() - massShift;
        
        
                List<Integer> highLimits = new ArrayList<>();
        List<Integer> lowLimits = new ArrayList<>();
        
        if(numIsotopes == 0) { // for low resolution, traditional sequest like  
            int highLimit =(int) ((prcMass + params.getHighPrecursorTolerance()/1000)*1000);
            int lowLimit =(int) ((prcMass - params.getLowPrecursorTolerance()/1000)*1000);
            lowLimits.add(lowLimit);
            highLimits.add(highLimit);
        } else if(numIsotopes == 1) { // for deisotoped high resolution data 
            double diffs = prcMass*acc;
            int highLimit = (int)((prcMass + diffs)*1000);
            int lowLimit = (int)((prcMass - diffs)*1000);
            lowLimits.add(lowLimit);
            highLimits.add(highLimit);
        } else { //for non deisotoped high resolution data
            
            double diffs = prcMass*acc*2;
            for(int i = 0; i < numIsotopes; i++) {
                //lowLimits[i] = prcMass - diffs/2 - i*MassSpecConstants.MASSDIFFC12C13; 
                //highLimits[i] = lowLimits[i] + diffs; 
                
                lowLimits.add((int)((prcMass - diffs/2 - i*MassSpecConstants.MASSDIFFC12C13)*1000));
                highLimits.add((int)((lowLimits.get(i)+diffs*1000)));
            }
            
        }
        
        
        
        
        int i = 0;
           

        //double minPrecMass = ppl.getSearchParams().getMinPrecursorMass();
        //double maxPrecMass = ppl.getSearchParams().getMaxPrecursorMass();

        do{
            float dPrecursorMass = (float) (prcMass - i*MassSpecConstants.MASSDIFFC12C13);
           // System.out.println(""+dPrecursorMass);
            if(dPrecursorMass<((startRange/1000)) || dPrecursorMass>((endRange/1000))) continue;
            //System.out.println(""+dbMap.keySet().size());
             List<Integer> keyset = new ArrayList<>(dbMap.keySet());
             Integer [] keys = keyset.toArray(new Integer[keyset.size()]);
             Arrays.sort(keys);
             //List<IndexedSequence> sequenceList = new ArrayList<>();
             Fasta fasta = null;
             //int startpoint = closest((int)((dPrecursorMass+0.0005)*1000-(float)precMassAccuracy), keys);
             //int endpoint = (int)((dPrecursorMass+0.0005)*1000+(float)precMassAccuracy);
             int startpoint = closest(lowLimits.get(i),keys);
             int endpoint = highLimits.get(i);
             //System.out.println("diff search\t"+dPrecursorMass);
             //System.out.println(keys[startpoint] +" \t"+endpoint);
             for(int j =startpoint;j<keys.length;j++){
                 if(keys[j] > keys[keys.length-1] || keys[j] > endpoint){
                     break;
                 }
               // System.out.println(""+j);
                 HashMap<String,IndexedSequence>seqMap = dbMap.get(keys[j]);
                 
                 
                 Iterator g = seqMap.values().iterator();
                 
                 while(g.hasNext()){
                     IndexedSequence iseq = (IndexedSequence) g.next();
                     List<Integer> index = iseq.getProteinIds();
                     
                     //System.out.println(""+iseq.getSequence());
                     
                     
                    //List<String> proList = new ArrayList<>();
                     fasta = new Fasta(iseq.getWholeCleanSequence());
                     for(int h : index){
                         fasta.addDefList(blazmass.io.Fasta.getSequestLikeAccession(proteinMap.get(h)));
                     }
                     getPeptides(sr, iseq, fasta, m);
                 }
             }
            
/*
            //System.out.println("1111\t" + dPrecursorMass +"\t====================\t" + precMassAccuracy + "\t"+ prcMass + "\t" + i  + "\t" + MassSpecConstants.MASSDIFFC12C13);
            float dPrecursorMass = prcMass - i*MassSpecConstants.MASSDIFFC12C13;

            if(dPrecursorMass<minPrecMass || dPrecursorMass>maxPrecMass) continue;

            // System.out.println(i + "\t" + numIsotopes + "\t" + dPrecursorMass +"\t====================\t" + precMassAccuracy + "\t"+ prcMass + "\t" + ppl.getSearchParams().getHighPrecursorTolerance());

 List<Integer> keyset = new ArrayList<>(dbMap.keySet());
             Integer [] keys = keyset.toArray(new Integer[keyset.size()]);
             Arrays.sort(keys);
             List<IndexedSequence> sequenceList = new ArrayList<>();
             Fasta fasta = null;
             int startpoint = closest((int)((dPrecursorMass*1000)-(float)precMassAccuracy), keys);
             int endpoint = (int)((dPrecursorMass*1000)+(float)precMassAccuracy);
             for(int j =startpoint;;j++){
                 if(keys[j] > endpoint){
                     break;
                 }
                 HashMap<String,IndexedSequence>seqMap = dbMap.get(keys[j]);
                 if(seqMap == null){
                     System.out.println("");
                 }
                 Iterator g = seqMap.values().iterator();
                 while(g.hasNext()){
                     IndexedSequence iseq = (IndexedSequence) g.next();
                     List<Integer> index = iseq.getProteinIds();
                     List<String> proList = new ArrayList<>();
                     fasta = new Fasta(iseq.getWholeCleanSequence());
                     for(int h : index){
                         fasta.addDefList(blazmass.io.Fasta.getSequestLikeAccession(proteinMap.get(h)));
                     }
                     getPeptides(sr, iseq, fasta);
                 }
             }*/
           
            /*List<IndexedSequence> sequenceList = dbIndexer.getSequences(dPrecursorMass, precMassAccuracy);
                               //   sequenceList = dbIndexer.getSequences(dPrecursorMass, precMassAccuracy);
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
            //prs[i].start();*/
        } while (++i < numIsotopes);
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

    private void diffModSearch_orig(ProcessedPeakList ppl, SearchResult sr, Modifications m) throws IOException {
        //Modifications m = ppl.getModifications();
        // temp solution, will use m.getMassShift()
        //double massShift = m.getDiffNTermMod();
        double massShift = m.getMassShift();

        SearchParams params = ppl.getSearchParams();
        int numIsotopes = params.getNumIsotopicPeaks();
        double precMassAccuracy = params.getPrecursorTolerance();
        double acc = precMassAccuracy/1000000.0f;
        double prcMass = ppl.getPrecursorMass() - massShift;
        //MassCalculator mc = ppl.getMassCalculator();
            
        if(numIsotopes == 0) { // for low resolution data, traditional sequest like 
            double highLimit = Math.round(prcMass + ppl.getSearchParams().getHighPrecursorTolerance()/1000);
            //double lowLimit = Math.round(prcMass - ppl.getSearchParams().getHighPrecursorTolerance()/1000);
            double lowLimit = Math.round(prcMass - ppl.getSearchParams().getLowPrecursorTolerance()/1000);
            getPeptideHits(sr, highLimit, lowLimit, m);
        } else if(numIsotopes == 1) { // for deisotoped high resolution data
            double diffs = prcMass*acc;
            double highLimit = prcMass + diffs;
            double lowLimit = prcMass - diffs;
            getPeptideHits(sr, highLimit, lowLimit, m);
        } else { // numIsotopes >= 2, for non deisotoped high resolution data

            double [] highLimits = new double[numIsotopes]; 
            double [] lowLimits = new double[numIsotopes]; 
            double diffs = prcMass*acc*2;
            for(int i = 0; i < numIsotopes; i++) {
                lowLimits[i] = prcMass - diffs/2 - i*MassSpecConstants.MASSDIFFC12C13; 
                highLimits[i] = lowLimits[i] + diffs; 
            }
            getPeptideHits(sr, highLimits, lowLimits, numIsotopes, m);
        } 
           
    }
//    public void getPeptideHits(SearchResult sr, MassCalculator mc, double highLimit, double lowLimit) {
    // precMassAccuracy in ppm
    public synchronized SearchResult search(ProcessedPeakList ppl) throws IOException {
        /*-------SearchResult sr = new SearchResult(ppl);
        SearchParams params = ppl.getSearchParams();
        int numIsotopes = params.getNumIsotopicPeaks();
        double precMassAccuracy = params.getPrecursorTolerance();
        double prcMass = ppl.getPrecursorMass();
        precMassAccuracy =prcMass*precMassAccuracy/1000f;
        prcMass -= ppl.getSearchParams().getStaticTerminalMods();-----*/
        SearchResult sr = new SearchResult(ppl);
        SearchParams params = ppl.getSearchParams();
        //isDeCharged = ppl.isDeCharged();
        int numIsotopes = params.getNumIsotopicPeaks();
        double precMassAccuracy = params.getPrecursorTolerance();
        double acc = precMassAccuracy/1000000.0f;
        double prcMass = ppl.getPrecursorMass() - params.getStaticTerminalMods();
        //double minPrecMass = ppl.getSearchParams().getMinPrecursorMass();
        //double maxPrecMass = ppl.getSearchParams().getMaxPrecursorMass();

        /*
        //MassCalculator mc = ppl.getMassCalculator();
        //double mass = prcMass; 
          
        if(numIsotopes == 0) { // for low resolution, traditional sequest like  
            double highLimit = (prcMass + params.getHighPrecursorTolerance()/1000);
            double lowLimit = (prcMass - params.getLowPrecursorTolerance()/1000);
            getPeptideHits(sr, highLimit, lowLimit);
        } else if(numIsotopes == 1) { // for deisotoped high resolution data 
            double diffs = prcMass*acc;
            double highLimit = prcMass + diffs;
            double lowLimit = prcMass - diffs;
            getPeptideHits(sr, highLimit, lowLimit);
        } else { //for non deisotoped high resolution data
            double [] highLimits = new double[numIsotopes]; 
            double [] lowLimits = new double[numIsotopes]; 
            double diffs = prcMass*acc*2;
            for(int i = 0; i < numIsotopes; i++) {
                lowLimits[i] = prcMass - diffs/2 - i*MassSpecConstants.MASSDIFFC12C13; 
                highLimits[i] = lowLimits[i] + diffs; 
            }

            //getPeptideHits(sr, highLimits, lowLimits, numIsotopes);

           // addPeptideHit(result, f, templeft, rightEnd);


        }

        double minPrecMass = ppl.getSearchParams().getMinPrecursorMass();
        double maxPrecMass = ppl.getSearchParams().getMaxPrecursorMass();
        
        for(Iterator<Modifications> it = ppl.getModifications(); it.hasNext();) {
            Modifications m = it.next(); 
            if(m != null && m.getDiffModsShift() != 0 ) {
 //System.out.println("modification searches, massShfit: " + m.getDiffModsShift());
                double prcmass = ppl.getPrecursorMass() - m.getMassShift();
                if(prcmass >= minPrecMass && prcmass <= maxPrecMass) {
                    diffModSearch(ppl, sr, m);
                }
//System.out.println("finished diffmodSearch\n\n\n\n");
            }
        }


*/
        List<Integer> highLimits = new ArrayList<>();
        List<Integer> lowLimits = new ArrayList<>();
        
        if(numIsotopes == 0) { // for low resolution, traditional sequest like  
            int highLimit =(int) ((prcMass + params.getHighPrecursorTolerance()/1000)*1000);
            int lowLimit =(int) ((prcMass - params.getLowPrecursorTolerance()/1000)*1000);
            lowLimits.add(lowLimit);
            highLimits.add(highLimit);
        } else if(numIsotopes == 1) { // for deisotoped high resolution data 
            double diffs = prcMass*acc;
            int highLimit = (int)((prcMass + diffs)*1000);
            int lowLimit = (int)((prcMass - diffs)*1000);
            lowLimits.add(lowLimit);
            highLimits.add(highLimit);
        } else { //for non deisotoped high resolution data
            
            double diffs = prcMass*acc*2;
            for(int i = 0; i < numIsotopes; i++) {
                //lowLimits[i] = prcMass - diffs/2 - i*MassSpecConstants.MASSDIFFC12C13; 
                //highLimits[i] = lowLimits[i] + diffs; 
                
                lowLimits.add((int)((prcMass - diffs/2 - i*MassSpecConstants.MASSDIFFC12C13)*1000));
                highLimits.add((int)((lowLimits.get(i)+diffs*1000)));
            }
            
        }

        int i = 0;

        do{

            float dPrecursorMass = (float) (prcMass - i*MassSpecConstants.MASSDIFFC12C13);

            if(dPrecursorMass<startRange/1000 || dPrecursorMass>endRange/1000) continue;
            //System.out.println(""+dbMap.keySet().size());
             List<Integer> keyset = new ArrayList<>(dbMap.keySet());
             Integer [] keys = keyset.toArray(new Integer[keyset.size()]);
             Arrays.sort(keys);
             //List<IndexedSequence> sequenceList = new ArrayList<>();
             Fasta fasta = null;
             int startpoint =closest(lowLimits.get(i),keys);
             int endpoint  = highLimits.get(i);
             
             //int startpoint = closest((int)((dPrecursorMass+0.0005)*1000-(float)precMassAccuracy), keys);
             //int endpoint = (int)((dPrecursorMass+0.0005)*1000+(float)precMassAccuracy);
             //System.out.println("normal search\t"+dPrecursorMass);
             //System.out.println(keys[startpoint] +" \t"+endpoint);
             for(int j =startpoint;j<keys.length;j++){
                 if(keys[j] > keys[keys.length-1] || keys[j] > endpoint){
                     break;
                 }
               // System.out.println(""+j);
                 HashMap<String,IndexedSequence>seqMap = dbMap.get(keys[j]);
                 Iterator g = seqMap.values().iterator();
               
                 while(g.hasNext()){
                     IndexedSequence iseq = (IndexedSequence) g.next();
                     List<Integer> index = iseq.getProteinIds();
                   
                     
                    //List<String> proList = new ArrayList<>();
                     fasta = new Fasta(iseq.getWholeCleanSequence());
                     //System.out.println(""+iseq.getWholeCleanSequence());
                     for(int h : index){
                         fasta.addDefList(blazmass.io.Fasta.getSequestLikeAccession(proteinMap.get(h)));
                     }
                     getPeptides(sr, iseq, fasta);
                 }
             }
            
           // List<IndexedSequence> sequenceList = dbIndexer.getSequences(dPrecursorMass, (float) precMassAccuracy);
        /*   Fasta fasta = null;

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
            }*/
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
    
    public int closest(int of, Integer[] in) {
    int min = Integer.MAX_VALUE;
    int closest = of;

    for (int i=0;i<in.length;i++) {
        final int diff = Math.abs(in[i] - of);

        if (diff < min) {
            min = diff;
            closest = i;
        }
    }

    return closest;
}
    

    private void getPeptides(SearchResult sr, IndexedSequence seq, Fasta fasta) throws IOException {
        //    System.out.println("processing " + proteinIndex[i] + " for " + sr);
        //sr.addPeptideHit(new PeptideHit(getFasta(proteinIndex[i]), peptideStartIndex[i], peptideEndIndex[i]));

        int size = seq.getSequence().length();
        int start = 3;
        int end = start+size-1;
        sr.addPeptideHit(fasta, start, end);
        /*
        if(sr.getProcessedPeakList().isDeCharged()) {
            sr.addPeptideHit(new DeChargedPeptideHit(fasta, start, end));

        } else {

            //sr.addPeptideHit(new PeptideHit(fasta, start, end));
        }*/
    }


    public SearchResult search_orig(ProcessedPeakList ppl) throws IOException {
        SearchResult sr = new SearchResult(ppl);
        SearchParams params = ppl.getSearchParams();
        precMasses = ppl.getMassCalculator().getPrecMasses();
        //isDeCharged = ppl.isDeCharged();
        int numIsotopes = params.getNumIsotopicPeaks();
        double precMassAccuracy = params.getPrecursorTolerance();
        double acc = precMassAccuracy/1000000.0f;
        double prcMass = ppl.getPrecursorMass() - params.getStaticTerminalMods();
        //MassCalculator mc = ppl.getMassCalculator();

        //double mass = prcMass;

        if(numIsotopes == 0) { // for low resolution, traditional sequest like
            double highLimit = (prcMass + params.getHighPrecursorTolerance()/1000);
            double lowLimit = (prcMass - params.getLowPrecursorTolerance()/1000);
            getPeptideHits(sr, highLimit, lowLimit);
        } else if(numIsotopes == 1) { // for deisotoped high resolution data
            double diffs = prcMass*acc;
            double highLimit = prcMass + diffs;
            double lowLimit = prcMass - diffs;
            getPeptideHits(sr, highLimit, lowLimit);
        } else { //for non deisotoped high resolution data
            double [] highLimits = new double[numIsotopes];
            double [] lowLimits = new double[numIsotopes];
            double diffs = prcMass*acc*2;
            for(int i = 0; i < numIsotopes; i++) {
                lowLimits[i] = prcMass - diffs/2 - i*MassSpecConstants.MASSDIFFC12C13;
                highLimits[i] = lowLimits[i] + diffs;
            }
            getPeptideHits(sr, highLimits, lowLimits, numIsotopes);
        }
//          System.out.println("diff: " + diff + "\thighLimit: " + highLimit + "\tlowLimit: " + lowLimit);
        /* // old
        Modifications m = ppl.getModifications();
        if(m != null && m.getDiffModsShift() != 0 ) {
            diffModSearch(ppl, sr);
            //System.out.println("finished diffmodSearch\n\n\n\n");
        }
        */
        double minPrecMass = ppl.getSearchParams().getMinPrecursorMass();
        double maxPrecMass = ppl.getSearchParams().getMaxPrecursorMass();

        for(Iterator<Modifications> it = ppl.getModifications(); it.hasNext();) {
            Modifications m = it.next();
            if(m != null && m.getDiffModsShift() != 0 ) {
                //System.out.println("modification searches, massShfit: " + m.getDiffModsShift());
                double prcmass = ppl.getPrecursorMass() - m.getMassShift();
                if(prcmass >= minPrecMass && prcmass <= maxPrecMass) {
                    diffModSearch(ppl, sr, m);
                }
//System.out.println("finished diffmodSearch\n\n\n\n");
            }
        }
        return sr;
    }


    // working on making the isotopic peak search faster
    protected void getPeptideHits(SearchResult result,  
                             double [] highLimits, double [] lowLimits, int numPeaks) {
        int num = 0; 
        double lowLimit = lowLimits[numPeaks-1];
        double highLimit = highLimits[0];

        for(Fasta f : sequences) {
           //num++;
        //System.out.println("numproteins: " + num);
            int lastIndex = f.getLength() - 1;
            byte [] seq = f.getSequenceAsBytes();
         // mass of H2O need to be added 
         // because the mass of the residues does not count H2O   
            double mass = MassSpecConstants.MASSH3O;
            double tempmass = 0;
            int templeft = 0;
            int rightEnd = -1;
            int leftEnd = 0; 
             
            while(leftEnd < lastIndex && rightEnd < lastIndex) { 
                while (mass < lowLimit && rightEnd < lastIndex) {
                    mass += precMasses[seq[++rightEnd]]; 
                }
                if (mass < lowLimit) {
                    //return;
                    break;
                }
                tempmass = mass;
                templeft = leftEnd;
                while(mass == tempmass) {
                    mass -= precMasses[seq[leftEnd++]]; 
                }
                // now the mass should be >= lowLimit
                if (tempmass <= highLimit) {  
                    for(int i = 0; i < numPeaks; i++) {
                        if(tempmass >= lowLimits[i] && tempmass <= highLimits[i]) {
//if(mass > 2040.070 && mass < 2040.072)
//System.out.println("Found!!!! mass: " + mass);
                            //result.addPeptideHit(new PeptideHit(f, templeft, rightEnd));
                            //result.addPeptideHit(createPeptideHit(new Peptide(f, templeft, rightEnd)));
                            //result.addPeptideHit(f, templeft, rightEnd);
                            addPeptideHit(result, f, templeft, rightEnd);
                            break;
                        }
                    }
                }

                // check the last residue in case removed the left the mass is within tolerance
                // because this is the last peptide, so no need to adjust the leftEnd and mass
                if(rightEnd == lastIndex && mass <= highLimit && mass >= lowLimit) {
                    for(int i = 0; i < numPeaks; i++) {
                        //if(tempmass >= lowLimits[i] && tempmass <= highLimits[i]) {
                        if(mass >= lowLimits[i] && mass <= highLimits[i]) {
                            addPeptideHit(result, f, leftEnd, rightEnd);
                            break;
                        }
                    }
                }
            }
        }
    }
    
    public SearchResult topDownSearch(ProcessedPeakList ppl) throws IOException {
        SearchResult sr = new SearchResult(ppl);
        SearchParams params = ppl.getSearchParams();
        for(Fasta f : sequences) {
            sr.addPeptideHit(new PeptideHit(f, 0, f.getLength() - 1));
        }
        //MassCalculator mc = ppl.getMassCalculator();
 
        //double mass = prcMass; 
          
        return sr;
    }
    protected void countIsolucineProteins() {
        int numI = 0;
        int numNoI = 0;
        for(Fasta f: sequences) {
            if(f.getSequence().indexOf("I") == -1) {
                numNoI++;
                System.out.println(f.getLength());
            } else {
                numI++;
            }
        }

        System.out.println("NumWithI " + numI + "\tNumNoI: " + numNoI);
    }
    protected void getTrypticPeptideDristribution() {
        int maxNumResidue = 5;
        int maxLength = 41;
        int[][] freq = new int[maxNumResidue][maxLength];
        int totalI = 0;
        int totalResidue = 0;
        for(Fasta f : sequences) {
            int lastIndex = f.getLength() - 1;
            byte [] seq = f.getSequenceAsBytes();
            totalResidue += seq.length;
            for(byte b : seq) {
                if(b == 'I') {
                    totalI++;
                }
            }
            int rightEnd = 6;
            int leftEnd = 0; 
            while(leftEnd < lastIndex && rightEnd <= lastIndex) { 
                if(leftEnd == 0 || isTryptic(seq[leftEnd-1])) {
                    rightEnd = leftEnd + 6;
                    while(rightEnd <= lastIndex && rightEnd - leftEnd < maxLength-1) {
                        if(rightEnd == lastIndex || isTryptic(seq[rightEnd])) {
                            int numI = 0;
                            int numMisCleavage = 0;
                            for(int i = leftEnd; i <= rightEnd; i++) {
                                byte b = seq[i];
                                if(b == 'I') {
                                    numI++;
                                }
                                if(b == 'R' || b == 'K') {
                                    numMisCleavage++;
                                }
                                
                            }
                            if(numI > 4) {
                                numI = 4;
                            }
                            if(numMisCleavage < 2) {
                                freq[numI][rightEnd-leftEnd+1]++;
                            }
                        }
                        rightEnd++;
                    }
                    
                }   
                leftEnd++;
            }
        }
        System.out.println("numResidue\t" + totalResidue + "\tnumI\t" + totalI);
        System.out.print("numI/length\t"); 
        for(int i = 7; i < maxLength; i++) {
            System.out.print(i + "\t");
        }
        System.out.println();
        for(int i = 0; i < maxNumResidue; i++) {
            StringBuffer sb = new StringBuffer(500);
            sb.append(i + "\t");
            for(int j = 7; j < maxLength; j++) {
                sb.append(freq[i][j]+ "\t"); 
            }
            System.out.println(sb);
        }
    }
 
    
    private boolean isTryptic(byte b) {
        return b == 'R' || b == 'K';    
    }
    // assuming the mass tolerance is smaller than the mass of any AA residue
    public static void main(String [] args) throws Exception {
        //ProteinDatabase pd = new ProteinDatabase(args[0]);
         //pd.getTrypticPeptideDristribution();
        // pd.countIsolucineProteins();
        /*
        for(Iterator<Fasta> it = pd.getFastas(); it.hasNext();) {
            Fasta f = it.next();
            String ac = f.getAccession();
            System.out.println("Seq from f: \n" + f.getSequence());
            System.out.println("Seq from hash: \n" + pd.accession2Fasta(ac).getSequence());
        } 
        */ 
    }
}
