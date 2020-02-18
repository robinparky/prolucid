
package blazmass.dbindex;

import blazmass.AssignMass;
import blazmass.Constants;
import blazmass.Enzyme;
import blazmass.dbindex.DBIndexStore.FilterResult;
import java.io.*;
import java.util.*;
import blazmass.io.Fasta;
import blazmass.io.FastaReader;
import blazmass.io.SearchParamReader;
import blazmass.io.SearchParams;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * Supports storing and retrieving proteins and sequences to/from index
 *
 * There can be different underlying store index implementations, implementing
 * DBIndexStore interface
 *
 * Supports modes: indexing only (store), search from indexed store, and search
 * without a store (cut on the fly) with temporary in-memory store
 *
 */
public class DBIndexerNoSQL {

    /**
     * Different operational modes supported
     */
    public enum IndexerMode {

        INDEX, SEARCH_INDEXED, SEARCH_UNINDEXED, 
    };
    /**
     * Different operational modes supported
     */
    public enum IndexType {

        INDEX_NORMAL {
            @Override
            public String toString() {
                return "Normal index (best for small and medium db)";
            }
            
        }, 
        INDEX_LARGE {
            @Override
            public String toString() {
                return "Large index (best for large db)";
            }
            
        },
    };
    private IndexerMode mode;
    private IndexType indexType;
   // private DBIndexStore indexStore;
    private blazmass.io.SearchParams sparam;
    //private float massTolerance;
    private boolean inited;
    private String indexName; //index name that is params-specific
    private long protNum; //protein number in sequence, starting at 1
    private static final Logger logger = Logger.getLogger(DBIndexer.class.getName());
    private static final int MAX_SEQ_LENGTH = 10000;
    private HashMap<Integer, HashMap<String,IndexedSequence>> dbMap = new HashMap<>();
    private HashMap<Integer,String> proteinMap = new HashMap<>();
    private int startRange;
    private int endRange;

    public static void main(String[] args) {

        if (args.length == 0) {
            System.out.println("Need a directory name");
            System.exit(1);
        }

    
    }

    public HashMap<Integer, String> getProteinMap() {
        return proteinMap;
    }

    public void setProteinMap(HashMap<Integer, String> proteinMap) {
        this.proteinMap = proteinMap;
    }
    



    public DBIndexerNoSQL(SearchParams sparam,int start,int end) throws DBIndexerException {
        this.sparam = sparam;

        boolean inMemoryIndex = true;
        this.startRange=start;
        this.endRange=end;
        init();
        run();

        //this.indexStore = new DBIndexStoreSQLiteMult(this.sparam, true);
        
        
        this.mode = IndexerMode.SEARCH_INDEXED;

        this.inited = false;
        protNum = -1;


        if (inited) {
            throw new IllegalStateException("Already inited");
        }

        //reset prot number
        protNum = -1;

        final String dbName = sparam.getDatabaseName();

        File dbFile = new File(dbName);
        if (!dbFile.exists() || !dbFile.canRead()) {
            throw new DBIndexerException("Cannot read (and index) the Fasta database file: " + dbName);
        }

        //get param specific name for the inde0x
        indexName = this.getFullIndexFileNameWithNoIndex();
        // indexName = this.getFullIndexFileName();

    

    }

    /**
     * Cut a Fasta protein sequence according to params spec and index the
     * protein and generated sequences
     *
     * Should be called once per unique Fasta sequence
     *
     * @param fasta fasta to index
     * @throws IOException
     */
    private void cutSeq(final Fasta fasta, int start, int end, int proteinIndex) throws IOException {
        final String protAccession = fasta.getSequestLikeAccession();
        final String protSeq = fasta.getSequence();
        fasta.getAccession();
        
        this.proteinMap.put(proteinIndex, fasta.getDefline());
        //cutSeq(protAccession, protSeq,start,end, proteinIndex);
        cutSeq(protAccession, protSeq,start, end,proteinIndex);
        //System.out.println("==========");
        //System.out.println("called");
    }

    /**
     * Cut a Fasta protein sequence according to params spec and index the
     * protein and generated sequences
     *
     * Should be called once per unique Fasta sequence
     *
     * @param fasta fasta to index
     * @throws IOException
     */
    
    
    
    
    
    
    
private void cutSeq(final String protAccession,final String protSeq,int startRange, int endRange, int proteinIndex) throws IOException {

        //Enzyme enz = sparam.getEnzyme();
        
        final int length = protSeq.length();
//
HashMap<String,IndexedSequence> obj = null;
        //AssignMass aMass = AssignMass.getInstance(true);

        final char[] pepSeq = new char[10000]; //max seq length
        int curSeqI = 0;

        final int maxMissedCleavages = sparam.getMaxMissedCleavages();
        int maxIntCleavage = sparam.getMaxInternalCleavageSites();
        HashMap<String,IndexedSequence> newobj;
        try {
           // long proteinId = indexStore.addProteinDef(++protNum, protAccession, protSeq);
//	    System.out.println(fasta.getSequestLikeAccession());
            //System.out.println(fasta.getDefline());
            
            for (int start = 0; start < length; ++start) {
                int end = start;

                //clear the preallocated seq byte array
                //Arrays.fill(seq, 0, curSeqI > 0?curSeqI-1:0, (byte) 0); //no need, we copy up to curSeqI nowu
                curSeqI = 0;

                //float precMass = Constants.H2O_PROTON_SCALED_DOWN;
                float precMass = Constants.H2O_PROTON;
                precMass += AssignMass.getcTerm();
                precMass += AssignMass.getnTerm();
                

                // System.out.println("===>>" + precMass + "\t" + Constants.MAX_PRECURSOR);
                //System.out.println("==" + j + " " + length + " " + (j < length));

                //int testC=0;
                int pepSize = 0;

                int intMisCleavageCount = -1;


                //while (precMass <= Constants.MAX_PRECURSOR_MASS && end < length) {
                while (precMass <= sparam.getMaxPrecursorMass() && end < length) 
                {
                    pepSize++;

                    final char curIon = protSeq.charAt(end);
                    pepSeq[curSeqI++] = curIon;
                    precMass += AssignMass.getMass(curIon);
                    

                  if (Enzyme.isEnzyme(protSeq.charAt(end))) {
                        intMisCleavageCount++;
                    }

                    final int cleavageStatus = Enzyme.checkCleavage(protSeq, start, end, sparam.getEnzymeNocutResidues());

                    //System.out.println("---" + String.valueOf(Arrays.copyOf(pepSeq, curSeqI)) + "\t" + intMisCleavageCount + "\t" + maxIntCleavage);

                    if (intMisCleavageCount > maxIntCleavage) {
                        break;
                    }

                    
                    //if (precMass > Constants.MAX_PRECURSOR_MASS) {
                   // if (precMass > sparam.getMaxPrecursorMass()) {
                    //    break;
                    //}
                  
                  

                    
                    //if (precMass > Constants.MAX_PRECURSOR_MASS) {
                    //System.out.println(""+precMass);
                 
                    //int newprecur = (int)((precMass+0.0005)*1000);
                    if (precMass>endRange/1000) {
                        //++end;
                      // continue;
                       break;
                    }
                    
                   // System.out.println(""+precMass);
                    

                    
                    
                    if (pepSize >= Constants.MIN_PEP_LENGTH && precMass >= startRange/1000) { //Constants.MIN_PRECURSOR ) {
                        
                        if (Enzyme.isEnzyme(protSeq.charAt(end))) {
                        intMisCleavageCount++;
                    }
                        //final int cleavageStatus = Enzyme.checkCleavage(protSeq, start, end, sparam.getEnzymeNocutResidues());
                        if (cleavageStatus >= maxMissedCleavages) {
                        //if (cleavageStatus == 2) {
                            //qualifies based on params

                            final String peptideSeqString = String.valueOf(Arrays.copyOf(pepSeq, curSeqI));
                           // if(peptideSeqString.equals("DDRTGYDDDLRRRAAHAS")){
                             //   System.out.println("");
                             // }
                            
                            
                           // if(peptideSeqString.equals("MSPAVNGLDTLAVLGICTPLK")){
                           //     System.out.println("");
                         //  }
                           // if(peptideSeqString.equals("LAEQAE")){
                            //    System.out.println("");
                          //  }
                            //System.out.println(peptideSeqString);
                          //  count++;
                           // if(peptideSeqString.equals("ABCDEK")){
                           //     System.out.println("");
                           // }
                            //System.out.println(""+peptideSeqString);
                            //check if index will accept it
                          //  final FilterResult filterResult = indexStore.filterSequence(precMass, peptideSeqString);

                         //   if (filterResult.equals(FilterResult.SKIP_PROTEIN_START) ) {
                               //bail out earlier as we are no longer interested in this protein starting at start
                         //       break; //move to new start position
                         //   }
                         //   else if (filterResult.equals(FilterResult.INCLUDE) ) {
                                final int resLeftI = start >= Constants.MAX_INDEX_RESIDUE_LEN ? start - Constants.MAX_INDEX_RESIDUE_LEN : 0;
                                final int resLeftLen = Math.min(Constants.MAX_INDEX_RESIDUE_LEN, start);
                                StringBuilder sbLeft = new StringBuilder(Constants.MAX_INDEX_RESIDUE_LEN);
                                for (int ii = 0; ii < resLeftLen; ++ii) {
                                    sbLeft.append(protSeq.charAt(ii + resLeftI));
                                }
                                final int resRightI = end + 1;
                                final int resRightLen = Math.min(Constants.MAX_INDEX_RESIDUE_LEN, length - end - 1);
                                StringBuilder sbRight = new StringBuilder(Constants.MAX_INDEX_RESIDUE_LEN);
                                if (resRightI < length) {
                                    for (int jj = 0; jj < resRightLen; ++jj) {
                                        sbRight.append(protSeq.charAt(jj + resRightI));
                                    }
                                }


                                //add -- markers to fill Constants.MAX_INDEX_RESIDUE_LEN length
                                final int lLen = sbLeft.length();
                                for (int c = 0; c < Constants.MAX_INDEX_RESIDUE_LEN - lLen; ++c) {
                                    sbLeft.insert(0, '-');
                                }
                                final int rLen = sbRight.length();
                                for (int c = 0; c < Constants.MAX_INDEX_RESIDUE_LEN - rLen; ++c) {
                                    sbRight.append('-');
                                }

                                final String resLeft = sbLeft.toString();
                                final String resRight = sbRight.toString();
                                int key = (int)((precMass+0.0005)*1000);
                                 //newprecur = (int)((precMass+0.0005)*1000);
                    //if (  newprecur <= startRange || newprecur>=endRange) {
                       // ++end;
                    //    break;
                       //break;
                   // }
                           // if (  newprecur <= startRange || newprecur>=endRange) {
                           // }else{
                               obj = this.dbMap.get(key);
                                
                                if(obj != null){
                                    
                                   IndexedSequence seq = obj.get(peptideSeqString);
                                   if(seq != null){

                                       List<Integer>proIds = seq.getProteinIds();
                                       //StringBuffer proList = seq.getProList();
                                       if(proIds.contains(proteinIndex)){                                           
                                       }else{
                                           proIds.add(proteinIndex);
                                       //if(proList.indexOf(","+proteinIndex+",") > -1 || proList.indexOf(proteinIndex+",") > -1 || proList.indexOf(","+proteinIndex) > -1)
                                       //      proList.append(proteinIndex+",");
                                       }
                                   }
                                   else{

                                       List<Integer> proList = new ArrayList<>();
                                       proList.add(proteinIndex);
                                       //StringBuffer proList = new StringBuffer();
                                       //proList.append(proteinIndex+",");
                                       IndexedSequence newseq = new IndexedSequence(precMass, peptideSeqString, resLeft, resRight, resLeftLen, proList);
                                       obj.put(peptideSeqString, newseq);
                                   }
                                }
                                else {
                                    newobj = new HashMap<>();
                                   List<Integer> prolist = new ArrayList<>();
                                   prolist.add(proteinIndex);
                                  // StringBuffer proList = new StringBuffer();
                                     //  proList.append(proteinIndex+",");
                                   IndexedSequence seq = new IndexedSequence(precMass, peptideSeqString, resLeft, resRight, resLeftLen,prolist);
                                   newobj.put(peptideSeqString,seq);
                                   dbMap.put(key, newobj);
                                }
                            }

                        }
                   // }


                    if (intMisCleavageCount == maxIntCleavage) {
                        break;
                    }

                    ++end;

                }
//
            }
        } catch (Exception e) {
            logger.log(Level.SEVERE, "Error writing sequence to db index store, ", e);
        }

    }



    
    /**
     * Initialize the indexer
     *
     * @throws IOException
     */
    public void init() throws DBIndexerException {
        if (inited) {
            throw new IllegalStateException("Already inited");
        }         

        //reset prot number
        protNum = -1;

        final String dbName = sparam.getDatabaseName();

        File dbFile = new File(dbName);
        if (!dbFile.exists() || !dbFile.canRead()) {
            throw new DBIndexerException("Cannot read (and index) the Fasta database file: " + dbName);
        }

        //get param specific name for the inde0x
        indexName = this.getFullIndexFileNameWithNoIndex();
       // indexName = this.getFullIndexFileName();

    }

    /**
     * Set protein cache if storage supports it Useful if database has already
     * been indexed Otherwise, the cache can be populated while indexing is done
     */
    private void setProteinCache() {


        final ProteinCache protCache = ProteinCache.getInstance();

        //ensure only 1 thread at time enters this
        synchronized (protCache) {

            logger.log(Level.INFO, "Populating protein cache");
            if (protCache.isPopulated() == false) {
                logger.log(Level.INFO, "Initializing protein cache");
                InputStream fis = null;
                try {
                    fis = new FileInputStream(sparam.getDatabaseName());
                } catch (FileNotFoundException ex) {
                    logger.log(Level.SEVERE, "Could not set protein cache", ex);
                    return;
                }
                Iterator<Fasta> itr = null;
                try {
                    itr = FastaReader.getFastas(fis);
                } catch (IOException ex) {
                    logger.log(Level.SEVERE, "Could not set protein cache", ex);
                    return;
                }

                while (itr.hasNext()) {
                    Fasta fasta = itr.next();
                    protCache.addProtein(fasta.getSequestLikeAccession(), fasta.getSequence());
                }
                
                logger.log(Level.INFO, "Done initializing protein cache");
            } else {
                //logger.log(Level.INFO, "Protein cache already populated, reusing.");
            }
        }

    }

    /**
     * Indexes the database if database index does not exists
     *
     * @throws IOException
     */
    public void run() throws DBIndexerException {
        /*if (!inited) {
            throw new IllegalStateException("Not initialized.");
        }*/

        //if search unindexed, do nothing, will use protein cache to cutSeq on the fly
       /* if (mode.equals(IndexerMode.SEARCH_UNINDEXED)) {
            return;
        }*/


        FileInputStream fis = null; //fasta reader

        //setup status writer
        String statusFilePath = blazmass.Util.getFileBaseName(sparam.getDatabaseName()) + "log";
        FileWriter statusWriter = null;
        String totalProteins = null;
        int indexedProteins = 0;
        try {
            statusWriter = new FileWriter(statusFilePath);
            totalProteins = Integer.toString(FastaReader.getNumberFastas(sparam.getDatabaseName()));
        } catch (IOException ex) {
            logger.log(Level.SEVERE, "Error initializing index progress writer for file path: " + statusFilePath, ex);
        }

        //start indexing
        try {
            //index, set protein cache on the fly
       //     indexStore.startAddSeq();

            fis = new FileInputStream(sparam.getDatabaseName());

            //create prot cache, in case searcher is kicked off in the same process
            //after indexing is done
            ProteinCache protCache = ProteinCache.getInstance();
  //          indexStore.setProteinCache(protCache);

            int proteinId=0;
            for (Iterator<Fasta> itr = FastaReader.getFastas(fis); itr.hasNext();) {
                Fasta fasta = itr.next();
                protCache.addProtein(fasta.getSequestLikeAccession(), fasta.getSequence());
                cutSeq(fasta,startRange,endRange,proteinId);
//                indexStore.lastBuffertoDatabase();
                      

              /*  ++indexedProteins;
                if (statusWriter != null) {
                    statusWriter.append(totalProteins).append("\t").append(Integer.toString(indexedProteins)).append("\n");
                    if (indexedProteins % 100 == 0) {
                        statusWriter.flush();
                    }
                }*/
                
                proteinId++;

           //    System.out.println(""+proteinId+"\t"+dbMap.size());

               
            }
          /* Iterator k = dbMap.keySet().iterator();
            while(k.hasNext()){
               // System.out.println("===========================");
             HashMap<String,IndexedSequence> hmap = this.dbMap.get(k.next());
                            
                            Iterator i = hmap.keySet().iterator();
                           // System.out.println("------------------------");
                            while(i.hasNext()){
                               
                                IndexedSequence iseq = hmap.get(i.next().toString());
                                System.out.println(iseq.getSequence());
                                 //System.out.print("proteinID:");
                                //for(int j :iseq.getProteinIds()){
                                //   System.out.print(j+"\t");
                              //  }
                               // System.out.println("");
                            }
            }*/
            //System.out.println("");
//               System.out.print("Printing the last buffer....");
//               indexStore.lastBuffertoDatabase();

        } catch (IOException ex) {
            logger.log(Level.SEVERE, "Error initializing and adding sequences", ex);
            throw new DBIndexerException("Error initializing and adding sequences", ex);
        }

    }

  /*  private List<IndexedSequence> cutAndSearch(List<MassRange> massRanges) {
        if (!mode.equals(IndexerMode.SEARCH_UNINDEXED)) {
            throw new RuntimeException("Cut and search only supported for SEARCH_UNINDEXED mode !");
        }

        //set up in memory temporary "index" that does the filtering
        ((MassRangeFilteringIndex) indexStore).init(massRanges);

        //reset prot number
        protNum = -1;

        //cut sequences in prot cache and check if qualify
/*
        final ProteinCache protCache = ProteinCache.getInstance();

        final int numProteins = protCache.getNumberProteins();

        for (int prot = 0; prot < numProteins; ++prot) {
            final String protDef = protCache.getProteinDef(prot);
            final String protSequence = protCache.getProteinSequence(prot);

            try {
                //stores sequences that match only in our temporary "index"
                cutSeq(protDef, protSequence);
            } catch (IOException ex) {
                logger.log(Level.SEVERE, "Error cutting sequence in no-index mode: " + protSequence, ex);
            }

        }

        List<IndexedSequence> sequences = null;
        try {
            sequences = indexStore.getSequences(massRanges);
            
        } catch (DBIndexStoreException ex) {
            logger.log(Level.SEVERE, "Error getting sequences from in-memory filtering index", ex);
        }


        return sequences;

    }*/

    /**
     * Get list of sequences in index for the mass, and using tolerance in
     * SearchParams. Mass dependant tolerance is already scaled/precalculated
     *
     * Requires call to init() and run() first, which might perform indexing, if
     * index for current set of sequest params does not exist
     *
     * @param precursorMass mass to select by, in ppm
     * @param massTolerance mass tolerance, already calculated for that mass as
     * it is mass-dependant
     * @return list of matching sequences
     */
   /* public List<IndexedSequence> getSequences(float precursorMass, float massTolerance) {
 
            try {
                return indexStore.getSequences(precursorMass, massTolerance);
            } catch (DBIndexStoreException ex) {
                logger.log(Level.SEVERE, "Error getting sequences for the mass.", ex);
                return null;
            }
        
    }*/

    /**
     * Get list of sequences in index for the ranges specified wit mass .
     * Requires call to init() and run() first, which might perform indexing, if
     * index for current set of sequest params does not exist
     *
     * @param massRanges mass ranges to query
     * @return list of matching sequences
     */
    /*public List<IndexedSequence> getSequences(List<MassRange> massRanges) {
        List<IndexedSequence> sequences = null;
        long t1 = System.currentTimeMillis();

            try {

                sequences = indexStore.getSequences(massRanges);
            } catch (DBIndexStoreException ex) {
                logger.log(Level.SEVERE, "Error getting sequences for the mass ranges: " + massRanges.toString(), ex);
                return null;
            }
        
        //logger.log(Level.INFO, "getSequences(): " + sequences.size());
        return sequences;
    }

    /**
     * Get proteins associated with the sequence from the index. Requires call
     * to init() and run() first, which might perform indexing, if index for
     * current set of sequest params does not exist
     *
     * @param seq peptide sequence
     * @return list of indexed protein objects associated with the sequence
     */
    /*public List<IndexedProtein> getProteins(IndexedSequence seq) {
        try {
            return indexStore.getProteins(seq);
        } catch (DBIndexStoreException ex) {
            logger.log(Level.SEVERE, "Error getting protein for the sequence", ex);
            return null;
        }
    }
    */
    /**
     * Get filename for index file base name The filename is specific to set of
     * params that affect the index
     *
     * @return file name that includes base name and unique token for the set of
     * params that affect index
     */
    
    private String getFullIndexFileName() {
        return sparam.getFullIndexFileName();
    }

    private String getFullIndexFileNameWithNoIndex() {
        return sparam.getFullFileNameWithNoIndex();
    }

    /*
     * 
     private String getFullIndexFileName(String baseName) {
     String uniqueIndexName = baseName + "_";

     //generate a unique string based on current params that affect the index
     final StringBuilder uniqueParams = new StringBuilder();
     //uniqueParams.append(sparam.getEnzyme().toString());
     //uniqueParams.append(sparam.getEnzymeNumber());                
     uniqueParams.append(sparam.getEnzymeOffset());
     uniqueParams.append(sparam.getEnzymeResidues());
     uniqueParams.append(sparam.getEnzymeNocutResidues());

     uniqueParams.append("\nCleav: ");
     uniqueParams.append(sparam.getMaxInternalCleavageSites());
     uniqueParams.append(sparam.getMaxMissedCleavages());

     uniqueParams.append(sparam.getMaxNumDiffMod());
     uniqueParams.append("\nMods:");
     for (final ModResidue mod : sparam.getModList()) {
     uniqueParams.append(mod.toString()).append(" ");
     }
     uniqueParams.append("\nMods groups:");
     for (final List<Float> modGroupList : sparam.getModGroupList()) {
     for (final Float f : modGroupList) {
     uniqueParams.append(f).append(" ");
     }
     }

     final String uniqueParamsStr = uniqueParams.toString();

     //logger.log(Level.INFO, "Unique params: " + uniqueParamsStr);

     }*/
    public void close() {
    }

    public boolean isInited() {
        return inited;
    }

    public HashMap<Integer, HashMap<String, IndexedSequence>> getDbMap() {
        return dbMap;
    }

    public void setDbMap(HashMap<Integer, HashMap<String, IndexedSequence>> dbMap) {
        this.dbMap = dbMap;
    }

    
    
}
