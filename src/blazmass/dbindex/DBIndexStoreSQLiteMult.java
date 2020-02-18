package blazmass.dbindex;

import blazmass.Constants;
import blazmass.dbindex.DBIndexer.IndexerMode;
import blazmass.io.SearchParams;
import java.io.File;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 * Using multiple databases to divide index into multiple manageable chunks
 *
 * @author Adam
 *
 */
public final class DBIndexStoreSQLiteMult implements DBIndexStore {

    private long totalSeqCount = 0;
    private boolean inited = false;
    private ProteinCache proteinCache;
    private String dbPathBase = null;
    private static final Logger logger = Logger.getLogger(DBIndexStoreSQLiteMult.class.getName());
    private final DBIndexStoreSQLiteByteIndexMerge[] buckets;
    private static final String IDX_SUFFIX = ".idx";
    //we divide into buckets based on this
    private static final int MAX_MASS = (int) Constants.MAX_PRECURSOR_MASS + 400;
    //adjust number of buckets to make sure no single bucket ends up with more than 100million sequences
    private static  int NUM_BUCKETS = 0; //will override in constr
    private static  int BUCKET_MASS_RANGE = 0; //will override in constr
    private static final int MASS_GROUP_FACTOR = 1; //1 ppm per row
    private static final int PPM_PER_ENTRY = 1;
    //commit the buffered row and clear the buffer after this many seq
    private static final int COMMIT_SEQUENCES = 20000;
    private static final int BYTE_PER_SEQUENCE = 4 * 4; //float mass, int offset, length, protein id
    private static final int MAX_MASS_RANGES = 24; //will use optimzed statements for first 24, and by hand for above 24


    private SearchParams sparam;
    private boolean inMemoryIndex = false;
    private IndexerMode indexerMode;
    

    public DBIndexStoreSQLiteMult(SearchParams sparam, boolean inMemoryIndex) { 
        inited = false;
        this.sparam = sparam;
        totalSeqCount = 0;
        
        this.inMemoryIndex = inMemoryIndex;
        
        NUM_BUCKETS = sparam.getIndexFactor();
        BUCKET_MASS_RANGE = MAX_MASS / NUM_BUCKETS;
        
        
        buckets = new DBIndexStoreSQLiteByteIndexMerge[NUM_BUCKETS];
        
        if (inMemoryIndex) {
          //  throw new IllegalArgumentException("In-memory index not supported");
        }
    }

    @Override
    public boolean indexExists() throws DBIndexStoreException {
        if (!inited) {
            throw new DBIndexStoreException("Not intialized");
        }

        File indexDir = new File(dbPathBase);
        if (!indexDir.exists()) {
            return false;
        }

        File[] indexFiles = indexDir.listFiles();
        if (indexFiles.length == 0) {
            return false;
        }

        return this.getNumberSequences() > 0;

    }

    @Override
    public void init(String databaseID) throws DBIndexStoreException {
        if (databaseID == null || databaseID.equals("")) {
            throw new DBIndexStoreException("Index path is missing, cannot initialize the indexer.");
        }

        if (inited) {
            throw new DBIndexStoreException("Already intialized");
        }

        File dbBase = new File(databaseID).getAbsoluteFile();
        String baseName = dbBase.getName();
        if (!baseName.endsWith(IDX_SUFFIX)) {
            baseName = baseName + IDX_SUFFIX;
        }
        File parentDir = dbBase.getParentFile();

        File indexDir = new File(parentDir.getAbsolutePath() + File.separator + baseName);
        if (indexDir.exists() && !indexDir.isDirectory()) {
            logger.log(Level.INFO, "Trying to delete old index file: " + indexDir.getAbsolutePath());
            indexDir.delete();

        }

        try {
            indexDir.mkdir();
        } catch (SecurityException e) {
        }

        this.dbPathBase = indexDir.getAbsolutePath();


        logger.log(Level.INFO, "Using database index dir: " + dbPathBase);

        if (!indexDir.exists()) {
            throw new DBIndexStoreException("Index dir does not exist: " + indexDir + ", cannot initialize the indexer.");
        }

        //initialize buckets
        for (int i = 0; i < NUM_BUCKETS; ++i) {
            //buckets[i] = new DBIndexStoreSQLiteByte(i); //version with merge during search
            buckets[i] = new DBIndexStoreSQLiteByteIndexMerge(inMemoryIndex, i); //version with merge during index
            String bucketId = dbPathBase + File.separator + i + IDX_SUFFIX;
            try {
                buckets[i].init(bucketId);
            } catch (DBIndexStoreException e) {
                logger.log(Level.SEVERE, "Error initializing bucket " + i, e);
                throw e;
            }
        }


        inited = true;


    }

    @Override
    public void startAddSeq() throws DBIndexStoreException {
        if (!inited) {
            throw new DBIndexStoreException("Indexer is not initialized");
        }

        logger.log(Level.INFO, "Starting adding sequences");

        for (int i = 0; i < NUM_BUCKETS; ++i) {
            buckets[i].startAddSeq();
        }
    }

    @Override
    public void stopAddSeq() throws DBIndexStoreException {
        if (!inited) {
            throw new DBIndexStoreException("Indexer is not initialized");
        }

        logger.log(Level.INFO, "Finalizing adding sequences");

        logger.log(Level.INFO, "Will flush sequences from cache, merge sequences, index and compact the index database.");

        for (int i = 0; i < NUM_BUCKETS; ++i) {
            buckets[i].stopAddSeq();
        }

        logger.log(Level.INFO, "Total sequences added: " + totalSeqCount);
    }

    @Override
    public long getNumberSequences() throws DBIndexStoreException {
        //NOTE, this only returns total number of mass entries, not really useful anymore
        //we don't know num of sequenes in index at search time, because we would have to readit all
        long total = 0;
        for (int i = 0; i < NUM_BUCKETS; ++i) {
            total += buckets[i].getNumberSequences();
        }
        return total;
    }

    private int getBucketForMass(float precMass) {
        return (int) precMass / BUCKET_MASS_RANGE;
    }

    /**
     * get min, max bucket inclusive for the orig. mass range
     *
     * @param minMass
     * @param maxMass
     * @return
     */
    private int[] getBucketsForMassRange(float minMass, float maxMass) {
        if (maxMass > MAX_MASS) {
            logger.log(Level.SEVERE, "Trying to get a sequence for more than supported mass: " + maxMass);
        }

        int[] ret = new int[2];
        ret[0] = (int) minMass / BUCKET_MASS_RANGE;
        ret[1] = (int) maxMass / BUCKET_MASS_RANGE;

        //if (ret[0] != ret[1]) {
        //checking how often we are between buckets
        //logger.log(Level.INFO, "Between buckets: masses: " + minMass + ", " + maxMass + ", buckets: "
        //      + ret[0] + ", " + ret[1]);
        //}
        return ret;
    }

    @Override
    public FilterResult filterSequence(float precMass, String sequence) {
        //no filtering, we add every sequence
        return FilterResult.INCLUDE;
    }

    @Override
    public void addSequence(float precMass, int seqOffset, int seqLength, String sequence, String resLeft, String resRight, long proteinId) throws DBIndexStoreException {
        if (!inited) {
            throw new DBIndexStoreException("Indexer is not initialized");
        }

        totalSeqCount++;
        if (totalSeqCount % 1000000 == 0) {
            logger.log(Level.INFO, "Total sequences so far: " + totalSeqCount);
        }

        int bucket = getBucketForMass(precMass);
        if (bucket > NUM_BUCKETS - 1) {
            logger.log(Level.SEVERE, "Cannot add to index, unsupported precursor mass: "
                    + precMass + " for sequence: "
                    + sequence + ". Max supported mass is: " + MAX_MASS);
            return;
        }
        //logger.log(Level.INFO, "Inserting into bucket: " + bucket);
        buckets[bucket].addSequence(precMass, seqOffset, seqLength, sequence, resLeft, resRight, proteinId);

    }

    @Override
    public ResidueInfo getResidues(IndexedSequence peptideSequence, IndexedProtein protein) throws DBIndexStoreException {
        //get the true residues for this peptide
        final String proteinSequence = this.proteinCache.getProteinSequence(protein.getId());

        int seqOffset = peptideSequence.getSequenceOffset();
        if (seqOffset == IndexedSequence.OFFSET_UNKNOWN) {
            String pepSeq = peptideSequence.getSequence();
            seqOffset = proteinSequence.indexOf(pepSeq);
        }
        if (seqOffset == -1) {
            throw new RuntimeException("Could not get subsequence, unexpected error: peptide: "
                    + peptideSequence + ", protein: " + protein);
        }
        int seqLen = peptideSequence.getSequenceLen();

        return Util.getResidues(peptideSequence, seqOffset, seqLen, proteinSequence);

    }

    @Override
    public List<IndexedSequence> getSequences(float precMass, float tolerance) throws DBIndexStoreException {
        if (!inited) {
            throw new DBIndexStoreException("Indexer is not initialized");
        }

        List<IndexedSequence> ret = new ArrayList<IndexedSequence>();
        //logger.log(Level.FINE, "Starting peptide sequences query");
        //long start = System.currentTimeMillis();

        float minMass = precMass - tolerance;

        if (minMass < 0) {
            minMass = 0;
        }
        float maxMass = precMass + tolerance;


//System.out.println("==============" + precMass + " ==" + tolerance + " " + minMass + " " + maxMass);
        int[] bucketRange = getBucketsForMassRange(minMass, maxMass);
        if (bucketRange[0] > NUM_BUCKETS - 1 || bucketRange[1] > NUM_BUCKETS - 1) {
            logger.log(Level.WARNING, "Cannot query, unsupported precursor mass: "
                    + precMass + ". Max supported mass is: " + MAX_MASS);
            return ret;
        }

        for (int bucket = bucketRange[0]; bucket <= bucketRange[1]; ++bucket) {
            List<IndexedSequence> sequencesPerBucket = buckets[bucket].getSequences(precMass, tolerance);
            ret.addAll(sequencesPerBucket);
        }

        //long end = System.currentTimeMillis();
        //logger.log(Level.INFO, "Peptide sequences query took: " + (end - start) + " milisecs, results: " + ret.size());

        return ret;
    }

    @Override
    public List<IndexedSequence> getSequences(List<MassRange> ranges) throws DBIndexStoreException {
        if (ranges.size() == 1) {
            //simple case - single range, use the single range version
            MassRange range = ranges.get(0);
            return getSequences(range.getPrecMass(), range.getTolerance());
        }

        //long start = System.currentTimeMillis();

        //handle multiple ranges
        if (!inited) {
            throw new DBIndexStoreException("Indexer is not initialized");
        }

        //convert mass ranges to intverals and merge overlapping mass ranges
        ArrayList<Interval> intervals = new ArrayList<Interval>();
        for (MassRange range : ranges) {
            Interval ith = Interval.massRangeToInterval(range);
            intervals.add(ith);
        }
        List<Interval> mergedIntervals = MergeIntervals.mergeIntervals(intervals);
        

        List<IndexedSequence> ret = new ArrayList<IndexedSequence>();

        //contruct composite queries per bucket
        //store mass ranges falling in every bucket
        //using non-generic list, because arrays don't support generic objects
        //and using [][] requires extra conversion later to adapt to API
        final List bucketIntervals[] = new ArrayList[buckets.length];

        for (Interval massInterval : mergedIntervals) {
            
            float minMass = massInterval.getStart();
            float maxMass = massInterval.getEnd();

            //figure out which buckets the range falls in
            int[] bucketRange = getBucketsForMassRange(minMass, maxMass);
            if (bucketRange[0] > NUM_BUCKETS - 1 || bucketRange[1] > NUM_BUCKETS - 1) {
                logger.log(Level.SEVERE, "Cannot query, unsupported precursor mass: "
                        + minMass + "-" + maxMass + ". Max supported mass is: " + MAX_MASS);
                return ret;
            }

            //add massrange to every bucket that qualifies
            for (int bucket = bucketRange[0]; bucket <= bucketRange[1]; ++bucket) {
                List<Interval> bucketR = bucketIntervals[bucket];
                if (bucketR == null) {
                    //lazy-initialize the bucket first time
                    bucketR = new ArrayList<Interval>();
                    bucketIntervals[bucket] = bucketR;
                }

                //append mass range
                if (!bucketR.contains(massInterval)) {
                    bucketR.add(massInterval);
                }
            }
        }


        //query the buckets that have mass ranges assigned
        for (int bucket = 0; bucket < buckets.length; ++bucket) {
            List<Interval> queryBucketIntervals = bucketIntervals[bucket];
            if (queryBucketIntervals == null) {
                continue; //skip, no intervals for this bucket
            }

            List<IndexedSequence> sequencesPerBucket = buckets[bucket].getSequencesIntervals(queryBucketIntervals);
            ret.addAll(sequencesPerBucket);
        }


        //long end = System.currentTimeMillis();
        //logger.log(Level.INFO, "Peptide sequences query took: " + (end - start) + " milisecs, results: " + ret.size());

        return ret;

    }

    @Override
    public boolean supportsProteinCache() {
        return true;
    }

    @Override
    public void setProteinCache(ProteinCache proteinCache) {
        this.proteinCache = proteinCache;
        for (int i = 0; i < NUM_BUCKETS; ++i) {
            buckets[i].setProteinCache(proteinCache);
        }
    }

    @Override
    public long addProteinDef(long num, String definition, String proteinSequence) throws DBIndexStoreException {
        //using protein cache entirely
        //nothing here
        return num;
    }

    @Override
    public List<IndexedProtein> getProteins(IndexedSequence sequence) throws DBIndexStoreException {
        List<IndexedProtein> ret = new ArrayList<IndexedProtein>();

        for (Integer protId : sequence.getProteinIds()) {
            IndexedProtein ip = new IndexedProtein(proteinCache.getProteinDef(protId), protId);
            ret.add(ip);
        }

        return ret;

    }

    /**
     * test driver
     *
     * @param args
     */
    public static void main(String[] args) {
        DBIndexStore store = new DBIndexStoreSQLiteMult(SearchParams.getInstance(), false);

        String protDef1 = "4R79.2 CE19650 WBGene00007067 Ras family status:Partially_confirmed TR:Q9XXA4 protein_id:CAA20282.1";
        String protDef2 = "Reverse_4R79.2  CE19650 WBGene00007067 Ras family status:Partially_confirmed TR:Q9XXA4 protein_id:CAA20282.1";

        try {
            logger.log(Level.INFO, "Testing init");
            File idx = new File("./test.fasta.idx");
            if (idx.exists() && idx.isDirectory()) {
                for (File f : idx.listFiles()) {
                    f.delete();
                }
                idx.delete();
            }

            ProteinCache protCache = ProteinCache.getInstance();
            store.init("./test.fasta");
            store.setProteinCache(protCache);

            //test a transaction
            logger.log(Level.INFO, "Testing adds 1");
            store.startAddSeq();

            long protId = store.addProteinDef(0, protDef1, "ABCDEFGHIJKL");
            protCache.addProtein("4R79.2", "ABCDEFGHIJKL");

            store.addSequence(1f, 0, 1, "A", null, null, protId);
            store.addSequence(2f, 0, 2, "AB", null, null, protId);
            store.addSequence(3f, 0, 3, "ABC", null, null, protId);
            store.addSequence(4f, 0, 4, "ABCD", null, null, protId);
            store.addSequence(6000.42323f, 0, 5, "ABCDE", null, null, protId);
            store.addSequence(6999.42323f, 0, 6, "ABCDEZ", null, null, protId);
            store.addSequence(3, 6, 3, "GHI", null, null, protId);
            //store.addSequence(3, 6, 3, "GHI", null, null, protId);
            //store.stopAddSeq();
        } catch (DBIndexStoreException ex) {
            logger.log(Level.SEVERE, null, ex);
        }

        try {
            //test more adds with the same statement in a new transaction
            logger.log(Level.INFO, "Testing adds 2");
            //store.startAddSeq();

            long protId = store.addProteinDef(1, protDef2, "GHIJKLMNOPR");
            ProteinCache protCache = ProteinCache.getInstance();
            protCache.addProtein("Reverse_4R79.2", "GHIJKLMN");

            store.addSequence(3, 1, 3, "HIJ", null, null, protId);
            store.addSequence(5, 2, 5, "IJKLM", null, null, protId);
            //store.addSequence(7, 2, 7, "IJKLMNO", null, null, protId);
            //store.addSequence(7.1f, 2, 8, "IJKLMNOO", null, null, protId);
            store.addSequence(3, 0, 3, "GHI", null, null, protId); //test dup
            store.addSequence(3, 0, 3, "GHI", null, null, protId); //test dup
            store.stopAddSeq();
        } catch (DBIndexStoreException ex) {
            logger.log(Level.SEVERE, null, ex);
        }

        try {
            //test gets
            logger.log(Level.INFO, "Testing gets 1");
            List<IndexedSequence> res1 = store.getSequences(10, 8.9f);
            for (IndexedSequence seq : res1) {
                System.out.println(seq);
                List<IndexedProtein> proteins = store.getProteins(seq);
                for (IndexedProtein protein : proteins) {
                    ResidueInfo res = store.getResidues(seq, protein);
                    System.out.println(protein);
                    System.out.println(res);
                }
            }
        } catch (DBIndexStoreException ex) {
            logger.log(Level.SEVERE, null, ex);
        }


        try {
            //test gets
            System.out.println("Testing range gets 1");
            MassRange r1 = new MassRange(2, 1f);
            MassRange r2 = new MassRange(6, 1f);
            MassRange r3 = new MassRange(6, 1f); //test redundant
            MassRange r4 = new MassRange(6, 1.2f); //test overlapping
            List<MassRange> ranges = new ArrayList<MassRange>();
            ranges.add(r2);
            ranges.add(r1);
            ranges.add(r3);
            ranges.add(r4);
            List<IndexedSequence> res1 = store.getSequences(ranges);
            for (IndexedSequence seq : res1) {
                System.out.println(seq);
                List<IndexedProtein> proteins = store.getProteins(seq);
                for (IndexedProtein protein : proteins) {
                    ResidueInfo res = store.getResidues(seq, protein);
                    System.out.println(protein);
                    System.out.println(res);
                }
            }
        } catch (DBIndexStoreException ex) {
            logger.log(Level.SEVERE, null, ex);
        }


        if (true) {
            return;
        }


        try {
            store = new DBIndexStoreSQLiteMult(SearchParams.getInstance(), false);
            store.init("EBI-IPI_Human_IPI_Human_3_85_06-29-2011_reversed.fasta");
            for (int i = 0; i < 1000; ++i) {
                logger.log(Level.INFO, "Testing getting sequence " + i);
                List<IndexedSequence> res = store.getSequences(i, 3f);
                if (res.size() > 0) {
                    IndexedSequence seq = res.get(0);
                    logger.log(Level.INFO, "Got " + res.size() + " sequences, first sequence: " + seq);
                    List<IndexedProtein> proteins = store.getProteins(seq);
                    logger.log(Level.INFO, "Got " + res.size() + " sequences, first sequence's proteins: ");
                    for (IndexedProtein protein : proteins) {
                        System.out.println(protein);
                    }
                }

            }

        } catch (DBIndexStoreException ex) {
            logger.log(Level.SEVERE, null, ex);
        }

        try {
            logger.log(Level.INFO, "Testing getting multiple sequence proteins ");
            store = new DBIndexStoreSQLiteMult(SearchParams.getInstance(), false);
            store.init("test_02.fasta");

            List<IndexedSequence> res = store.getSequences(1000, 3000);
            for (IndexedSequence seq : res) {
                System.out.println(seq);
                List<IndexedProtein> prots = store.getProteins(seq);
                for (IndexedProtein prot : prots) {
                    System.out.println("\t" + prot);
                }
            }

        } catch (DBIndexStoreException ex) {
            logger.log(Level.SEVERE, null, ex);
        }

    }

    @Override
    public void lastBuffertoDatabase() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    /**
     * Simple table implementation with: blazmass_sequences integer mass (PRI
     * KEY, ASC and DESC index for quick range lookup) byte data : encoded
     * peptide info (prot id, offset, length)
     *
     * Meant to be used by DBIndexStoreSQLiteMult only, internal private class
     *
     * Caches sequences in intermediate structure and commits to db every 200k
     * sequences.
     *
     */
    private static class DBIndexStoreSQLiteByte extends DBIndexStoreSQLiteAbstract {

        private PreparedStatement addSeqStatement;
        protected PreparedStatement updateSeqStatement;
        private PreparedStatement updateAppendSeqStatement;
        protected PreparedStatement getSeqStatement;
        protected PreparedStatement[] getSeqMassRanges2_Statements;
        private PreparedStatement getSeqDataStatement;
        private PreparedStatement getSeqExistsStatement;
        private final String blazmass_sequences_pri_key = "precursor_mass_key";
        private final DynByteBuffer data[] = new DynByteBuffer[MAX_MASS]; //we really only need MAX_MASS / NUM_BUCKETS, but need to shift
        private int bucketId;
        //single sequence commit mode coutners
        private int curCommit = 0;
        private static final int SEQ_UPDATE_INTERVAL = 20000;
        //seq count since last commit of entire cache
        private int seqCount = 0;
        //total sequence count indexed
        private int totalSeqCount = 0;
        //full cache commit interval
        //NOTE: optimized for 7GB heap space, adjsut this setting accordingly
        private static final int FULL_CACHE_COMMIT_INTERVAL = (100 * 1000 * 1000) / NUM_BUCKETS;

        ;

        DBIndexStoreSQLiteByte(int bucketId) {
            super();
            this.bucketId = bucketId;
        }

        DBIndexStoreSQLiteByte(boolean inMemory, int bucketId) {
            super(SearchParams.getInstance(), inMemory);
            this.bucketId = bucketId;
        }

        @Override
        protected void initStatementsCommon() throws DBIndexStoreException {
        }

        @Override
        protected void initStatements() throws DBIndexStoreException {
            try {

                addSeqStatement = con.prepareStatement(
                        "INSERT INTO " + getIndexTableName() + " (precursor_mass_key, data) "
                        + "VALUES (?, ?);");

                updateSeqStatement = con.prepareStatement(
                        "UPDATE " + getIndexTableName() + " SET data = ? "
                        + "WHERE precursor_mass_key = ?;");

                updateAppendSeqStatement = con.prepareStatement(
                        "UPDATE " + getIndexTableName() + " SET data = (data || ?) "
                        + "WHERE precursor_mass_key = ?;");


                getSeqStatement = con.prepareStatement(
                        "SELECT precursor_mass_key, data "
                        + "FROM " + getIndexTableName() + " "
                        + "WHERE precursor_mass_key BETWEEN ? AND ?;");


                getSeqMassRanges2_Statements = new PreparedStatement[MAX_MASS_RANGES];

                StringBuilder stSb = new StringBuilder();
                stSb.append("SELECT DISTINCT precursor_mass_key, data ");
                stSb.append("FROM ").append(getIndexTableName()).append(" WHERE ");
                for (int st = 0; st < MAX_MASS_RANGES; ++st) {
                    if (st > 0) {
                        stSb.append(" OR ");
                    }
                    stSb.append("precursor_mass_key BETWEEN ? AND ? ");
                    getSeqMassRanges2_Statements[st] = con.prepareStatement(stSb.toString());
                }


                getSeqDataStatement = con.prepareStatement(
                        "SELECT data "
                        + "FROM " + getIndexTableName() + " "
                        + "WHERE precursor_mass_key = ?;");

                getSeqExistsStatement = con.prepareStatement(
                        "SELECT precursor_mass_key "
                        + "FROM " + getIndexTableName() + " "
                        + "WHERE precursor_mass_key = ? LIMIT 1;");


            } catch (SQLException e) {
                logger.log(Level.SEVERE, "Error initializing statements in db, path: " + dbPath, e);
                throw new DBIndexStoreException("Error initializing statements in db, path: " + dbPath, e);
            }
        }

        @Override
        protected String getSequencesTablePriKey() {
            return blazmass_sequences_pri_key;
        }

        /**
         * Get sequence data for long representation of precursor mass or null
         * if not present
         *
         * @param massKey
         * @return encoded sequence data string for multiple peptide, or null
         * @throws SQLException
         */
        private byte[] getSequenceData(int massKey) throws SQLException {
            byte[] ret = null;
            getSeqDataStatement.setInt(1, massKey);
            final ResultSet rs = getSeqDataStatement.executeQuery();
            if (rs.next()) {
                ret = rs.getBytes(1);
            }
            rs.close();

            return ret;
        }

        /**
         * check if sequence exists by mass
         *
         * @param massKey
         * @return
         * @throws SQLException
         */
        private boolean getSequenceExists(int massKey) throws SQLException {
            boolean ret = false;
            getSeqExistsStatement.setInt(1, massKey);
            final ResultSet rs = getSeqExistsStatement.executeQuery();
            if (rs.next()) {
                ret = true;
            }
            rs.close();

            return ret;
        }

        /**
         * Add the new sequence to in memory cached pre-commit sequence hash
         *
         * @param precMass
         * @param seqOffset
         * @param seqLength
         * @param proteinId
         */
        private void updateCachedData(float precMass, int seqOffset, int seqLength, int proteinId) {
            int rowId = (int) (precMass * MASS_GROUP_FACTOR);

            DynByteBuffer byteBuffer = data[rowId];
            if (byteBuffer == null) {
                byteBuffer = new DynByteBuffer();
                data[rowId] = byteBuffer;
            }



            //store in Little-Endian order
            byte[] seqMassB = DynByteBuffer.toByteArray(precMass);
            byteBuffer.add(seqMassB);

            byte[] seqOffsetB = DynByteBuffer.toByteArray(seqOffset);
            byteBuffer.add(seqOffsetB);

            byte[] seqLengthB = DynByteBuffer.toByteArray(seqLength);
            byteBuffer.add(seqLengthB);

            byte[] proteinIdB = DynByteBuffer.toByteArray(proteinId);
            byteBuffer.add(proteinIdB);

            //commit this mass after COMMIT_SEQUENCES sequences
            if (true) {
                if (byteBuffer.getSize() > COMMIT_SEQUENCES * BYTE_PER_SEQUENCE) {
                    try {
                        insertSequence(rowId, byteBuffer);
                        //clear since we wrote it to db
                        byteBuffer.clear();
                    } catch (SQLException e) {
                        logger.log(Level.SEVERE, "Error commiting mass buffer");
                    }
                }
            }
        }

        /**
         * updates db rows with the sequence, optionally commits it
         *
         * @param rowId
         * @param buf
         * @throws SQLException
         */
        private void insertSequence(int rowId, DynByteBuffer buf) throws SQLException {
            //logger.log(Level.INFO, bucketId + ": Commiting cached sequence data for rowId: " + rowId);
            if (getSequenceExists(rowId)) {
                //merge, bytes in db with the bytes in cache
                //build up updated string

                updateAppendSeqStatement.setBytes(1, buf.getData());
                updateAppendSeqStatement.setInt(2, rowId);
                updateAppendSeqStatement.executeUpdate();
            } else {
                //write a brand new row for the mass

                //insert 
                addSeqStatement.setInt(1, rowId);
                addSeqStatement.setBytes(2, buf.getData());
                addSeqStatement.executeUpdate();
            }

            if (curCommit++ == SEQ_UPDATE_INTERVAL) {
                curCommit = 0;
                // con.commit(); //do not db commit every sequence for now, let the big cache commit handle it every FULL_CACHE_COMMIT_INTERVAL
            }

        }

        @Override
        protected final void commitCachedData() throws SQLException {
            //called at the end of indexing to commit all buffers
            //logger.log(Level.INFO, "Bucket " + bucketId + ": Commiting cached data");
            for (int massKey = 0; massKey < MAX_MASS; ++massKey) {
                //see if already in index

                DynByteBuffer cached = data[massKey];
                if (cached == null || cached.getSize() == 0) {
                    continue;
                }

                if (getSequenceExists(massKey)) {
                    //merge, bytes in db with the bytes in cache
                    //build up updated string


                    updateAppendSeqStatement.setBytes(1, cached.getData());
                    updateAppendSeqStatement.setInt(2, massKey);
                    updateAppendSeqStatement.executeUpdate();
                } else {
                    //write a brand new row for the mass

                    //insert 
                    addSeqStatement.setInt(1, massKey);
                    addSeqStatement.setBytes(2, cached.getData());
                    addSeqStatement.executeUpdate();
                }

            }

            //clear cached data
            Arrays.fill(data, null);


            //logger.log(Level.INFO, "Commit start");
            con.commit();
            //logger.log(Level.INFO, "Commit end");


            // logger.log(Level.INFO, bucketId + ": Done commiting cached data");



        }

        @Override
        public FilterResult filterSequence(float precMass, String sequence) {
            //no filtering, we add every sequence
            return FilterResult.INCLUDE;
        }

        @Override
        public void addSequence(float precMass, int seqOffset, int seqLength, String sequence, String resLeft, String resRight, long proteinId) throws DBIndexStoreException {
            if (!inited) {
                throw new DBIndexStoreException("Indexer is not initialized");
            }

            totalSeqCount++;

            updateCachedData(precMass, seqOffset, seqLength, (int) proteinId);

            if (true) {
                //despite commiting each sequence, commit and clear all after 100mil
                //to free up memory, optimize buffers, and ensure a commit
                if (seqCount++ == FULL_CACHE_COMMIT_INTERVAL) {
                    //logger.log(Level.INFO, "Commiting all");
                    seqCount = 0;
                    try {
                        commitCachedData();
                    } catch (SQLException ex) {
                        logger.log(Level.SEVERE, "Error commiting cached data", ex);
                    }
                }
            }

        }

        @Override
        public ResidueInfo getResidues(IndexedSequence peptideSequence, IndexedProtein protein) throws DBIndexStoreException {
            //implemented by outter class
            return null;

        }

        @Override
        public List<IndexedSequence> getSequences(float precMass, float tolerance) throws DBIndexStoreException {

            if (!inited) {
                throw new DBIndexStoreException("Indexer is not initialized");
            }

            //logger.log(Level.FINE, "Starting peptide sequences query");
            //long start = System.currentTimeMillis();

            float toleranceFactor = tolerance; //precMass * tolerance;
            float minMassF = precMass - toleranceFactor;
            if (minMassF < 0f) {
                minMassF = 0;
            }
            float maxMassF = precMass + toleranceFactor;

            int precMassInt = (int) (precMass * MASS_GROUP_FACTOR);
            //long toleranceInt = (long) (tolerance * MASS_STORE_MULT);
            int toleranceInt = (int) (toleranceFactor * MASS_GROUP_FACTOR); //Robin fixed the ppm calculation

            int minMass = precMassInt - toleranceInt;
            //System.out.println(precMassInt + " " + toleranceInt);

            if (minMass < 0) {
                minMass = 0;
            }
            int maxMass = precMassInt + toleranceInt;

            List<IndexedSequence> ret = new ArrayList<IndexedSequence>();

            ResultSet rs = null;
            try {
                getSeqStatement.setInt(1, minMass - PPM_PER_ENTRY);
                getSeqStatement.setInt(2, maxMass + PPM_PER_ENTRY);
                rs = getSeqStatement.executeQuery();

                while (rs.next()) {
                    //final int massKey = rs.getInt(1);
                    final byte[] seqData = rs.getBytes(2);

                    parseAddPeptideInfo(seqData, ret, minMassF, maxMassF);

                }


            } catch (SQLException e) {
                String msg = "Error getting peptides ";
                logger.log(Level.SEVERE, msg, e);
                throw new DBIndexStoreException(msg, e);
            } finally {
                if (rs != null) {
                    try {
                        rs.close();
                    } catch (SQLException ex) {
                        logger.log(Level.SEVERE, null, ex);
                    }
                }
            }

            //long end = System.currentTimeMillis();
            //logger.log(Level.INFO, "Peptide sequences query took: " + (end - start) + " milisecs, results: " + ret.size());

            return ret;
        }

        @Override
        public List<IndexedSequence> getSequences(List<MassRange> ranges) throws DBIndexStoreException {
            if (ranges.size() == 1) {
                MassRange range = ranges.get(0);
                return getSequences(range.getPrecMass(), range.getTolerance());
            }

            throw new DBIndexStoreException("Multiple mass ranges not implemented yet!");
        }

        /**
         * Parse the encoded sequences and return them in toInsert
         *
         * @param data
         * @param toInsert
         * @param minMass
         * @param maxMass
         */
        private void parseAddPeptideInfo(byte[] data, List<IndexedSequence> toInsert, float minMass, float maxMass) {

            //to collapse multiple sequences into single one, with multproteins
            Map<String, List<IndexedSeqInternal>> temp = new HashMap<String, List<IndexedSeqInternal>>();

            int dataLength = data.length;


            //if (dataLength % 4 != 0) {
            //  throw new RuntimeException("Unexpected number of peptide items: " + dataLength);
            //}

            for (int i = 0; i < dataLength; i += 16) {
                final int first = i + 4;
                final int second = i + 8;
                final int third = i + 12;
                final int fourth = i + 16;

                byte[] slice = Arrays.copyOfRange(data, i, first);
                float seqMass = DynByteBuffer.toFloat(slice);

                if (seqMass < minMass || seqMass > maxMass) {
                    //skip the sequence, it qualified the bucket, but not the actual mass 
                    continue;
                }

                slice = Arrays.copyOfRange(data, first, second);
                int offset = DynByteBuffer.toInt(slice);
                slice = Arrays.copyOfRange(data, second, third);
                int length = DynByteBuffer.toInt(slice);
                slice = Arrays.copyOfRange(data, third, fourth);
                int proteinId = DynByteBuffer.toInt(slice);

                String peptideSequence = null;
                //System.out.print("prot id: " + proteinId + " offset: " + offset + " len: " + length);

                peptideSequence = proteinCache.getPeptideSequence(proteinId, offset, length);
                //System.out.println("pep: " + peptideSequence);

                //we are cheating and supplying protein id instead of peptide id
                //to set it temporarily, before we merge them into a list
                final IndexedSeqInternal tempSequence = new IndexedSeqInternal(seqMass, offset, length, proteinId, peptideSequence);


                List<IndexedSeqInternal> sequences = temp.get(peptideSequence);
                if (sequences == null) {
                    sequences = new ArrayList<IndexedSeqInternal>();
                    temp.put(peptideSequence, sequences);
                }
                sequences.add(tempSequence);

            }

            //group the same peptides from many proteins into single peptide with protein id list
            for (String pepSeqKey : temp.keySet()) {

                List<IndexedSeqInternal> sequences = temp.get(pepSeqKey);

                Set<Integer> proteinIds = new HashSet<Integer>();
                for (IndexedSeqInternal tempSeq : sequences) {
                    proteinIds.add((int) tempSeq.proteinId);
                }

                IndexedSeqInternal firstSeq = sequences.get(0);
                IndexedSequence mergedSequence = new IndexedSequence(0, firstSeq.mass, pepSeqKey, "", "");
                mergedSequence.setProteinIds(new ArrayList<Integer>(proteinIds));
                //set residues
                final String protSequence = proteinCache.getProteinSequence(firstSeq.proteinId);
                ResidueInfo residues = Util.getResidues(null, firstSeq.offset, firstSeq.length, protSequence);
                mergedSequence.setResidues(residues);
                toInsert.add(mergedSequence);
            }


        }

        @Override
        public boolean supportsProteinCache() {
            return true;
        }

        @Override
        public void setProteinCache(ProteinCache proteinCache) {
            this.proteinCache = proteinCache;
        }

        @Override
        public long addProteinDef(long num, String definition, String proteinSequence) throws DBIndexStoreException {
            //implemented by outter class
            return num;

        }

        @Override
        public List<IndexedProtein> getProteins(IndexedSequence sequence) throws DBIndexStoreException {
            //implemented by outter class
            return null;

        }

        @Override
        protected void createTempIndex() {
            //not needed
        }

        @Override
        protected void deleteTempIndex() {
            //not needed
        }

        protected String getIndexTableName() {
            return "blazmass_sequences";
        }

        @Override
        protected void createTables() throws DBIndexStoreException {
            try {

                executeStatement("CREATE TABLE IF NOT EXISTS " + getIndexTableName() + " "
                        + "(precursor_mass_key INTEGER PRIMARY KEY, "
                        + "data BINARY"
                        + ");");

            } catch (SQLException e) {
                logger.log(Level.SEVERE, "Error creating tables, ", e);
                throw new DBIndexStoreException("Error creating tables. ", e);
            }
        }

        @Override
        protected void createIndex() throws DBIndexStoreException {
            try {
                //for fast lookup by mass range 
                // logger.log(Level.INFO, "Creating precursor_mass_key_index_asc ");
                // executeStatement("CREATE INDEX IF NOT EXISTS precursor_mass_key_index_asc ON blazmass_sequences (precursor_mass_key ASC);");
                //dsc index speeds up the range queries another 10x more
                //logger.log(Level.INFO, "Creating precursor_mass_key_index_dsc ");
                executeStatement("CREATE INDEX IF NOT EXISTS precursor_mass_key_index_dsc ON " + getIndexTableName() + " (precursor_mass_key DESC);");
            } catch (SQLException e) {
                logger.log(Level.SEVERE, "Error creating index, ", e);
                throw new DBIndexStoreException("Error creating index, ", e);
            }
        }

        @Override
        protected void closeConnection() {
            try {

                if (addSeqStatement != null) {
                    addSeqStatement.close();
                    addSeqStatement = null;
                }

                if (updateSeqStatement != null) {
                    updateSeqStatement.close();
                    updateSeqStatement = null;
                }

                if (updateAppendSeqStatement != null) {
                    updateAppendSeqStatement.close();
                    updateAppendSeqStatement = null;
                }

                if (getSeqStatement != null) {
                    getSeqStatement.close();
                    getSeqStatement = null;
                }

                for (int i = 0; i < MAX_MASS_RANGES; ++i) {
                    if (getSeqMassRanges2_Statements[i] != null) {
                        getSeqMassRanges2_Statements[i].close();
                        getSeqMassRanges2_Statements[i] = null;
                    }
                }

                if (getSeqDataStatement != null) {
                    getSeqDataStatement.close();
                    getSeqDataStatement = null;
                }

                if (getSeqExistsStatement != null) {
                    getSeqExistsStatement.close();
                    getSeqExistsStatement = null;
                }

                if (con != null) {
                    con.close();
                    con = null;
                }
            } catch (SQLException ex) {
                logger.log(Level.SEVERE, "Error closing SQLite connection", ex);
            }
        }

        @Override
        public long getNumberSequences() throws DBIndexStoreException {
            if (!inited) {
                throw new DBIndexStoreException("Indexer is not initialized");
            }

            ResultSet rs = null;
            try {
                //we don't really now, return number of rows for now
                rs = this.executeQuery("SELECT COUNT(*) FROM " + getIndexTableName());
                return rs.getLong(1);
            } catch (SQLException ex) {
                logger.log(Level.SEVERE, "Error executing query to get number of sequences.", ex);
                throw new DBIndexStoreException("Error executing query to get number of sequences.", ex);
            } finally {
                if (rs != null) {
                    try {
                        rs.close();
                    } catch (SQLException ex) {
                        logger.log(Level.SEVERE, null, ex);
                    }
                }
            }


        }

        @Override
        public void lastBuffertoDatabase() {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
    }

    /**
     * a variation that performs merge at index time instead of merge at search
     * time
     */
    private static final class DBIndexStoreSQLiteByteIndexMerge extends DBIndexStoreSQLiteByte {

        private static final Logger logger = Logger.getLogger(DBIndexStoreSQLiteByteIndexMerge.class.getName());
        private static final int SEQ_SEPARATOR_INT = Integer.MAX_VALUE;
        private static final byte[] SEQ_SEPARATOR = DynByteBuffer.toByteArray(SEQ_SEPARATOR_INT); //read protein ids until hitting the separator
        private PreparedStatement getMassDataStatement;

        DBIndexStoreSQLiteByteIndexMerge(int bucketId) {
            super(bucketId);
        }

        DBIndexStoreSQLiteByteIndexMerge(boolean inMemory, int bucketId) {
            super(inMemory, bucketId);
        }

        @Override
        protected void initStatements() throws DBIndexStoreException {
            try {
                super.initStatements();

                getMassDataStatement = con.prepareStatement(
                        "SELECT precursor_mass_key, data "
                        + "FROM " + getIndexTableName());
            } catch (SQLException ex) {
                throw new DBIndexStoreException("Error initializing get mass data st", ex);
            }
        }

        @Override
        protected void createIndex() throws DBIndexStoreException {
            //using the createIndex hook to merge peptides at end of the stop()

            logger.log(Level.INFO, "Merging sequences and compacting the index db: " + this.dbPath);
            mergePeptides();

            try {
                //vacuum, should trim 30-50% size
                //also consider TEMP tables
                executeStatement("VACUUM;");
            } catch (SQLException ex) {
                logger.log(Level.SEVERE, "Error compacting the index after the merge", ex);
            }


            //now can create the index
            super.createIndex();

            logger.log(Level.INFO, "Done merging sequences and compacting the index db: " + this.dbPath);
        }

        private void mergePeptides() throws DBIndexStoreException {
            // logger.log(Level.INFO, "Merging sequence index start");
            ResultSet massKeysRs = null;
            try {
                con.setAutoCommit(false);

                massKeysRs = getMassDataStatement.executeQuery();

                //go over each row and merge peptides
                while (massKeysRs.next()) {
                    final int massKey = massKeysRs.getInt(1);
                    final byte[] seqData = massKeysRs.getBytes(2);

                    final byte[] seqDataMerged = getMergedData(seqData);
                    updateSeqStatement.setBytes(1, seqDataMerged);
                    updateSeqStatement.setInt(2, massKey);
                    updateSeqStatement.execute();

                }

                con.commit();

            } catch (SQLException ex) {
                logger.log(Level.SEVERE, "Error merging sequences", ex);
            } finally {
                try {
                    if (massKeysRs != null) {
                        massKeysRs.close();
                    }
                    con.setAutoCommit(true);
                } catch (SQLException ex) {
                    logger.log(Level.SEVERE, "Error restoring autocommit", ex);
                }
            }
            //logger.log(Level.INFO, "Merging sequence index end");
        }

        @Override
        protected void closeConnection() {

            if (getMassDataStatement != null) {
                try {
                    getMassDataStatement.close();
                    getMassDataStatement = null;
                } catch (SQLException ex) {
                    logger.log(Level.SEVERE, "Error closing get mass data st.", ex);
                }
            }

            super.closeConnection();

        }

        @Override
        public List<IndexedSequence> getSequences(float precMass, float tolerance) throws DBIndexStoreException {
            if (!inited) {
                throw new DBIndexStoreException("Indexer is not initialized");
            }

            //logger.log(Level.FINE, "Starting peptide sequences query");
            //long start = System.currentTimeMillis();

            float toleranceFactor = tolerance; //precMass * tolerance;
            float minMassF = precMass - toleranceFactor;
            if (minMassF < 0f) {
                minMassF = 0f;
            }
            float maxMassF = precMass + toleranceFactor;

            int precMassInt = (int) (precMass * MASS_GROUP_FACTOR);
            //long toleranceInt = (long) (tolerance * MASS_STORE_MULT);
            int toleranceInt = (int) (toleranceFactor * MASS_GROUP_FACTOR); //Robin fixed the ppm calculation

            int minMass = precMassInt - toleranceInt;
            //System.out.println(precMassInt + " " + toleranceInt);

            if (minMass < 0) {
                minMass = 0;
            }
            int maxMass = precMassInt + toleranceInt;

            List<IndexedSequence> ret = new ArrayList<IndexedSequence>();

            ResultSet rs = null;
            try {
                getSeqStatement.setInt(1, minMass - PPM_PER_ENTRY);
                getSeqStatement.setInt(2, maxMass + PPM_PER_ENTRY);
                rs = getSeqStatement.executeQuery();

                while (rs.next()) {
                    //final int massKey = rs.getInt(1);
                    final byte[] seqData = rs.getBytes(2);

                    parseAddPeptideInfo(seqData, ret, minMassF, maxMassF);

                }


            } catch (SQLException e) {
                String msg = "Error getting peptides ";
                logger.log(Level.SEVERE, msg, e);
                throw new DBIndexStoreException(msg, e);
            } finally {
                if (rs != null) {
                    try {
                        rs.close();
                    } catch (SQLException ex) {
                        logger.log(Level.SEVERE, "Error closing result set after getting sequences", ex);
                    }
                }
            }

            //long end = System.currentTimeMillis();
            //logger.log(Level.INFO, "Peptide sequences query took: " + (end - start) + " milisecs, results: " + ret.size());

            return ret;
        }

        @Override
        public List<IndexedSequence> getSequences(List<MassRange> ranges) throws DBIndexStoreException {

            int numRanges = ranges.size();

            if (numRanges == 1) {
                //use the standard method for the single range, as it is slightly more optimized
                MassRange range = ranges.get(0);
                return getSequences(range.getPrecMass(), range.getTolerance());
            }

            ArrayList<Interval> intervals = new ArrayList<Interval>();

            for (MassRange range : ranges) {
                Interval ith = Interval.massRangeToInterval(range);
                intervals.add(ith);
            }

            return getSequencesIntervals(intervals);

        }

        public List<IndexedSequence> getSequencesIntervals(List<Interval> ranges) throws DBIndexStoreException {

            int numRanges = ranges.size();

            List<IndexedSequence> ret = new ArrayList<IndexedSequence>();

            PreparedStatement rangesStatement = null;
            boolean tempStatement = false;
            if (numRanges <= MAX_MASS_RANGES) {
                rangesStatement = getSeqMassRanges2_Statements[numRanges - 1];
            } else {
                //fallback to build non-optimized statement on the fly

                StringBuilder stSb = new StringBuilder();
                stSb.append("SELECT DISTINCT precursor_mass_key, data ");
                stSb.append("FROM ").append(getIndexTableName()).append(" WHERE ");

                for (int st = 0; st < ranges.size(); ++st) {
                    Interval range = ranges.get(st);
                    if (st > 0) {
                        stSb.append(" OR ");
                    }

                    float minMassF = range.getStart();
                    float maxMassF = range.getEnd();
                    stSb.append("precursor_mass_key BETWEEN ").append(minMassF).append(" AND ").append(maxMassF);
                }
                try {
                    rangesStatement = con.prepareStatement(stSb.toString());
                } catch (SQLException ex) {
                    logger.log(Level.SEVERE, "Error preparing a temp statement: " + stSb.toString(), ex);
                }
                tempStatement = true;

            }

            //setup the query
            //supporting 2, 3 or 4+ ranges per query
            int curRange = 0;
            //record actual float min and max masses to reuse when filtering out sequences within a row
            float[] minMassesF = new float[numRanges];
            float[] maxMassesF = new float[numRanges];

            for (Interval range : ranges) {

                //float precMass = range.getPrecMass();
                //float tolerance = range.getTolerance();

                //float toleranceFactor = tolerance; //precMass * tolerance;
                float minMassF = range.getStart();
                float maxMassF = range.getEnd();

                //int precMassInt = (int) (precMass * MASS_GROUP_FACTOR);
                //long toleranceInt = (long) (tolerance * MASS_STORE_MULT);
                //int toleranceInt = (int) (toleranceFactor * MASS_GROUP_FACTOR); //Robin fixed the ppm calculation

                //figure out the rows to query
                int minMass = (int) minMassF;
                //System.out.println(precMassInt + " " + toleranceInt);
                int maxMass = (int) (maxMassF + PPM_PER_ENTRY);

                int startCol = 2 * curRange + 1;
                try {
                    if (tempStatement == false) {
                        //use optimized prep statements
                        rangesStatement.setInt(startCol, minMass); // - PPM_PER_ENTRY); //TODO do we need to look in another +-1 row ? I don't think so
                        rangesStatement.setInt(startCol + 1, maxMass); // + PPM_PER_ENTRY); //since BETWEEN query is inclusive
                    }
                } catch (SQLException ex) {
                    logger.log(Level.SEVERE, "Error setting up range query for range: " + curRange, ex);
                    throw new DBIndexStoreException("Error setting up range query for range: " + curRange, ex);
                }

                minMassesF[curRange] = minMassF;
                maxMassesF[curRange] = maxMassF;

                ++curRange;
            }

            //run the query
            ResultSet rs = null;

            try {
                rs = rangesStatement.executeQuery();
                while (rs.next()) {
                    //final int massKey = rs.getInt(1);
                    final byte[] seqData = rs.getBytes(2);
                    parseAddPeptideInfo(seqData, ret, minMassesF, maxMassesF);
                }
            } catch (SQLException e) {
                String msg = "Error getting peptides ";
                logger.log(Level.SEVERE, msg, e);
                throw new DBIndexStoreException(msg, e);
            } finally {
                if (rs != null) {
                    try {
                        rs.close();
                    } catch (SQLException ex) {
                        logger.log(Level.SEVERE, "Error closing result set after getting sequences", ex);
                    }
                }
                if (rangesStatement != null && tempStatement == true) {
                    try {
                        //destroy if temp statement
                        rangesStatement.close();
                    } catch (SQLException ex) {
                        logger.log(Level.SEVERE, "Error closing temp statement after getting sequences", ex);
                    }
                }
            }

            return ret;
        }

        /**
         * Parse the encoded sequences and return them in toInsert
         *
         * @param data
         * @param toInsert
         * @param minMass used to eliminate not matching masses in the 1ppm row
         * @param maxMass used to eliminate not matching masses in the 1ppm row
         */
        private void parseAddPeptideInfo(byte[] data, List<IndexedSequence> toInsert, float minMass, float maxMass) {

            int dataLength = data.length;

            //if (dataLength % 4 != 0) {
              //  throw new RuntimeException("Unexpected number of peptide items: " + dataLength);
            //}

            for (int i = 0; i < dataLength;) {
                final int zeroeth = i;
                i += 4;
                final int first = i;
                i += 4;
                final int second = i;
                i += 4;
                final int third = i;
                i += 4;
                final int fourth = i;

                byte[] slice = Arrays.copyOfRange(data, zeroeth, first);
                float seqMass = DynByteBuffer.toFloat(slice);

                //since it's sorted
                //if current mass great than maxmass, break - optimize
                if (seqMass > maxMass) {
                    break;
                }

                if (seqMass < minMass) {
                    //skip the sequence, it qualified the row, but not the actual peptide by precise mass 

                    //skip to next sequence
                    while (true) {
                        slice = Arrays.copyOfRange(data, i, i + 4);
                        i += 4;
                        int proteinId = DynByteBuffer.toInt(slice);
                        if (proteinId == SEQ_SEPARATOR_INT) {
                            break;
                        }
                    }
                    //go to next seq
                    continue;
                }


                //create sequence object to return

                slice = Arrays.copyOfRange(data, first, second);
                int offset = DynByteBuffer.toInt(slice);
                slice = Arrays.copyOfRange(data, second, third);
                int length = DynByteBuffer.toInt(slice);

                slice = Arrays.copyOfRange(data, third, fourth);
                int proteinId = DynByteBuffer.toInt(slice);

                //System.out.print("prot id: " + proteinId + " offset: " + offset + " len: " + length);

                String peptideSequence = proteinCache.getPeptideSequence(proteinId, offset, length);
                //System.out.println("pep: " + peptideSequence);

                List<Integer> proteinIds = new ArrayList<Integer>();
                proteinIds.add(proteinId);
                IndexedSequence firstSequence = new IndexedSequence(0, seqMass, peptideSequence, "", "");
                firstSequence.setProteinIds(proteinIds);
                //set residues
                final String protSequence = proteinCache.getProteinSequence(proteinId);
                ResidueInfo residues = Util.getResidues(null, offset, length, protSequence);
                firstSequence.setResidues(residues);

                toInsert.add(firstSequence);

                //loop over more proteins ids until separator
                while (true) {
                    slice = Arrays.copyOfRange(data, i, i + 4);
                    i += 4;
                    proteinId = DynByteBuffer.toInt(slice);

                    if (proteinId == SEQ_SEPARATOR_INT) {
                        break;
                    } else {
                        proteinIds.add(proteinId);
                    }
                }

            } //end for every byte



        }

        /**
         * Parse the encoded sequences and return them in toInsert This is
         * version with multiple mass ranges (slightly slower as we need to
         * check every range if sequence qualifies)
         *
         * @param data the data row for all peptides grouped by a single key,
         * i.e. 1 ppm value
         * @param toInsert parsed sequences that qualify
         * @param minMasses used to eliminate not matching masses in the 1ppm
         * row
         * @param maxMasses used to eliminate not matching masses in the 1ppm
         * row
         */
        private void parseAddPeptideInfo(byte[] data, List<IndexedSequence> toInsert, float[] minMasses, float[] maxMasses) {

            int dataLength = data.length;

           // if (dataLength % 4 != 0) {
             //   throw new RuntimeException("Unexpected number of peptide items: " + dataLength);
            //}

            final int numRanges = minMasses.length;

            //go over every sequence in the data row
            for (int i = 0; i < dataLength;) {
                final int zeroeth = i;
                i += 4;
                final int first = i;
                i += 4;
                final int second = i;
                i += 4;
                final int third = i;
                i += 4;
                final int fourth = i;

                byte[] slice = Arrays.copyOfRange(data, zeroeth, first);
                float seqMass = DynByteBuffer.toFloat(slice);

                //check how the actual sequence mass fits in all mass ranges requested
                //we should bail out if pass all the ranges
                boolean greaterThanMax = true;
                //boolean lesserThanMin = true;
                boolean qualifiesRange = false; //check if qualifies in any range supplied
                for (int range = 0; range < numRanges; ++range) {
                    if (seqMass < maxMasses[range]) {
                        greaterThanMax = false;
                    }

                    //if (seqMass > minMasses[range]) {
                    //  lesserThanMin = false;
                    //}

                    if (seqMass >= minMasses[range]
                            && seqMass <= maxMasses[range]) {
                        qualifiesRange = true;
                    }

                    //if fits, bail out of this check earlier, otherwise check all ranges
                    if (qualifiesRange == true) {
                        break;
                    }

                    //check if makes any range

                }

                //since it's sorted
                //if current mass great than maxmass, break early - optimize
                if (greaterThanMax && !qualifiesRange) {
                    break;
                }

                //if (lesserThanMin && !qualifiesRange) {
                if (!qualifiesRange) {
                    //skip the sequence, it qualified the row, but not the actual peptide by precise mass 

                    //skip to next sequence
                    while (true) {
                        slice = Arrays.copyOfRange(data, i, i + 4);
                        i += 4;
                        int proteinId = DynByteBuffer.toInt(slice);
                        if (proteinId == SEQ_SEPARATOR_INT) {
                            break;
                        }
                    }
                    //go to next seq
                    continue;
                }


                //create sequence object to return

                slice = Arrays.copyOfRange(data, first, second);
                int offset = DynByteBuffer.toInt(slice);
                slice = Arrays.copyOfRange(data, second, third);
                int length = DynByteBuffer.toInt(slice);

                slice = Arrays.copyOfRange(data, third, fourth);
                int proteinId = DynByteBuffer.toInt(slice);

                //System.out.print("prot id: " + proteinId + " offset: " + offset + " len: " + length);

                String peptideSequence = proteinCache.getPeptideSequence(proteinId, offset, length);
                //System.out.println("pep: " + peptideSequence);

                List<Integer> proteinIds = new ArrayList<Integer>();
                proteinIds.add(proteinId);
                IndexedSequence firstSequence = new IndexedSequence(0, seqMass, peptideSequence, "", "");
                firstSequence.setProteinIds(proteinIds);
                //set residues
                final String protSequence = proteinCache.getProteinSequence(proteinId);
                ResidueInfo residues = Util.getResidues(null, offset, length, protSequence);
                firstSequence.setResidues(residues);

                toInsert.add(firstSequence);

                //loop over more proteins ids until separator
                while (true) {
                    slice = Arrays.copyOfRange(data, i, i + 4);
                    i += 4;
                    proteinId = DynByteBuffer.toInt(slice);

                    if (proteinId == SEQ_SEPARATOR_INT) {
                        break;
                    } else {
                        proteinIds.add(proteinId);
                    }
                }

            } //end for every byte

        }

        private byte[] getMergedData(byte[] seqData) {
            DynByteBuffer seqDataMerged = new DynByteBuffer();

            //to collapse multiple sequences into single one, with multproteins
            Map<String, List<IndexedSeqInternal>> temp =
                    new HashMap<String, List<IndexedSeqInternal>>();

            int dataLength = seqData.length;


           // if (dataLength % 4 != 0) {
             //   throw new RuntimeException("Unexpected number of peptide items: " + dataLength);
            //}

            for (int i = 0; i < dataLength;) {
                final int zeroeth = i;
                i += 4;
                final int first = i;
                i += 4;
                final int second = i;
                i += 4;
                final int third = i;
                i += 4;
                final int fourth = i;

                byte[] slice = Arrays.copyOfRange(seqData, zeroeth, first);
                float seqMass = DynByteBuffer.toFloat(slice);
                slice = Arrays.copyOfRange(seqData, first, second);
                int offset = DynByteBuffer.toInt(slice);
                slice = Arrays.copyOfRange(seqData, second, third);
                int length = DynByteBuffer.toInt(slice);
                slice = Arrays.copyOfRange(seqData, third, fourth);
                int proteinId = DynByteBuffer.toInt(slice);

                String peptideSequence = null;
                //System.out.print("prot id: " + proteinId + " offset: " + offset + " len: " + length);

                peptideSequence = proteinCache.getPeptideSequence(proteinId, offset, length);
                //System.out.println("pep: " + peptideSequence);

                final IndexedSeqInternal tempSequence = new IndexedSeqInternal(seqMass, offset, length, proteinId, null);

                List<IndexedSeqInternal> sequences = temp.get(peptideSequence);
                if (sequences == null) {
                    sequences = new ArrayList<IndexedSeqInternal>();
                    temp.put(peptideSequence, sequences);
                }
                sequences.add(tempSequence);

            }

            //sort sequences by masses  for search optimization
            List<IndexedSeqMerged> sortedMerged = new ArrayList<IndexedSeqMerged>();

            //group the same peptides from many proteins into single peptide with protein id list
            //sort them by mass for optimization
            for (String pepSeqKey : temp.keySet()) {
                //for each sequence str

                List<IndexedSeqInternal> sequences = temp.get(pepSeqKey);

                List<Integer> proteinIds = new ArrayList<Integer>();
                for (IndexedSeqInternal tempSeq : sequences) {
                    proteinIds.add((int) tempSeq.proteinId);
                }

                //make sure the 1st protein id is that of the first sequence
                IndexedSeqInternal firstSeq = sequences.get(0);

                IndexedSeqMerged merged = new IndexedSeqMerged(firstSeq.mass, firstSeq.offset, firstSeq.length, proteinIds);

                sortedMerged.add(merged);

            } //end of this sequence

            Collections.sort(sortedMerged);

            //write out grouped merged peptides to byte buf
            for (IndexedSeqMerged merged : sortedMerged) {
                //for each sequence str
                //System.out.println(merged);

                byte[] seqMassB = DynByteBuffer.toByteArray(merged.mass);
                seqDataMerged.add(seqMassB);

                byte[] seqOffsetB = DynByteBuffer.toByteArray(merged.offset);
                seqDataMerged.add(seqOffsetB);

                byte[] seqLengthB = DynByteBuffer.toByteArray(merged.length);
                seqDataMerged.add(seqLengthB);

                for (int proteinId : merged.proteinIds) {
                    byte[] proteinIdB = DynByteBuffer.toByteArray(proteinId);
                    seqDataMerged.add(proteinIdB);
                }
                //add separator
                seqDataMerged.add(SEQ_SEPARATOR);

            } //end of this sequence



            return seqDataMerged.getData();
        }
    }
}
/**
 * internal temporarly representation of indexed sequence before sequences are
 * merged into IndexedSequence
 *
 * two IndexedSeqInternal are equal if they sequence strings are equal
 *
 * two IndexedSeqInternal sequences are comparable by their mass
 *
 */
class IndexedSeqMerged implements Comparable<IndexedSeqMerged> {

    public IndexedSeqMerged(float mass, int offset, int length, List<Integer> proteinIds) {
        this.mass = mass;
        this.offset = offset;
        this.length = length;
        this.proteinIds = proteinIds;

    }
    float mass;
    int offset;
    int length;
    List<Integer> proteinIds;

    @Override
    public int compareTo(IndexedSeqMerged o) {
        if (mass < o.mass) {
            return -1;
        } else if (mass > o.mass) {
            return 1;
        } else {
            return 0;
        }
    }

    @Override
    public String toString() {
        return "IndexedSeqMerged{" + "mass=" + mass + ", offset=" + offset + ", length=" + length + ", proteinIds=" + proteinIds + '}';
    }
}
