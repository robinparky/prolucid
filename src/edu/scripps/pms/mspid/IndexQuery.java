package edu.scripps.pms.mspid;

import blazmass.dbindex.DBIndexer;
import blazmass.dbindex.DBIndexerException;
import blazmass.dbindex.IndexedProtein;
import blazmass.dbindex.IndexedSequence;
import blazmass.io.SearchParamReader;
import edu.scripps.pms.util.seq.Fasta;
import org.jdom.JDOMException;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

/**
 * Created by rpark on 3/20/16.
 */

public class IndexQuery {

    private static final Logger logger = Logger.getLogger(IndexQuery.class.getName());
    public static void main(String[] args) throws IOException, JDOMException, DBIndexerException {
        System.out.println("aaaa");

        String path = "/home/rpark/test_data/search/small";
        String filename = "search.xml";
        //SearchParams sp = new SearchParams(param);
        SearchParamReader paramReader = new SearchParamReader(path, filename);
        blazmass.io.SearchParams param = paramReader.getParam();
        //param.setIndexType(DBIndexer.IndexerMode.SEARCH_INDEXED);


        //final DBIndexer.IndexerMode indexerMode = sParam.isUseIndex()? DBIndexer.IndexerMode.SEARCH_INDEXED: DBIndexer.IndexerMode.SEARCH_UNINDEXED;
        final DBIndexer indexer = new DBIndexer(param);

        List<IndexedSequence> sequenceList = indexer.getSequences(2000f, 500f);

        for(Iterator<IndexedSequence> itr=sequenceList.iterator(); itr.hasNext(); ) {
            IndexedSequence seq = itr.next();

  //          if(seq.getWholeCleanSequence().contains("WHLKTEI"))
//                System.out.println("@@@@@============\t" + seq.getMass() + "++\t" + seq.getWholeCleanSequence() + "\t" + seq.getSequence() + "\t" + seq.getResLeft() + "\t" + seq.getResLeft());


        }



/*
        System.out.println(sequenceList.get(0));
        sequenceList = indexer.getSequences(1234.2f, 0.8f);
        System.out.println(sequenceList.get(0));

        sequenceList = indexer.getSequences(1234.2f, 0.8f);
        System.out.println(sequenceList.get(0));

        sequenceList = indexer.getSequences(1234.2f, 0.8f);
        System.out.println(sequenceList.get(0));

        sequenceList = indexer.getSequences(1234.2f, 0.8f);
        System.out.println(sequenceList.get(0));
*/

        /*
        for(Iterator<IndexedSequence> itr=sequenceList.iterator(); itr.hasNext(); ) {
            System.out.println(itr.next());
        }
        */



    }
}
