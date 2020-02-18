/**
 * @file BlindSearchEngine.java
 * This is the source file for edu.scripps.pms.blindptm.BlindSearchEngine
 * @author Tao Xu
 * @date $Date: 2007/10/10 23:07:56 $
 */
package edu.scripps.pms.blindssptm;

import edu.scripps.pms.util.seq.Fasta;

import java.net.*;
import java.util.*;
import java.io.*;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.TimeUtils;
import edu.scripps.pms.util.io.SpectrumReader;
import edu.scripps.pms.util.io.MzxmlSpectrumReader;
import org.jdom.JDOMException;

class BlindSearchEngine {

    //public static final int DEFAULTNUMPEPTIDEHIT = 1000;
    //public static final int DEFAULTPEPTIDELENGTH = 40;
    private SearchParams params;
    private Searchable pdb;
    private MassCalculator mc;
    //private ProcessedPeakList ppl = null;
    public static final int HYPERGEOMETRYSCORE = 1;
    private static final String USAGE = "java BlindSearchEngine dbfile ms2file";
    private String hostName;
    private long startTime;
    private long endTime;
    private double maxPrcMass;
    private double minPrcMass;

    private static StringBuffer HEADER = new StringBuffer(1200);
      /*  "H\tSQTGeneratorVersion\t 0.1\n"+
        "H\tComment\tThe Department of Cell Biology\n"+
        "H\tComment\tThe Scripps Research Institute\n" +
        "H\tComment\tProbability Score and Confidence developed\n" +
        "H\tComment\tby Rovshan Sadygov and JR Yates, III\n" +
        "H\tComment\tRef. Anal. Chem., 2003\n"+
        "H\tComment\tCan use Cross-Correlation of J.Eng and J.Yates\n" +
        "H\tComment Ref. J. Am. Soc. Mass Spectrom., 1994, 4, p. 976\n"+
        "H\tComment Paralellization program is prolucid_submit\n";
      */
    private static DistributionCalculator dc;

    public BlindSearchEngine(SearchParams sp, Searchable db) { 
        this.params = sp;
        minPrcMass = params.getMinPrecursorMass(); 
        maxPrcMass = params.getMaxPrecursorMass();
        pdb = db;
        mc = new MassCalculator(sp);
        hostName = getHostName();
        getHeader();
    }

    private void getHeader() {
        HEADER.append("H\tSQTGenerator\tProLuCID\n");
        HEADER.append("H\tSQTGeneratorVersion\t0.1\n");
        HEADER.append("H\tDatabase\t" + params.getDbName()+ "\n");
        HEADER.append("H\tPrecursorMasses\t" + params.getPrecursorIsotope() + "\n");
        HEADER.append("H\tFragmentMasses\t" + params.getFragmentIsotope() + "\n");
    }
    public static String getHostName() {
        String hostName = null;
        try {
            hostName  = InetAddress.getLocalHost().getHostName();
            System.out.println("host name " + hostName);
            if(hostName.indexOf(".") != -1) {
                hostName = hostName.substring(0, hostName.indexOf("."));
            }
        } catch (UnknownHostException e) {
            e.printStackTrace();
        }
        return hostName;

    }
    public List<SearchResult> search(PeakList peaks) throws IOException {
        ArrayList<SearchResult> results = new ArrayList<SearchResult>();
        //PeptideHit.setFragMasses(mc.getFragMasses()); makes no big difference
        for(Iterator<Zline> zlines = peaks.getZlines(); zlines.hasNext();) {
            Zline z = zlines.next();
            double precMass = z.getM2z();
            int chargeState = z.getChargeState(); 
            //int numPeaks = peaks.numPeaks();
            startTime = System.currentTimeMillis();
      
            ProcessedPeakList ppl = new ProcessedPeakList(peaks, z, params, mc); 
            if(precMass >= minPrcMass && precMass <= maxPrcMass && ppl.getPTrue() < 0.5) {
                //SearchResult sr = ppl.search(pdb);//search(peaks, params, z);
                SearchResult sr = pdb.search(ppl);//search(peaks, params, z);
         
            //System.out.println("TimeUsed: " + timeUsedInSeconds);
                sr.calcScores();
                sr.setHostName(hostName);
            //System.out.println(sr.outputResults());
                //endTime = System.currentTimeMillis();
                //long timeUsedInSeconds = (endTime - startTime);
                sr.setSearchTime(System.currentTimeMillis() - startTime);
                results.add(sr);
            }
            //System.out.println("ChargeState: " + chargeState + "\tPrecursorM2z: "
            //       + precMass + "\tNumber of peptideHits: " + sr.getNumPeptideHit());
       	}
        if(params.getChargeDisambiguation()) {
            results = chargeDisambiguation(results);
        }
        
        return results;
    }

    public ArrayList<SearchResult> chargeDisambiguation(ArrayList<SearchResult> results) {

        ArrayList<SearchResult> newResults = new ArrayList<SearchResult>();
        SearchResult topPrimary = null;
        SearchResult topSp = null;

        for(Iterator<SearchResult> it = results.iterator(); it.hasNext();) {
            SearchResult sr = it.next();
            if(topSp == null) {
                topSp = sr;
                topPrimary = sr;
            } else {
                if(sr != null) {
                    if(topSp.getTopHit().getSecondaryScore() < sr.getTopHit().getSecondaryScore()) {
                        topSp = sr;
                    }
                    if(topPrimary.getTopHit().getPrimaryScore() < sr.getTopHit().getPrimaryScore()) {
                        topPrimary = sr;
                    }
                }
            }
        }
        if(topSp != null) {
            newResults.add(topPrimary);
            if(topSp != topPrimary) {
                newResults.add(topSp);
            } 
        }

        return newResults;

    }
    private static Iterator<PeakList> getSpectra(String file) throws IOException, JDOMException, Exception  {
        Iterator<PeakList> peaklists = null;
        if(file.endsWith(".mzXML")) {
            MzxmlSpectrumReader sr = new MzxmlSpectrumReader(file);
            ArrayList<PeakList> spectra = new ArrayList<PeakList>(20000);
            for(Iterator<MzxmlPeakList> it = sr.getSpectra(2); it.hasNext();) {
                spectra.add(it.next()); 
            }
            peaklists = spectra.iterator();
            sr.closeDataFile();
            sr = null;
        } else {
            SpectrumReader sr = new SpectrumReader(file, "ms2");
            peaklists = sr.getSpectraList().iterator();
            sr.closeDataFile();
            sr = null;
        }
    
        return peaklists;
    }

    public static void main(String args[]) throws Exception {
        try {
            String msFileName = args[0];
            SearchParams sp = new SearchParams("search.xml");
           
            // protein database 
            Searchable pd = null;
            if(sp.isDatabaseIndexed()) {
                pd = new IndexedProteinDatabase(sp.getDbName());
            } else {

                pd = new ProteinDatabase(sp.getDbName());
            }
            pd = new BlindProteinDatabase(sp.getDbName(), sp);
            dc = new DistributionCalculator();
            //ProteinDatabase pd = new IndexedProteinDatabase(args[0]);
            BlindSearchEngine se = new BlindSearchEngine(sp, pd);
            //SpectrumReader sr = new SpectrumReader(msFileName, "ms2");
            //String outFileName = msFileName.substring(0, msFileName.lastIndexOf(".ms2")) + ".sqt";
            String outFileName = msFileName.substring(0, msFileName.lastIndexOf(".")) + ".sqt";
            PrintWriter outFile = new PrintWriter(new BufferedWriter(new FileWriter(outFileName)));
 
            outFile.print(HEADER);
            //int numSearch = Integer.parseInt(args[2]);
            //System.out.println();
            int numPeaksMin = sp.getMinNumSpectra(); 
            int numPeaksMax = sp.getMaxNumSpectra(); 
            long beginTime = System.currentTimeMillis();
            int i = 0;
            TimeUtils timer = new TimeUtils();
            timer.startTiming();
            int numSpectraSearched = 0;
            int numSpectraIgnored = 0;
            int totalNumSpectra = 0;
            //for (Iterator<PeakList> it = sr.getSpectraList().iterator(); it.hasNext();) {
            for (Iterator<PeakList> it = getSpectra(msFileName); it.hasNext();) {
                i++;
                PeakList pl = it.next(); 
 
                totalNumSpectra++;
                int numPeaks = pl.numPeaks();
                if(numPeaks >= numPeaksMin && numPeaks <= numPeaksMax) {
                //timer.stopTiming(); 
                //System.out.println("Time used to read spectrum number " + i + ": " + timer.getTimeUsedMillis());
                //double prcMass = pl.getPrecursorMass();
                    List<SearchResult> results = se.search(pl);
                    for(Iterator<SearchResult> resultIt = results.iterator(); resultIt.hasNext();) {
                       
                        SearchResult rest = resultIt.next();
                        if(rest != null) { 
                            outFile.print(rest.outputResults());
                        }
                      
                        numSpectraSearched++;
                    //outFile.flush();
                    }
               
                    System.gc();
                } else {
                    numSpectraIgnored++;
                } 


                //timer.startTiming();
                //if (i == numSearch) {
                //    break;
               // }
            }
            long endTime = System.currentTimeMillis();
            long timeUsed = endTime - beginTime;
            outFile.close();
            System.out.println("NumSpectraSearched: " + numSpectraSearched + "\t numSpectraIgnored: " + numSpectraIgnored + "\ttotalNumSpectra: " + totalNumSpectra);
            System.out.println("Time used for search: " + timeUsed);
            System.out.println();
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
        }
    }
}


