/**
 * @file SearchEngine.java
 * This is the source file for edu.scripps.pms.mspid.SearchEngine
 * @author Tao Xu
 * @date $Date: 2014/07/08 23:52:13 $
 */
package edu.scripps.pms.mspid;

import edu.scripps.pms.util.seq.Fasta;

import java.net.*;
import java.util.*;
import java.io.*;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.enzyme.Protease;
import edu.scripps.pms.util.TimeUtils;
import edu.scripps.pms.mspid.db.IndexedProteinDatabase;
import edu.scripps.pms.util.io.SpectrumReader;
import edu.scripps.pms.util.io.MzxmlSpectrumReader;
import org.jdom.JDOMException;


// original single process prolucid
class SearchEngine {

    //public static final int DEFAULTNUMPEPTIDEHIT = 1000;
    //public static final int DEFAULTPEPTIDELENGTH = 40;
    private SearchParams params;
    private String spectrumFile;
    private ProteinDatabase pdb;
    private MassCalculator mc;
    //private ProcessedPeakList ppl = null;
    public static final int HYPERGEOMETRYSCORE = 1;
    private static final String USAGE = "Usage: java SearchEngine ms2fileOrmzXMLfile";
    private String hostName;
    private long startTime;
    private long endTime;
    private double maxPrcMass;
    private double minPrcMass;

    //private static StringBuffer HEADER = new StringBuffer(1200);
    private static DistributionCalculator dc;
    private static boolean [] searchedScans = new boolean[1000000];
    protected static StringBuffer logbuffer = new StringBuffer(2000);

    public SearchEngine(SearchParams sp, ProteinDatabase db, String msFileName) { 
        this.params = sp;
        spectrumFile = msFileName;
        minPrcMass = params.getMinPrecursorMass(); 
        maxPrcMass = params.getMaxPrecursorMass();
        pdb = db;
        mc = new MassCalculator(sp);
        hostName = getHostName();
    }

    private String getHeader() {
        StringBuffer HEADER = new StringBuffer(1200);
        HEADER.append("\n\nH\tSQTGenerator\tProLuCID\n");
        HEADER.append("H\tSQTGeneratorVersion\t1.2.4\n");
        HEADER.append("H\tComment ProLuCID is developed in the Yates laboratory at The Scripps Research Institute, La Jolla, CA\n");
        HEADER.append("H\tComment ProLuCID ref. Xu T, Venable JD, Park SK, Cociorva D, Lu B, Liao L, Wohlschlegel J, Hewel J, Yates JR 3rd\n");
        HEADER.append("H\tComment ProLuCID ref. ProLuCID, a fast and sensitive tandem mass spectra-based protein identification program.\n");
        HEADER.append("H\tComment ProLuCID ref. MOL CELL PROTEOMICS vol. 5(10): S174-S174 671 Suppl. S OCT 2006\n");
        HEADER.append("H\tComment Paralellization Program using PBS is submit_prolucid\n");
        HEADER.append("H\tComment Please send bug report or comments to Tao Xu by email taoxu@scripps.edu\n");
        HEADER.append("H\tHostName\t" + getHostName() + "\n");
        HEADER.append("H\tSpectrumFile\t" + spectrumFile + "\n");
        HEADER.append("H\tDatabase\t" + params.getDbName()+ "\n");
        HEADER.append("H\tDBLocusCount\t" + pdb.getNumSequences() + "\n");
        HEADER.append("H\tNumOutput\t" + params.getNumOutput() + "\n");
        HEADER.append("H\tPrecursorMasses\t" + params.getPrecursorIsotope() + "\n");
        HEADER.append("H\tFragmentMasses\t" + params.getFragmentIsotope() + "\n");
        int numIsotopicPeaks = params.getNumIsotopicPeaks();
        if(numIsotopicPeaks == 0) {
            HEADER.append("H\tHighPrecursorMassTolerance\t" + params.getHighPrecursorTolerance() + "\n");
            HEADER.append("H\tLowPrecursorMassTolerance\t" + params.getLowPrecursorTolerance() + "\n");
        } else {
            HEADER.append("H\tNumPrecursorIsotopicPeaks\t" + numIsotopicPeaks + "\n");
            HEADER.append("H\tPrecursorMassTolerance\t" + params.getPrecursorTolerance() + "\n");
        }
        HEADER.append("H\tFragmentMassTolerance\t" + params.getFragmentTolerance() + "\n");
        if(params.getStaticNTermMod() != 0) {
            HEADER.append("H\tNTermStaticMod\t" + params.getStaticNTermMod() + "\n");
        }
        if(params.getStaticCTermMod() != 0) {
            HEADER.append("H\tCTermStaticMod\t" + params.getStaticCTermMod() + "\n");
        }
        for(Iterator<Modification> it = params.getStaticMods(); it.hasNext();) {
            Modification m = it.next();
            HEADER.append("H\tStaticMod\t" + m.getResidue()+ "=" + mc.getPrecursorMass(m.getResidue()) + "\n");
        }
        
        if(params.getNumNTermDiffMods() > 0) {
            HEADER.append("H\tNTermDiffMod");
            for(Iterator<TerminalModification> it =  params.getNTermDiffMods(); it.hasNext();) {
                TerminalModification m = it.next();
                HEADER.append("\t" + m.getMassShift());
            }
            HEADER.append("\n");
        }
 
        if(params.getNumCTermDiffMods() > 0) {
            HEADER.append("H\tCTermDiffMod");
            for(Iterator<TerminalModification> it =  params.getCTermDiffMods(); it.hasNext();) {
                TerminalModification m = it.next();
                HEADER.append("\t" + m.getMassShift());
            }
            HEADER.append("\n");
        }

        HEADER.append("H\tMaxNumInternalDiffModsPerPeptide\t" + params.getMaxAlter() + "\n");
        if(params.getMaxAlter() > 0) {
            for(Iterator<DiffMod> it = params.getDiffMods(); it.hasNext();) {
                DiffMod m = it.next();
                //HEADER.append("H\tInternalDiffMod\t" + m.toString() + "\n");
                HEADER.append("H\tDiffMod\t" + m.toString() + "\n");
            }
        }
        int enzymespec = params.getEnzymeSpecificity();
        if(enzymespec > 0) {
            HEADER.append("H\tEnzymeSpecificity\t" + params.getEnzymeSpecificity() + "\n");
            Protease p = params.getEnzyme();
            boolean isC = p.getType();

            HEADER.append("H\tEnzymeName\t" + p.getName() + "\n");
            if(isC) {
                HEADER.append("H\tEnzymeEnd\tCTerm\n");
            } else {
                HEADER.append("H\tEnzymeEnd\tNTerm\n");
            }
            HEADER.append("H\tEnzymeResidues\t");
            for(int i = 0; i < MassSpecConstants.NUMCHARS; i++) {
                if(p.isDigestable((char)i)) {
                    HEADER.append((char)i); 
                }
            }
            HEADER.append("\n");

        } else {
            HEADER.append("H\tEnzymeSpecificity\tNo_Enzyme\n");
        }
        return HEADER.toString();
    }
    public static String getHostName() {
        String hostName = null;
        try {
            hostName  = InetAddress.getLocalHost().getHostName();
            //logbuffer.append("host name " + hostName + "\n");
            if(hostName.indexOf(".") != -1) {
                hostName = hostName.substring(0, hostName.indexOf("."));
            }
        } catch (UnknownHostException e) {
            //e.printStackTrace();
        }
        if(hostName == null) {
            hostName =  "UnknownHost";
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
            if(chargeState < params.getMinPrecursorCharge() || chargeState > params.getMaxPrecursorCharge()) {
                continue;
            }
            startTime = System.nanoTime();
      
            ProcessedPeakList ppl = new ProcessedPeakList(peaks, z, params, mc); 
            
            //logbuffer.append("FNumPeaks: " + ppl.getFinalNumPeaks() +  "\tScan#: " + peaks.getLoscan() + "\t+" + chargeState + "\tNumPeaks: " + peaks.numPeaks() + "\tm2z: " + precMass + "\tptrue: " + ppl.getPTrue() );
            logbuffer.append("FNumPeaks: " + ppl.getFinalNumPeaks() + "\t+" + chargeState + "\tNumPeaks: " + peaks.numPeaks() + "\tM+H+: " + precMass + "\tptrue: " + ppl.getPTrue() );

            if(precMass >= minPrcMass && precMass <= maxPrcMass && ppl.getPTrue() < 0.5) {
                //SearchResult sr = ppl.search(pdb);//search(peaks, params, z);
                SearchResult sr = pdb.search(ppl);//search(peaks, params, z);
         
            //System.out.println("TimeUsed: " + timeUsedInSeconds);
                sr.calcScores();
                sr.setHostName(hostName);
            //System.out.println(sr.outputResults());
                //endTime = System.nanoTime();
                //long timeUsedInSeconds = (endTime - startTime);
                sr.setSearchTime((System.nanoTime() - startTime)/1000000);
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
    private static ArrayList<PeakList> getSpectra(String file) throws IOException, JDOMException, Exception  {
        ArrayList<PeakList> peaklists = null;
        if(file.endsWith(".mzXML")) {
            MzxmlSpectrumReader sr = new MzxmlSpectrumReader(file);
            ArrayList<PeakList> spectra = new ArrayList<PeakList>(20000);
            for(Iterator<MzxmlPeakList> it = sr.getSpectra(2); it.hasNext();) {
                spectra.add(it.next()); 
            }
            peaklists = spectra;
            sr.closeDataFile();
            sr = null;
        } else {
            SpectrumReader sr = new SpectrumReader(file, "ms2");
            peaklists = sr.getSpectraList();
            sr.closeDataFile();
            sr = null;
        }
    
        return peaklists;
    }

    // currently not used, because if the ms2 file names before the _ is not unique,
    // it will cause trouble
    private static void markSearchedScans(String ms2file) throws IOException {
        //System.out.println("Checking searched spectra for " + ms2file);
        if(ms2file == null) return;
        int posdot = ms2file.lastIndexOf(".");
        int posunderscore = ms2file.lastIndexOf("_");
        String longsqtfile = "";
        String shortsqtfile = "";
    
        if(posdot > 0) longsqtfile = ms2file.substring(0, posdot) + ".sqt";
        if(posunderscore > 0) shortsqtfile = ms2file.substring(0, posunderscore) + ".sqt"; 
                
        File sqtfile = null;
        File shortfile = new File(shortsqtfile);
        File longfile = new File(longsqtfile);
         
        if(longfile.exists()) {
            sqtfile = longfile;
        } else if(shortfile.exists()) {
            sqtfile = shortfile;
        }
        if(sqtfile != null) {
            BufferedReader br = new BufferedReader(new FileReader(sqtfile));
            String line = null;
            while((line = br.readLine()) != null) {
                if(line.startsWith("S\t")) {
                    String [] arr = line.split("\t");
                    int scannum = Integer.parseInt(arr[2]);
                    searchedScans[scannum] = true;
                    //System.out.println("Searched scan: " + scannum);
                }
            }
            br.close();
        }
    }
    public static void main(String args[]) throws Exception {
        try {
            String msFileName = args[0];
            String searchxml = "search.xml";
            if(args.length > 1) {
                searchxml = args[1];
            }
            
            // if file name before the last _ is not unique, may cause trouble
            //markSearchedScans(msFileName); //it is more likely to cause trouble. not use any more

            SearchParams sp = new SearchParams(searchxml);
            SearchResult.setNumOutput(sp.getNumOutput());           
 
            ProteinDatabase pd = null;
            if(sp.isDatabaseIndexed()) {
                pd = new IndexedProteinDatabase(sp.getDbName());
            } else {

                pd = new ProteinDatabase(sp.getDbName());
            }

            pd.setMaxMisCleavage(sp.getMaxInternalMisCleavage());
            pd.setProtease(sp.getProtease());
            pd.setEnzymeSpecificity(sp.getEnzymeSpecificity());
/*
            long totalmemory = Runtime.getRuntime().totalMemory();
            long usedmemory  = totalmemory - Runtime.getRuntime().freeMemory();

            logbuffer.append("Total memory\t" + totalmemory + "\n");
            logbuffer.append("Memory used after loading the database\t" + usedmemory + "\n");
            logbuffer.append("Host name " + getHostName() + "\n");*/
            //ProteinDatabase pd = new IndexedProteinDatabase(args[0]);
            SearchEngine se = new SearchEngine(sp, pd, msFileName);
            //SpectrumReader sr = new SpectrumReader(msFileName, "ms2");
            //String outFileName = msFileName.substring(0, msFileName.lastIndexOf(".ms2")) + ".sqt";
            String outFileName = msFileName.substring(0, msFileName.lastIndexOf(".")) + ".sqt";
            // default buffer size is 8192, 40960 is 5X the default buffer size
            //PrintWriter outFile = new PrintWriter(new BufferedWriter(new FileWriter(outFileName, true), 40960));
            PrintWriter outFile = new PrintWriter(new BufferedWriter(new FileWriter(outFileName)));
 
            outFile.print(se.getHeader());
            //int numSearch = Integer.parseInt(args[2]);
            //System.out.println();
            int numPeaksMin = sp.getMinNumSpectra(); 
            int numPeaksMax = sp.getMaxNumSpectra(); 
            long beginTime = System.nanoTime();
            int i = 0;
            TimeUtils timer = new TimeUtils();
            timer.startTiming();
            int numSpectraSearched = 0;
            int numSpectraIgnored = 0;
            int totalNumSpectra = 0;
            //for (Iterator<PeakList> it = sr.getSpectraList().iterator(); it.hasNext();) {
            ArrayList<PeakList> spectra = getSpectra(msFileName);
            int numSpectraInFile =  spectra.size();
            /*
            logbuffer.append("Total Memory after loading the search engine\t" + (Runtime.getRuntime().totalMemory() + "\tFree memory\t" + Runtime.getRuntime().freeMemory()) + "\n");
            logbuffer.append("Memory used after loading the search engine\t" + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) + "\n");
            logbuffer.append("TotalNumSpectra\t" + numSpectraInFile + "\n");
            //logbuffer.append("host name " + getHostName() + "\n");

            logbuffer.append("NumSearched\t" + totalNumSpectra + "\t" + numSpectraInFile + "\n" );
            System.out.print(logbuffer.toString()); 
            logbuffer = new StringBuffer(2000);
            */

            
            dc = new DistributionCalculator();

            Iterator<PeakList> it = spectra.iterator();
            while(it.hasNext()) {
                i++;
                PeakList pl = it.next(); 

            int tempnumprev = pl.numPeaks();
//pl.processHcdSpectrum();
//System.out.println("With HCD processing, number of peaks before: " + tempnumprev + " and after: " +  pl.numPeaks());
 
                totalNumSpectra++;
                //String logline = "NumSearched\t" + totalNumSpectra + "\t" + numSpectraInFile + "\t" ;
              //  String logline = "NumSearched\t" + totalNumSpectra + "\t" + numSpectraInFile + "\tScan#: " + pl.getLoscan() + "\tNumPeaks: " + tempnumprev + "\t";
              //  logbuffer.append(logline);
                /*if(totalNumSpectra%20 == 0) {

                    //System.out.print("NumSearched\t" + totalNumSpectra + "\t" + numSpectraInFile + "\t" );
                    System.out.print(logbuffer.toString()); 
                    logbuffer = new StringBuffer(2000);
                } */

                int numPeaks = pl.numPeaks();
                int scan = pl.getLoscan();
                //if(numPeaks >= numPeaksMin && numPeaks <= numPeaksMax && !searchedScans[scan]) {
                if(numPeaks >= numPeaksMin && numPeaks <= numPeaksMax) {
                //timer.stopTiming(); 
                //System.out.println("Time used to read spectrum number " + i + ": " + timer.getTimeUsedMillis());
                //double prcMass = pl.getPrecursorMass();
                    List<SearchResult> results = se.search(pl);
                    if(results.size() > 0) {
                        numSpectraSearched++;
                    } else {
                        numSpectraIgnored++;
                    }

                    for(Iterator<SearchResult> resultIt = results.iterator(); resultIt.hasNext();) {
                       
                        SearchResult rest = resultIt.next();
                        if(rest != null) { 
                            outFile.print(rest.outputResults());
                        }
                      
                        //outFile.flush();
                    }
               
                } else {
                    numSpectraIgnored++;
                } 

                //System.out.println(""); // end the log line for this spectrum
             //   logbuffer.append("\n");
                System.gc();
                //timer.startTiming();
                //if (i == numSearch) {
                //    break;
               // }
            } 

            long endTime = System.nanoTime();
            long timeUsed = (endTime - beginTime)/1000000;
            outFile.close();
            /*
            logbuffer.append("Total Memory after the search\t" + (Runtime.getRuntime().totalMemory() + "\tFree memory\t" + Runtime.getRuntime().freeMemory()) + "\n");
            logbuffer.append("Memory used after the search\t" + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) + "\n");
            logbuffer.append("NumSpectraSearched: " + numSpectraSearched + "\t numSpectraIgnored: " + numSpectraIgnored + "\ttotalNumSpectra: " + totalNumSpectra + "\n");
            logbuffer.append("Time used for search: " + timeUsed + "\n");
            System.out.println(logbuffer.toString());*/
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
            System.out.println("Problematic run on " + getHostName());
        }
    }
}


