/**
 * @file SearchResult.java
 * This is the source file for edu.scripps.pms.blindptm.SearchResult
 * @author Tao Xu
 * @date $Date
 */
package edu.scripps.pms.blindssptm;




import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.TimeUtils;
import edu.scripps.pms.util.enzyme.Protease;
import java.util.*;
import java.net.*;
import java.text.DecimalFormat;

public class SearchResult {
    /* for binomial and poisson
    private final int NUM = 36;
    private int [] freq = new int[NUM];
    private int [] freqall = new int[NUM];
    */
    private static int NUMSCORED;
    private static final int NUMFINALRESULT = 5;
    private static final char DELIMITER = '\t';
    //private LinkedList<PeptideHit> hits = new LinkedList<PeptideHit>(); 
    private ArrayList<PeptideHit> hits; 
    private List<ScoredPeptideHit> scoredHits;
    private List<ScoredPeptideHit> finalResult = null; 
    private boolean isSorted = true;
//    private int totalPeptideLength = 0;
    private int numPeaks = 0;  // total number of peaks from all matched peptides
    private int totalPeptideLength = 0;
    private int numPeaksMatched = 0;
    private final int minNumPeaksMatched;
    private final int minPeptideLength;
    private int chargeState;
    //private final int effectChargeState; // chargeState of ppl
    private ProcessedPeakList ppl;
    private SearchParams params;
    private int numPeptidesMatched = 0;
    private MassCalculator mc;
    private long searchTime;
    private String hostName;
    private PeakList peaks;
    private Protease protease;
    private int enzymeSpecificity;
    public static DecimalFormat threeDigits = new DecimalFormat("0.000");
    private static DecimalFormat fourDigits = new DecimalFormat("0.0000");
    private static DecimalFormat twoDigits = new DecimalFormat("0.00");
    private ScoredPeptideHit bestSecondaryScore;
    private boolean sortByXcorr = false;
    private int primaryScoreType;     
    private int secondaryScoreType;     
    private double primaryScoreDeviation = 1;
    private double primaryScoreMean = 0;
//    private static DistributionCalculator dc = new DistributionCalculator();

    public SearchResult(ProcessedPeakList ppl) {
        this.ppl = ppl;
        params = ppl.getSearchParams();
        protease = params.getProtease();
        
        enzymeSpecificity = params.getEnzymeSpecificity();
        chargeState  = ppl.getZline().getChargeState();       
        //effectChargeState  = ppl.getZline().getChargeState() > 2? 3 : 2;       
        peaks = ppl.getPeakList();
        mc = ppl.getMassCalculator();
        minNumPeaksMatched = params.getMinMatch();
        minPeptideLength = ppl.getSearchParams().getMinimumPeptideLength() - 2;
        NUMSCORED = ppl.getSearchParams().getCandidatePeptideThreshold();
        hits = new ArrayList<PeptideHit>(NUMSCORED+1); 
        scoredHits = new ArrayList<ScoredPeptideHit>(NUMSCORED);
        primaryScoreType = ppl.getSearchParams().getPrimaryScoreType(); 
        secondaryScoreType = ppl.getSearchParams().getSecondaryScoreType(); 
        sortByXcorr = (primaryScoreType != 0); // ie, primaryScoreType not probability score
    }
    public void setHostName(String host) {
        hostName = host;
    }
    
    public void setSearchTime(long timeUsedInMilliSeconds) {
        searchTime = timeUsedInMilliSeconds;
    }
    protected void ignoreModifiedPeptideHit() {
        ScoredPeptideHit  topModified = null;
        ScoredPeptideHit topUnModified = null;
        for(Iterator<ScoredPeptideHit> it = scoredHits.iterator(); it.hasNext();) {
            ScoredPeptideHit hit = it.next();
            if(hit.isModified()) {
                if(topModified == null) {
                    topModified = hit;
                }
            } else {
                if(topUnModified == null) {
                    topUnModified = hit;
                }            
            }
            if(topModified != null && topUnModified != null) {
                break;
            }
        } 
        if(topModified != null && topUnModified != null) {
            //if(topModified.isIgnorable(topUnModified)) {}
            // need to implement the criteria to determine ignore or not
        } 
    } 
    public List<ScoredPeptideHit> getFinalResult() {
       

        ArrayList<ScoredPeptideHit> finalList = new ArrayList<ScoredPeptideHit>();
        if(scoredHits.size() == 0) {
            return finalList;
        }
        String topseq = scoredHits.get(0).getOriginalSequence();

        int numAdded = 0;
        for(int i = 0; i < scoredHits.size(); i++) {
            ScoredPeptideHit sph = scoredHits.get(i);
            finalList.add(scoredHits.get(i));
            numAdded++;
//System.out.println(topseq + "\t" + sph.getOriginalSequence() + "\t" + NUMFINALRESULT + "\t" + sph.getPrimaryRank() + "\t" + !topseq.equals(sph.getOriginalSequence()));
//System.out.println(NUMFINALRESULT + "\t" + i + "\t" + (i < NUMFINALRESULT) + "\t" +  !topseq.equals(sph.getOriginalSequence()));
            if(numAdded >= NUMFINALRESULT && (!topseq.equals(sph.getOriginalSequence()))) {
                break; 
            }
        }
        //finalResult.add(bestSecondaryScore);
        return finalList;

        // the following was trying to do something special with modifcation hits
        /*
        ArrayList<ScoredPeptideHit> finalList = new ArrayList<ScoredPeptideHit>();
        for(int i = 0; finalList.size()< NUMFINALRESULT && i < scoredHits.size(); i++) {
            ScoredPeptideHit hit = scoredHits.get(i);
            if(!ignoreModifiedPeptideHit() || !hit.isModified()) {
                finalList.add(scoredHits.get(i));
            } 
        }
        //finalResult.add(bestSecondaryScore);
        if(finalList.size() > 0) {
            finalList.get(0).setPrimaryRank(1);
        }
        return finalList;
        */
    }
    public double getDeltaCn(ScoredPeptideHit sph) {
        /*
        if(tscore) {
            // use t score to replace deltaCN
            return (sph.getPrimaryScore() - primaryScoreMean)/primaryScoreDeviation;
        } else {
            double bestScore = finalResult.get(0).getPrimaryScore();       
            if(bestScore == 0) {
                return 0;
            }
            return (bestScore - sph.getPrimaryScore())/bestScore;
        }
        */
        double bestScore = finalResult.get(0).getPrimaryScore();       
        if(bestScore == 0) {
            return 0;
        }
        return (bestScore - sph.getPrimaryScore())/bestScore;

    } 
    public ScoredPeptideHit getTopHit() {
        if(finalResult != null && finalResult.size() > 0) {
            return finalResult.get(0);
        } else {
            return null;
        }
    }
    public int getNumPeptidesMatched() {
        return numPeptidesMatched;
    }

    public ProcessedPeakList getProcessedPeakList() {
        return ppl;
    }
    public String outputResults() {
        //calcScores();
        
        StringBuffer result = new StringBuffer(8000); 
        appendSline(result);
        appendMlines(result); // L lines are also appended
        result.append('\n');
/*
System.out.println("numPeaks: " + numPeaksMatched + "\tnumPeptides: " + numPeptidesMatched);
double lamda = numPeaksMatched/(numPeptidesMatched+0.0); 
System.out.println("lamda: \t" + lamda);
System.out.println("numPeaksMatched\tlength " + NUM/2 + "\tallLength\t");

for(int i = 0; i < freq.length; i++) {
    System.out.println(i + "\t" + freq[i] + "\t" + freqall[i]);
}
*/
        return result.toString();
        
    }
    

    private void appendMlines(StringBuffer result) {
        if(finalResult == null || finalResult.size() < 1) { return; }
        for(ScoredPeptideHit p : finalResult) {
            appendMline(result, p);
        }
    }   
    private void appendMline(StringBuffer result, ScoredPeptideHit sph) {
        // Format of M line:
        // "M\tRankByXcorr\tRankBySp\tCalculatedMass\tDeltaCN\tXcorr\t
        // Sp\tNumMatchedIons\tNumExpectedIons\tPeptideSequence\tValidateState\n"
        if (sph == null) return; 
        result.append("M");
        result.append(DELIMITER);
        result.append(sph.getPrimaryRank());
        result.append(DELIMITER);
        result.append(sph.getSecondaryRank());
        result.append(DELIMITER);
        
        result.append(fourDigits.format(sph.getTheorMass()));
        result.append(DELIMITER);
        result.append(fourDigits.format(getDeltaCn(sph))); //deltaCn 
        result.append(DELIMITER);
        result.append(fourDigits.format(sph.getPrimaryScore())); // XCorr
        result.append(DELIMITER);
        result.append(threeDigits.format(sph.getSecondaryScore())); // XCorr
        result.append(DELIMITER);
        result.append(sph.getNumPeaksMatched()); //numPeaksMatched 
        result.append(DELIMITER);
        result.append(sph.getNumPeaks()); //numPeaks
        result.append(DELIMITER);
        result.append(sph.getExtendedSequence()); // Peptide sequence
        result.append(DELIMITER);

        result.append("U\n"); // Validation state

        int locusType = params.getLocusType(); 
        //append Llines 
        for(PeptideHit p : sph.getPeptideHits()) {
            if(locusType == 1) {
                result.append("L\t" + p.getSequestLikeAccession() + "\n");
            } else {
                result.append("L\t" + p.getAccession() + "\n");
            }
           
        }

    }
    private void appendSline(StringBuffer result) {

        // Format of S line:
        // "S\tLowScan\tHighScan\tChargeState\tProcessTime\tSever\t
        // ObservedMass\tTotalIntensity\tLowestSp\tNumSequencesMatched
        result.append("S");
        result.append(DELIMITER);
        result.append(ppl.getPeakList().getHiscan());
        result.append(DELIMITER);
        result.append(ppl.getPeakList().getLoscan());
        result.append(DELIMITER);
        result.append(chargeState);
        result.append(DELIMITER);
        result.append(searchTime);
        result.append(DELIMITER);
        result.append(hostName);
        result.append(DELIMITER);
        result.append(fourDigits.format(ppl.getZline().getM2z()));
        result.append(DELIMITER);
        result.append(twoDigits.format(peaks.getTotalIntensity()));
        result.append(DELIMITER);
        result.append(fourDigits.format(ppl.getPTrue())); // for ptrue 
        result.append(DELIMITER);
        result.append(numPeptidesMatched);
        result.append("\n");
    }
    public void calcScores() {
        
        if(hits.size() > 0) {
            scoredHits = getTopHits(NUMSCORED);

            //numPeaks = ppl.getNumFragBins();
            //numPeaksMatched = ppl.getNumTrues();

            int prob = (int)(ppl.getPTrue()*100 + 0.5);
        //System.out.println("prob: " + prob);
            for(ScoredPeptideHit p : scoredHits) {
/*
                float score = ScoreCalculator.hypergeometry(
                        numPeaks, numPeaksMatched, 
                   p.getNumPeaks(), p.getNumPeaksMatched()); 
                score =  -(float)Math.log(score);
                p.setScore(score);
                p.setSecondaryScore(score);

                int n = p.getNumPeaks();
                int k = p.getNumPeaksMatched();
                double score = DistributionCalculator.getBinomialSum(prob, n, k); 
                score = -Math.log10(score);
*/
                
                double score = -Math.log10(DistributionCalculator.getBinomialSum
                                 (prob, p.getNumPeaks(), p.getNumPeaksMatched()));
 
                p.setPScore(score); // set probability score
                ppl.correlation(p);
 
            }
            
            //ignoreModifiedPeptideHit();
 
            // currently sorted by probability score
            setPScoreRank(scoredHits);
            if(sortByXcorr) {
                bestSecondaryScore = scoredHits.get(0);
                Collections.sort(scoredHits);
                finalResult = getFinalResult(); // after sort by XCorr
                 
            } else {
                finalResult = getFinalResult();  // before sort by SCorr
                Collections.sort(scoredHits);
                bestSecondaryScore = scoredHits.get(0);
            }
            // now sorted by XCorr
            setXCorrRank(scoredHits); 
            finalResult.add(bestSecondaryScore);

            calcPrimaryScoreDeviation(); //temperally added to use tscore as primary score 

            for(ScoredPeptideHit p : scoredHits) {
               //System.out.println(p.getXCorr());
                p.setZscore((p.getPrimaryScore() - primaryScoreMean)/primaryScoreDeviation);
              //p.setXCorr((p.getPrimaryScore() - primaryScoreMean)/primaryScoreDeviation);
            }
        }
    }
    private void calcPrimaryScoreDeviation() {
        List<ScoredPeptideHit> resultList = scoredHits; 
        //List<ScoredPeptideHit> resultList = finalResult; 
        Iterator<ScoredPeptideHit> it = resultList.iterator();
        int numElement = 0;
        String topseq = null;
        while(it.hasNext()) {
            ScoredPeptideHit p = it.next();
            if(topseq == null) {
                topseq = p.getOriginalSequence();
            }
            if(p.getPrimaryRank() != 1 && (!topseq.equals(p.getOriginalSequence()))) {
                primaryScoreMean += p.getPrimaryScore();
//System.out.println("PrimaryScore: " + p.getPrimaryScore());
                numElement++;
            }
        }
        if(numElement > 2) { 
            primaryScoreMean /= numElement;
            it = resultList.iterator();
            while(it.hasNext()) {
                ScoredPeptideHit p = it.next();
                if(p.getPrimaryRank() > 1 && (!topseq.equals(p.getOriginalSequence()))) {
                    double diff = p.getPrimaryScore() - primaryScoreMean;
                    primaryScoreDeviation  += diff*diff;
                }
            }
            primaryScoreDeviation = Math.sqrt(primaryScoreDeviation/(numElement-1));
            if(primaryScoreDeviation == 0) {
                primaryScoreDeviation = 1;
            }
        } 

//       System.out.println("Number of Elements for zscore: " + numElement);        
       // System.out.println("PrimaryScoreMean: " + primaryScoreMean + "\tPrimaryScoreDeviation: " + primaryScoreDeviation);
            
    }
    private void setPScoreRank(List<ScoredPeptideHit> scoredList) { 

        int rank = 0;
        double lastScore = 100;
        for(ScoredPeptideHit p : scoredList) {
            double pScore = p.getPScore();
            if(pScore < lastScore) {
                rank++;
                lastScore = pScore; 
            } 
            p.setPScoreRank(rank);
        }
    }
    private void setXCorrRank(List<ScoredPeptideHit> scoredList) {

        int rank = 0;
        double lastScore = 100;
        for(ScoredPeptideHit p : scoredList) {
            double xcorr = p.getXCorr();
            if(xcorr < lastScore) {
                rank++;
                lastScore = xcorr; 
            } 
            p.setXCorrRank(rank);
        }

    }
    private void chebyshev() {
        int numPeptideHits = hits.size();
        // j in CHEBYSHEV
        int avgNumPeaksMatched = Math.round(numPeaksMatched/(float)numPeptideHits);
        // k in CHEBYSHEV
        int avgPeptideLength = Math.round(totalPeptideLength/(float)numPeptideHits); 
                            
        int avgPeaks = 2*(avgPeptideLength-1)*(chargeState-1);
        float avgHypergeometry = ScoreCalculator.hypergeometry(
                 numPeaks, numPeaksMatched, avgPeaks, avgNumPeaksMatched);
        float avgScore = -(float)Math.log(avgHypergeometry);


        // calculate the best five scores here 
    }

    public void addPeptideHit(ModifiedPeptideHit p) {

        if(isValidHit(p)) {
            // check enzyme specificity 
            for(Iterator<ModifiedPeptideHit> it = p.getAllModifiedPeptideHits(); it.hasNext();) {
           //System.out.println("in add modified peptide hit"); 
                addFinalPeptideHit(it.next());
            }
        }
    }

    // for semi-blind modification search
    public void addPeptideHit(Peptide p, DiffMod m) {
        
        if(isValidHit(p)) {
            int length = p.getLength();
            for(int i = 0; i < length; i++) {
                ModifiedPeptideHit mph = new ModifiedPeptideHit(p);
                if(p.byteAt(i) == 'C') {  // temp added for SS search
                    mph.setDiffMod(i, m);
                    addFinalPeptideHit(mph);
                }
            }
        }
    }

    private boolean isValidHit(Peptide p) {
        int start = p.getStart();
        int end = p.getEnd();
        if((end-start) > minPeptideLength) {
            // check enzyme specificity 
            if(protease == null || protease.checkEnzymeSpecificity(p.getParent(), start, end) >= enzymeSpecificity) { 
                return true;
            }
        }
        return false;

    }
    private boolean isValidHit(PeptideHit p) {
        int start = p.getStart();
        int end = p.getEnd();
        if((end-start) > minPeptideLength) {
            // check enzyme specificity 
            if(protease == null || protease.checkEnzymeSpecificity(p.getParent(), start, end) >= enzymeSpecificity) { 
                return true;
            }
        }
        return false;

    }
    public void addPeptideHit(PeptideHit p) {
        if(isValidHit(p)) {
            addFinalPeptideHit(p);
        }
    }
    /**
     * @param cutoff - minimum number of peaks matched specified in search parameters
     */ 
    private void addFinalPeptideHit(PeptideHit p) {
        //numPeaks += (p.getLength()-1) * 2; 
        int start = p.getStart();
        int end = p.getEnd();
                //PeptideHit p = new PeptideHit(f, start, end);
                //int numMatched = ppl.calcNumPeaksMatched(p);

        p.setProcessedPeakList(ppl);
        p.calcNumPeaksMatched();
             
        int numMatched = p.getNumPeaksMatched();
        numPeptidesMatched++;
        //numPeaks += counts[1]; 
        numPeaks += p.getNumPeaks(); 
        totalPeptideLength += (end - start + 1);
        numPeaksMatched += p.getNumPeaksMatched();
        if (numMatched >= minNumPeaksMatched) {
        //    addPeptideHit(p);
            double prob = p.getProbability();
            int lastIndex = hits.size() - 1;
            int index = findIndex(prob, 0, lastIndex); 
            if(index < NUMSCORED) {
                hits.add(index, p);
                //p.setProcessedPeakList(ppl);
                if(hits.size() > NUMSCORED) {
                    hits.remove(hits.size()-1);
                }
            }
       
        }
        //System.out.println("new peptide added: " + p.getSequence() + ": " + p.getNumPeaksMatched()); 
    }
    public int getNumPeptideHit() {
        return hits.size();
    }
    private ArrayList<ScoredPeptideHit> getTopHits(int numTopHits) {
        /*
        if (!isSorted) {
            Collections.sort(hits);
            isSorted = true;
        }
        */
        // group PeptideHits by sequence string
        HashMap<String, ScoredPeptideHit> topHits = 
                    new HashMap<String, ScoredPeptideHit>(numTopHits);
        ArrayList<ScoredPeptideHit> topHitList = new ArrayList<ScoredPeptideHit>(numTopHits);
        for (int i = 0; i < hits.size(); i++) {
            PeptideHit ph = hits.get(i);
            //String seq = ph.getExtendedSequence();
            //String seq = ph.getSequence();
            String seq = ph.getExactSequence();
//if(seq.indexOf("DDIAALVVD") != -1)
//System.out.println(seq + "\t" + hits.size());
            ScoredPeptideHit s = topHits.get(seq);
            if(s == null) {
                s = new ScoredPeptideHit(seq, primaryScoreType, secondaryScoreType);
                topHitList.add(s);
                topHits.put(seq, s);
            }
            s.addPeptideHit(ph);
            if(topHitList.size() == numTopHits) {
                break;
            }
        }
        //return topHits.values(); 
        return topHitList; 
    }

    // using binary search
    private int findIndex(double p, int low, int high) {
        int lastIndex = hits.size() - 1;
        if(lastIndex == -1 || p <= hits.get(0).getProbability()) {
            return 0;
        } else if(high == low || p >= hits.get(lastIndex).getProbability()) {
            return lastIndex + 1;
        }
        
        int mid = 0;
        while (high > low) {
            //mid = low + (high - low) / 2;
            mid = (high + low) / 2;
            double midP = hits.get(mid).getProbability();
            if (p < midP)
                high = mid;
            else if (p > midP)
                low = mid;
            else
                return mid;
            if((high-low) == 1) {
               if(p <= hits.get(low).getProbability()) {
                   return low;
               }
               if(p > hits.get(high).getProbability()) {
                   return high + 1;
               } else {
                   return high;
               }
           }
        }
        return mid;
    }
 
}
