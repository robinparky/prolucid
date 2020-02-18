/**
 * @file Signal2NoiseRatioCalculator.java
 * This is the source file for edu.scripps.pms.util.spectrum.Signal2NoiseRatioCalculator
 * @author Tao Xu
 * @date $Date: 2012/09/28 00:04:16 $
 */



import edu.scripps.pms.util.io.*;

//import edu.scripps.pms.util.io.DTASelectFilterReader; 
import edu.scripps.pms.util.dtaselect.Peptide; 
import edu.scripps.pms.util.dtaselect.Protein; 
import java.io.*;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import edu.scripps.pms.util.spectrum.*;

public class Signal2NoiseRatioCalculator {
    public static final String USAGE = "java Signal2NoiseRatioCalculator ms1File ms2File percentage(20/30/40) intensityThreshold\n" +
                                       "!!! ms1, ms2 and DTASelect-filter.txt files must be in the same folder!!!";
    public static final int MILLION = 1000000;
    public static final float MASSDIFFC12C13 = 1.003354826f;
    private String experimentName;
    private String ms1FileName;
    private String ms2FileName;
    private String sqtFileName;
    private int percentage;
    private PeakList [] ms1Spectra = new PeakList[1000000];
    private int numMs2Scans = 0;
    private int [] scanNumbers = new int[1000000]; 
    private double [] signal2NoiseRatio = new double[1000000]; 
    private double [] zScores = new double[1000000]; 
    private double [] intensities = new double[1000000]; 
    private double [] means = new double[1000000]; 
    private double [] sigmas = new double[1000000]; 
    private HashMap<Integer, Double> scan2DeltaMass = new HashMap<Integer, Double>();
//    private double [] xcorrs = new double[1000000];
//    private double [] deltaMasses = new double[1000000];
    private int intensityThreshold; // for intensity threshold
    public Signal2NoiseRatioCalculator(String ms1, String ms2, int percentage, int threshold)throws IOException {
        ms1FileName = ms1; 
        ms2FileName = ms2; 
        sqtFileName = ms2.substring(0, ms2.indexOf(".")) + ".sqt";
        //System.out.println("sqt file: " + sqtFileName);
        this.percentage = percentage; 
        intensityThreshold = threshold;
        System.out.println("Start reading DTASelect-filter.txt file");
        calcDeltaMasses();
        System.out.println("Finished reading DTASelect-filter.txt file");
        readMs1Spectra();
        System.out.println("Finished reading ms1 file");
     
        calcSignal2NoiseRatio();
        System.out.println("Finished calculating signal2NoiseRatio");
        output();
        
    }
    
    public void output() throws IOException {
         
        PrintWriter resultFile = new PrintWriter(new BufferedWriter(new FileWriter(percentage + "percent"+intensityThreshold+ms2FileName+".ratio")));
        PrintWriter problemFile = new PrintWriter(new BufferedWriter(new FileWriter(percentage + "percent"+intensityThreshold+ms2FileName+".problem")));
        resultFile.println("ScanNumber\tIntensity\tSignal2NoiseRatio\tzScore\tmean\tsigma\tdeltaMass(ppm)");
        for(int i = 0; i < numMs2Scans; i++) {
            Double d =  scan2DeltaMass.get(new Integer(scanNumbers[i]));

            double deltaMass = d == null? -1000 : d.doubleValue();
            String s = (scanNumbers[i] + "\t" + intensities[i] + "\t" + signal2NoiseRatio[i] + "\t" + zScores[i] + "\t" + means[i] + "\t" + sigmas[i] +"\t" + deltaMass);
            resultFile.println(s);
            if(intensities[i] < intensityThreshold) {
                problemFile.println(s);    
            }
            
        }
        resultFile.close();
        problemFile.close();
    }
    private void readMs1Spectra() throws IOException {
        SpectrumReader sr = new SpectrumReader(ms1FileName, "ms1");
        for(Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList pl = it.next();
            ms1Spectra[pl.getLoscan()] = pl;
        }
    }

    private void calcDeltaMasses() throws IOException {
        DTASelectFilterReader reader = new DTASelectFilterReader("DTASelect-filter.txt");
        for (Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {
            Protein protein = itr.next();

            for(Iterator<Peptide> pepItr=protein.getPeptides(); pepItr.hasNext(); ) {
                Peptide peptide = pepItr.next();
                scan2DeltaMass.put(new Integer(peptide.getLoScan()), new Double(getDeltaMassInPpm(peptide)));
            }
        }

    }
    private double getDeltaMassInPpm(Peptide peptide) {
        double ppm = -10000;
        double calcMass = Double.parseDouble(peptide.getCalcMHplus()); 
        double deltaMass = Double.parseDouble(peptide.getMhPlus()) - calcMass; 
        if(deltaMass < 0) {
            ppm = deltaMass/calcMass*MILLION;
        } else {
            int isotop = (int)(deltaMass+0.5);   
            deltaMass -= isotop*MASSDIFFC12C13;    
            ppm = deltaMass/calcMass*MILLION;
        }
        //System.out.println("calcMass: " + calcMass + "\tobserveredMass: " + peptide.getMhPlus() + "\tppm: " + ppm);
        return ppm;
    }
    private void calcSignal2NoiseRatio() throws IOException {

        SpectrumReader sr = new SpectrumReader(ms2FileName, "ms2");
        for(Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList pl = it.next(); 
            scanNumbers[numMs2Scans] = pl.getLoscan(); 
            int prcScan = pl.getPrecursorScan();
            double prcM2z = pl.getPrecursorMass(); // get precursor m2z
            PeakList plms1 = ms1Spectra[prcScan];
            Peak prcPeak = plms1.getPeakByM2z(prcM2z);  
            intensities[numMs2Scans] = prcPeak.getIntensity();  
            // need to get low intensity peaks and calculate the mean and variance 
            // as well as signal to noise
            double sum = 0;
            ArrayList<Peak> leastIntensePeaks = plms1.getLeastIntensePeaks(percentage); 
            int numPeaks = leastIntensePeaks.size();
            for(Iterator<Peak> ms1It = leastIntensePeaks.iterator(); ms1It.hasNext();) {
                Peak p = ms1It.next();
                sum += p.getIntensity(); 
               // System.out.println(p.getIntensity());
            }
            double mean = sum/leastIntensePeaks.size();

//System.out.print("numLeastPeaks: " + leastIntensePeaks.size() + "\tmean: " + mean);
            means[numMs2Scans] = mean;
            signal2NoiseRatio[numMs2Scans] = prcPeak.getIntensity()/mean; 
       
            double sigma = 0;
            for(Iterator<Peak> ms1It = leastIntensePeaks.iterator(); ms1It.hasNext();) {
                double diff = mean - ms1It.next().getIntensity(); 
                 
                sigma += (diff*diff);
            }
            if(numPeaks == 1) {
                sigma = mean; 
            } else {
                sigma = Math.sqrt(sigma/(numPeaks-1)); 
            } 
//System.out.print("\tsigma: " + sigma);
            sigma = sigma == 0? mean : sigma;
//System.out.println("\tsigma after: " + sigma);
            sigmas[numMs2Scans] = sigma; 
            zScores[numMs2Scans] = (intensities[numMs2Scans]-mean)/sigma; 
            numMs2Scans++;

        }
    }
    public static void main(String args[]) throws Exception {
        try {
            int threshold = Integer.parseInt(args[3]);
            int percentage = Integer.parseInt(args[2]);

            Signal2NoiseRatioCalculator snc = new Signal2NoiseRatioCalculator(args[0], args[1], percentage, threshold); 

        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("Finished");
    }
}
