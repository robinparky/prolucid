/**
 * @file InclusionList.java
 * This is the source file for edu.scripps.pms.util.spectrum.InclusionList
 * @author Tao Xu
 * @date $Date
 */



import edu.scripps.pms.util.io.*;

//import edu.scripps.pms.util.io.DTASelectFilterReader; 
import edu.scripps.pms.util.dtaselect.Peptide; 
import edu.scripps.pms.util.dtaselect.Protein; 
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.PmsUtil; 
import edu.scripps.pms.util.stat.StatCalc; 
import java.io.*;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import java.util.HashSet;
import edu.scripps.pms.util.spectrum.*;

// command line processor
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

// for generating the inclusion list and exclusion list 
// -s segment length in minutes
public class InclusionList {
    public static final String USAGE = "\n\n!!! USAGE: InclusionList -i input -r ret_tolerance -s segment -p ppm_tolerance -f frequence_threshold -d use_dtaselect_result !!!";
    public static final String FOLDERNAME = "inclusionlists/";
    private static ArrayList<HashSet<String>> peptideSets = new ArrayList<HashSet<String>>();
//    private static HashSet<String> peptidesInCommon = new HashSet<String>();
    private static String inputFile = null;
    // segment length in minutes
    private static double segmentlength = 10;
    private static HashMap<String, Peptide []> filename2identifiedscans = new HashMap<String, Peptide []>(1000);   
    private static PrintStream allidentified; 
    private static double overlap = 2.5; // retention time tolerance
    private static double minummz = 400;
    private static double maxmz = 2000;
    public static double accuracyfactor = 200;
    private static boolean identifiedmzonly = false; 
    private static boolean useDtaselectResult = false; 
    private static int minfreq = 2; // minimum frequence - 1 
    private static int MAXNUMPEAKS = 500; //numpeaks per segment 
    public static final double MASSDIFFC12C13 = 1.003354826;
    public static final double TWOMASSDIFFC12C13 = MASSDIFFC12C13*2;
    public static final double THREEMASSDIFFC12C13 = MASSDIFFC12C13*3;

 
    public static ArrayList<Protein> getDTASelectProteins(String dtaselectFilterFile) throws IOException {

        ArrayList<Protein> result = new ArrayList(10000);
        if(useDtaselectResult) {
            DTASelectFilterReader reader = new DTASelectFilterReader(dtaselectFilterFile);
            for (Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {
                result.add(itr.next());
            }
            reader.close(); 
        }
        return result;
    }
    // need to use --extra -DM -t 0 options in DTASelect2
    public static void getIdentifiedSpectra(ArrayList<Protein> dtaselectProteins) throws IOException {
        
        for(Protein p : dtaselectProteins) {
            int peptidecount = 1;
            String acc = p.getAccession();
            if(!acc.startsWith("Re")) {
                //HashSet<String> peptideadded = new HashSet<String>(100000); 
                for(Iterator<Peptide> pepItr=p.getPeptides(); pepItr.hasNext(); ) {
              
                    Peptide peptide = pepItr.next();
                    String ms2file = peptide.getFileName();
                    int scan = new Integer(peptide.getLoScan()).intValue();
                    Peptide [] scans = filename2identifiedscans.get(ms2file);
                    if(scans == null) {
                        scans = new Peptide[100000];
                //        scans[scan] = peptide;
                        filename2identifiedscans.put(ms2file, scans);
                    }
                    scans[scan] = peptide;;

                }
            }
        }
        
    }

    public static void outputInfo(String ms2file) throws IOException {
        SpectrumReader ms2sr = new SpectrumReader(ms2file+".ms2", "ms2");
        SpectrumReader ms1sr = new SpectrumReader(ms2file+".ms1", "ms1");
        Peptide [] identifiedscans = filename2identifiedscans.get(ms2file);
        if(identifiedscans == null) {
            return;
        }
        //PrintStream identified = new PrintStream(new File("identified/" + ms2file+"_identified.ms2"));
        PrintStream identified = allidentified;
        PrintStream unidentified = new PrintStream(new File("identified/" + ms2file+"_unidentified.ms2"));
        ArrayList<PeakList> ms2spectra = ms2sr.getSpectraList();
        ArrayList<PeakList> ms1spectra = ms1sr.getSpectraList();
        ms2sr.closeDataFile();
        ms1sr.closeDataFile();

        for(Iterator<PeakList> it = ms1spectra.iterator(); it.hasNext();) {
            PeakList spectrum = it.next();
            int scan = spectrum.getLoscan();
            double ret = spectrum.getRetentionTime();            
            double precmz = spectrum.getPrecursorMass();            
            if(identifiedscans != null && identifiedscans[scan] != null) {
                identified.println(precmz + "\t" + spectrum.getLoscan() + "\t" + ret + "\t" + spectrum.getPrecursorScanNumber() + "\t" + spectrum.getIonInjectionTime() + "\t" + spectrum.getPrecursorInt());
            } else {
                unidentified.println(precmz + "\t" + spectrum.getLoscan() + "\t" + ret + "\t" + spectrum.getPrecursorScanNumber() + "\t" + spectrum.getIonInjectionTime() + "\t" + spectrum.getPrecursorInt());

            }
        }
        //identified.close();
        unidentified.close();
    }
    private static String filename2SaltStep(String spectrumFile) {
        // System.out.println("spectrum file " + spectrumFile);
        String [] arr = spectrumFile.split("\\_");
        String saltstep = arr[arr.length-1];
        return saltstep;
    }
    public static void outputMzLists(String ms2file) throws IOException {
        System.out.println("ms2 file: " + ms2file);
        String saltstep = filename2SaltStep(ms2file); 
        SpectrumReader ms2sr = new SpectrumReader(ms2file+".ms2", "ms2");
        SpectrumReader ms1sr = new SpectrumReader(ms2file+".ms1", "ms1");
        Peptide [] identifiedscans = filename2identifiedscans.get(ms2file);
        //if(identifiedscans == null) {
        //    return;
        //}
        //PrintStream identified = new PrintStream(new File("identified/" + ms2file+"_identified.ms2"));
        PrintStream identified = allidentified;
        System.out.println("salt step " + saltstep);
        //PrintStream unidentified = new PrintStream(new File("identified/" + ms2file+"_unidentified.ms2"));
        System.out.println("start reading ms2: " + ms2file);
        ArrayList<PeakList> ms2spectra = ms2sr.getSpectraList();
        System.out.println("finished reading ms2: " + ms2file);
        System.out.println("start reading ms1: " + ms2file);
        ArrayList<PeakList> ms1spectra = ms1sr.getSpectraList();
        System.out.println("finished reading ms1: " + ms2file);
        ms2sr.closeDataFile();
        ms1sr.closeDataFile();

        double maxretentiontime = ms1spectra.get(ms1spectra.size()-1).getRetentionTime(); 
        int numsegments = (int)(maxretentiontime/segmentlength) + 1;

        // process the first step
        for(int i = 0; i < numsegments; i++) {
            System.out.println("start " + i + " segment");
            double starttime = i*segmentlength - overlap; 
            double endtime = (i+1)*segmentlength;
            ArrayList<Double> identifiedmzs = new ArrayList<Double>(5000); // identified peaks in ms2
            ArrayList<Double> unidentifiedmzs = new ArrayList<Double>(5000); // unidentified peaks in ms2
            for(Iterator<PeakList> it = ms2spectra.iterator(); it.hasNext();) {
                PeakList spectrum = it.next();
                double ret = spectrum.getRetentionTime();            
                if(ret < starttime) {
                    continue;
                }
                if(ret > endtime) {
                    break;
                }

                int scan = spectrum.getLoscan();
                int chargestate = spectrum.getFirstChargeState();
                double precmz = spectrum.getPrecursorMass();            
                if(identifiedscans != null && identifiedscans[scan] != null) {
                    identifiedmzs.add(new Double(precmz));
                    identifiedmzs.add(new Double(precmz + MASSDIFFC12C13/chargestate)); // for M+1 peak
                    identifiedmzs.add(new Double(precmz + TWOMASSDIFFC12C13/chargestate)); // for M+2 peak
                    identifiedmzs.add(new Double(precmz + THREEMASSDIFFC12C13/chargestate)); // for M+3 peak
                } else {
                    unidentifiedmzs.add(new Double(precmz));
                }

            }
            System.out.println("In " + saltstep + "_" + (i+1));
            ArrayList<Peak> exclusionlist = new ArrayList<Peak>(1000);
            ArrayList<Peak> unidedpeaks = getUnidentifiedPeaks(starttime, endtime, ms1spectra.iterator(), 
                           identifiedmzs.iterator(), unidentifiedmzs.iterator(), exclusionlist);
            System.out.println("Number of unidentified peaks in " + saltstep + "_" + (i+1) + " is " + unidedpeaks.size());
            System.out.println("Number of peaks to be excluded in " + saltstep + "_" + (i+1) + " is " + exclusionlist.size());
            String segname = "";
            int iplus1 = i + 1;
            segname += iplus1 > 9? ("" + iplus1) : ("0" + iplus1);
            int startrt = (int)(i*segmentlength*10);
            int endrt = (int)(iplus1*segmentlength*10);
            //String rtrange = "-rt" + (startrt/10.0) + "-" + (endrt/10.0);
            String rtrange = "-rt" + ((int)(segmentlength));
            //rtrange = rtrange.replaceAll("\\.", "_");
            String segfilename = saltstep + "_" + segname + rtrange  + ".csv";

            PrintStream inclusionfile = new PrintStream(new File(FOLDERNAME + "inclst_" + segfilename));
            Collections.sort(unidedpeaks, new PeakComparator(true));
            Collections.reverse(unidedpeaks);
            for(int ii = 0; ii < MAXNUMPEAKS && ii < unidedpeaks.size(); ii++) {
                Peak p = unidedpeaks.get(ii);
                inclusionfile.println(p.getM2z() + ",");
                //inclusionfile.print(p.getM2z() + "\t" + p.getIntensity() + "\n");
            }
            //inclusionfile.println();
            inclusionfile.close(); 
            
             // now for the exlusion list
            PrintStream exclusionfile = new PrintStream(new File(FOLDERNAME + "exlst_" + segfilename));
            Collections.sort(exclusionlist, new PeakComparator(true));
            Collections.reverse(exclusionlist);
            for(int ii = 0; ii < MAXNUMPEAKS && ii < exclusionlist.size(); ii++) {
                Peak p = exclusionlist.get(ii);
                exclusionfile.println(p.getM2z() + ",");
                //inclusionfile.print(p.getM2z() + "\t" + p.getIntensity() + "\n");
            }
            //exclusionfile.println();
            exclusionfile.close(); 
        }

    }
    
    private static ArrayList<Peak> getUnidentifiedPeaks(double starttime, 
            double endtime, Iterator<PeakList> peaks, 
            Iterator<Double> identifiedmzs, Iterator<Double> unidentifiedmzs, ArrayList<Peak> exclusionlist) {
         
        
        ArrayList<Peak> ms1peaks = new ArrayList<Peak>(10000);
        int timed = (int)(endtime - starttime)*60;
        int mzd = (int)((maxmz+10)*accuracyfactor); // max mz * accuraacy factor
        int timeoffset = (int)starttime*60;
       
        boolean [] identified = new boolean[mzd];
        while(identifiedmzs.hasNext()) {
           int idedindex = (int)(identifiedmzs.next().doubleValue()*accuracyfactor);
           identified[idedindex] = true;
           
           //identified[idedindex+1] = true;
        }
        if(!identifiedmzonly) {
            while(unidentifiedmzs.hasNext()) {
                int idedindex = (int)(unidentifiedmzs.next().doubleValue()*accuracyfactor);
                identified[idedindex] = true;
                //identified[idedindex+1] = true;
            }
        }
 
        int [] maxFreq = new int[mzd]; 
        int [] currentFreq = new int[mzd]; 
        boolean [] fps = new boolean[mzd]; // found in prevous scan or not
        //boolean [][] mzs = new boolean[timed][mzd];
        while(peaks.hasNext()) {
            PeakList ps = peaks.next();
            double retentiontime = ps.getRetentionTime();
            if(retentiontime < starttime) {
                continue;
            }
            if(retentiontime > endtime) {
                break;
            }
            ListIterator<Peak> pit = ps.getPeaks();
            boolean [] fcs = new boolean[mzd]; // found in current scan 
            int [] freq = new int[mzd]; 
            while(pit.hasNext()) {
                Peak p = pit.next();
                int mzindex = (int)(p.getM2z()*accuracyfactor);
                fcs[mzindex] = true;
                freq[mzindex] = currentFreq[mzindex] + 1;
                if(fps[mzindex]) {
                    currentFreq[mzindex]++;
                    if(freq[mzindex] > maxFreq[mzindex]) {
                         maxFreq[mzindex] = freq[mzindex];
                    }
                }
            }
            fps = fcs;
            currentFreq = freq;
        }
        //System.out.println("massFreq length: " + maxFreq.length);
        // remove low frequence peaks
        for(int i = 0; i < maxFreq.length; i++) {
            if(maxFreq[i] <= minfreq) {
                maxFreq[i] = 0;
            } else {
                //System.out.println(i + "\t" + maxFreq[i]);
            }
        }
        return getUnidentifiedPeaks(maxFreq, identified, exclusionlist);
    }
    private static void removeIsotopicPeak(int [] freq, int mzindex, int charge, int numisotope) {
        for(int i = 1; i <= numisotope; i++) {
            int index = mzindex + (int)(i*MASSDIFFC12C13*accuracyfactor/charge); 
            freq[index] = 0;
            freq[index+1] = 0;
            freq[index-1] = 0;
        }
    }
    // remove the isotopic peaks and return the charge state of the peak
    // return -1 if the charge state cannot be determined based on isotopic peaks
    private static int removeIsotopicPeaks(int [] freq, int index) {
        int charge = -1;
        int charge2index = index + (int)(MASSDIFFC12C13*accuracyfactor/2); 
        int charge3index = index + (int)(MASSDIFFC12C13*accuracyfactor/3); 
        int charge1index = index + (int)(MASSDIFFC12C13*accuracyfactor); 
//System.out.println("index: " + index + "\tcharge1 index: " + charge1index + "\tcharge2index: " + charge2index + "\tcharge3index: " + charge3index);
        if(freq[charge2index] > 0 || freq[charge2index+1] > 0 || freq[charge2index-1] > 0) {
            charge = 2;
        } else if(freq[charge3index] > 0 || freq[charge3index+1] > 0|| freq[charge3index-1] > 0) {
            charge = 3;

        } else if(freq[charge1index] > 0 || freq[charge1index+1] > 0 || freq[charge1index-1] > 0) {
            charge = 1;
        }
      
        switch(charge) {  
            case 1: 
                 removeIsotopicPeak(freq, index, 1, 3);
                 break; 
            case 2:
                 removeIsotopicPeak(freq, index, 2, 3);
                 break; 
            case 3:
                 removeIsotopicPeak(freq, index, 3, 3);
                 break;
        }
         
        return charge;
    }
    private static ArrayList<Peak> getUnidentifiedPeaks(int [] maxfreq, 
                  boolean [] identified, ArrayList<Peak> exclusionlist) {
        ArrayList<Peak> peaks = new ArrayList<Peak>(10000);
        int maxindex = (int)(maxmz*accuracyfactor);
        for(int i = (int)(minummz*accuracyfactor); i < maxindex; i++) {
            if(maxfreq[i] <= minfreq) {
                continue;
            }
            if(identified[i] || identified[i+1] || identified[i-1]) {
                //System.out.println("identified\t" + i/accuracyfactor + "\t" + maxfreq[i]);
                exclusionlist.add(new Peak(i/accuracyfactor, maxfreq[i]));
                
                maxfreq[i] = 0;
                maxfreq[i+1] = 0;
            } else {
                int charge = removeIsotopicPeaks(maxfreq, i);
                if(charge != -1) { 
                    //System.out.println("unidentified\t" + i/accuracyfactor + "\t" + maxfreq[i]);
                    peaks.add(new Peak(i/accuracyfactor, maxfreq[i]));
                    maxfreq[i] = 0;
                    maxfreq[i+1] = 0;
                } else {
                    //System.out.println("noisotope\t" + i/accuracyfactor + "\t" + maxfreq[i]);
                    maxfreq[i] = 0;
                }
            }
        } 
        return peaks;
    }
    public static ArrayList<String> getMs2Files(String dir) {

         ArrayList<String> ms2Files = new ArrayList<String>();
         File currentDir = new File(dir);
         for(String s : currentDir.list()) {
             if(s.endsWith(".ms2")) {
                 //System.out.println("got ms2 file: " + s);
                 ms2Files.add(s.split(".ms2")[0]);
             }
         }
         return ms2Files;
    }


    public static void main(String args[]) throws Exception {

        // Command line processing
        Options opts = new Options();
        Option inputOpt = new Option
            ("i", "input", true, "Input file name - DTASelect-filer.txt run w/o -DM option");
        Option outOpt = new Option
            ("r", "ret_tolerance", true,  "retntion time tolerance, i.e., overlap between segments");
        Option segmentOpt = new Option
            ("s", "segment", true,  "time segment length in mimutes");
        Option ppmOpt = new Option
            ("p", "ppm", true,  "mass accuracy in ppm");
        Option freqOpt = new Option
            ("f", "freq", true,  "frequence threshold");
        Option dOpt = new Option
            ("d", "dtaselect", false,  "use dtaselect result");
        //inputOpt.setRequired(true);
        //pValueOpt.setRequired(true);
        //hprdOpt.setRequired(false);
        //opts.addOption(dbNameOpt);
        opts.addOption(inputOpt);
        opts.addOption(outOpt);
        opts.addOption(segmentOpt);
        opts.addOption(ppmOpt);
        opts.addOption(freqOpt);
        opts.addOption(dOpt);

        BasicParser cliParser = new BasicParser();
        CommandLine cli = null;

        try {
            cli = cliParser.parse(opts, args);
            
        } catch (Exception e) {
            HelpFormatter helpFormatter = new HelpFormatter();
            helpFormatter.printHelp(USAGE, opts);
            e.printStackTrace();
            System.out.println
                ("\nPlease use \"java -Xmx1500M ... \" when processing a large dataset.");
            System.exit(1);
        }
        inputFile = cli.getOptionValue("i");
        if(inputFile == null) {
            inputFile = "DTASelect-filter.txt";
        }
        String seglength = cli.getOptionValue("s");
        if(seglength != null) {
            segmentlength = Double.parseDouble(seglength);
        }
        String overlapopt = cli.getOptionValue("r");
        if(overlapopt != null) {
            overlap = Double.parseDouble(overlapopt);
        }
        
        String ppm = cli.getOptionValue("p");
        if(ppm != null) {
            accuracyfactor = (int)(1000/Double.parseDouble(ppm));
        }

        String freq = cli.getOptionValue("f");
        if(freq != null) {
            minfreq = (int)(Double.parseDouble(freq)) - 1;
        }
        String usedtaselect = cli.getOptionValue("d");
        if(cli.hasOption('d')) {
            useDtaselectResult = true;
        } 

        try {
            //System.out.println("p value cutoff: " + minPValue);
            getIdentifiedSpectra(getDTASelectProteins(inputFile));
            ArrayList<String> ms2files = getMs2Files(".");
            //Runtime.getRuntime().exec("mkdir identified");
            //allidentified = new PrintStream(new File("identified/ided.txt"));
            new File(FOLDERNAME).mkdir(); 
            //allidentified = new PrintStream(new File("ided.txt"));
            for(Iterator<String> it = ms2files.iterator(); it.hasNext();) {
                String ms2file = it.next();
                System.out.print("Processing " + ms2file + ".ms2 file ...");
                //outputInfo(ms2file);
                outputMzLists(ms2file);
                System.out.println();
            }

        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
        }
    }


}
