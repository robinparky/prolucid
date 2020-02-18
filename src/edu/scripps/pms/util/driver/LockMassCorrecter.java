/**
 * @file LockMassCorrecter.java
 * This is the source file for edu.scripps.pms.util.spectrum.LockMassCorrecter
 * @author Tao Xu
 * @date $Date: 2008/10/24 22:39:14 $
 */



import edu.scripps.pms.util.io.*;

import edu.scripps.pms.mspid.MassSpecConstants;
import edu.scripps.pms.util.PmsUtil; 
import edu.scripps.pms.util.io.FileUtils; 
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

// adjust the ms1 peak intensity according to maximum ion injection time
public class LockMassCorrecter {
    public static final String USAGE = "\n\n!!! USAGE: lockmasscorrecter -i ms1file -p ppmvalue -m 1or2 -w movingavgwindowsize -t maxioninjectiontime!!!";
    private static ArrayList<HashSet<String>> peptideSets = new ArrayList<HashSet<String>>();
//    private static HashSet<String> peptidesInCommon = new HashSet<String>();
    public static int AVGMODE = 1;
    public static int MOVINGAVGMODE = 2;
    private static double [] lockmasses = null;
    private static double [] tolerances = null;
    private static double ppmcutoff = 20; 
    private static int avgwindow = 50;
    private static int mode = 1; // 1 for average, 2 for moving average
    private static double maxinjectiontime = 1; //do not adjust the intensity 

    public static double [] getCorrectionFactors(String ms1file) throws IOException {
        
        double [] corrfactors = new double[100000];
        int [] numlockfound = new int[100000];
        double defaultValue = -100000;
        double [][] deltamasses = new double[lockmasses.length][100000]; 
        for(int i = 0; i < lockmasses.length; i++) {
            for(int j = 0; j < 100000; j++) {
                deltamasses[i][j] = defaultValue;
            }
        }
        int numScans = 0; 
        SpectrumReader sr = new SpectrumReader(ms1file+".ms1", "ms1");
        for(Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList spectrum = it.next();
            int scan = spectrum.getLoscan();
            numScans = scan;
//System.out.print(scan + "\t");
            for(int i = 0; i < lockmasses.length; i++) {
                Peak p = spectrum.getMinimumDeltaMassPeak(lockmasses[i], tolerances[i]);
                if(p != null) {
                    double deltam = (p.getM2z()-lockmasses[i])/lockmasses[i]*1000000;
//System.out.print(lockmasses[i] + "\t" + p.getM2z() + "\t" + deltam + "\t");
                    deltamasses[i][scan] = deltam;
                }
            }
//System.out.println();
        }    
        numScans += 20; // for MS/MS scans
        double previouscorrfactor = 0;
        int numNoLockMass = 0;
        for(int i = 0; i < numScans; i++) {
            double sum = 0;
            int count = 0;
            for(int j = 0; j < lockmasses.length; j++) {
                double value = deltamasses[j][i];

                if(value != defaultValue) {
                    sum += value;
                    count++;
                }
                if(count > 0) {
                    corrfactors[i] = sum/count;
                    numlockfound[i] = count;
                    previouscorrfactor = corrfactors[i];
                } else {
                    corrfactors[i] = previouscorrfactor;
                    numNoLockMass++;
                }
            }
            //System.out.println(i + "\t" + corrfactors[i] + "\t" + numlockfound[i]);
        }
        double [] tempcorrfactors = new double[100000];
        if(mode == MOVINGAVGMODE) {
            for(int i = 0; i < numScans; i++) {
                tempcorrfactors [i] = getMovingAvg(i, corrfactors, numScans-1);
            }
            corrfactors = tempcorrfactors;
        }
        sr.closeDataFile(); 
        return corrfactors;
    }
   
    private static double getMovingAvg(int index, double [] values, int maxindex) {
        int start = index - avgwindow/2;
        start = start < 0? 0 : start;
        int end = index + avgwindow/2;
        end = end > maxindex? maxindex : end;
        double total = 0;
        for(int i = start; i <= end; i++) {
            total += values[i];
        }
        
        return total/(end-start);
        
    }
    private static void correctMs1(String ms1file, double [] corrfactors) throws IOException {
        PrintStream lockedms1 = new PrintStream(new File("locked/" + ms1file+"_locked.ms1"));
        SpectrumReader sr = new SpectrumReader(ms1file+".ms1", "ms1");
        boolean isFirstSpectrum = true;

        System.out.println("maximum ion injection time: " + maxinjectiontime);
        for(Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList spectrum = it.next();
            int scan = spectrum.getLoscan();
            double ppm = corrfactors[scan];

            double ioninjectiontime = spectrum.getIonInjectionTime();
            double intensityfactor = maxinjectiontime == 1 ? 1 : maxinjectiontime/ioninjectiontime;
            intensityfactor = intensityfactor < 1? 1 : intensityfactor;

            System.out.println("ioninjectiontime: " + ioninjectiontime + "\tintensityadjustment: " + intensityfactor);
            for(Iterator<Peak> itp = spectrum.getPeaks(); itp.hasNext();) {
               Peak p = itp.next();
               double mz = p.getM2z();
               double correctedmz = mz - mz*ppm/1000000;
               p.setM2z(correctedmz);
               p.setIntensity(p.getIntensity()*intensityfactor);
            }
            if(isFirstSpectrum) {
                lockedms1.print(spectrum.getSpectrumWithHlines());
                isFirstSpectrum = false;
            } else {
                lockedms1.print(spectrum.getSpectrumWithoutHlines());
            }
        }
        lockedms1.close();
    } 
    public static void outputSpectra(String ms1file) throws IOException {
        String ms2file = ms1file + ".ms2";
        String ms2output = null;
        String ms2input = null;
        if(ms1file.contains("ftms") && !new File(ms2file).exists()) {
            ms2input = ms2file.replaceAll("ftms", "itms");
            ms2output = ms1file.replaceAll("ftms", "itms");
            
        } else {
            ms2input = ms2file;
            ms2output = ms1file;
        } 
        
        System.out.println("ms2input: " + ms2input + "\t ms2output: " + ms2output);
        double [] corrfactors = getCorrectionFactors(ms1file);
        SpectrumReader sr = new SpectrumReader(ms2input, "ms2");
        PrintStream lockedms2 = new PrintStream(new File("locked/" + ms2output +"_locked.ms2"));
        //PrintStream lockedms1 = new PrintStream(new File("locked/" + ms1file+"_locked.ms1"));
        boolean isFirstSpectrum = true;
        
        for(Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList spectrum = it.next();
            
            int scan = spectrum.getLoscan();
            double ppm = corrfactors[scan];
            double precmz = spectrum.getPrecursorMass();
            double correctedmz = precmz - precmz*ppm/1000000;
//System.out.println("precmz: " + precmz + "\tcorrectedmz: " + correctedmz);
            spectrum.setPrecursorMass(correctedmz);
            spectrum.addIline("I\tMZ\t" + precmz);
            spectrum.addIline("I\tDeltaMassPPM\t" + ppm);
            ArrayList<Zline> zlines = new ArrayList<Zline>();
            for(Iterator<Zline> itz = spectrum.getZlines(); itz.hasNext();) {
                Zline z = itz.next();
                int chargestate = z.getChargeState();
                double mplush = correctedmz*chargestate - chargestate*MassSpecConstants.MASSH + MassSpecConstants.MASSH;
                zlines.add(new Zline(chargestate, mplush));
                spectrum.addIline("I\tPZ\t" + chargestate + "\t" + z.getM2z());
            }
            spectrum.resetZlines(zlines);
            if(isFirstSpectrum) {
                lockedms2.print(spectrum.getSpectrumWithHlines());
                isFirstSpectrum = false;
            } else {
                lockedms2.print(spectrum.getSpectrumWithoutHlines());
            }
        }    
        lockedms2.close();
        correctMs1(ms1file, corrfactors);
        //lockedms1.close();
    }

    private static void readParameters(String file) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(file));
        HashSet<Double> mzs = new HashSet<Double>();
        ArrayList<Double> locks = new ArrayList<Double>();
        String line = null;
        while((line = br.readLine()) != null) {
            line = line.trim();
            if(line.length() > 0) {
                Double d = Double.parseDouble(line);
                if(!mzs.contains(d)) {
                    mzs.add(d);
                    locks.add(d);
                }
//System.out.println(line);
            }
        }
        lockmasses = new double[locks.size()];
        tolerances = new double[locks.size()];
        int count = 0;
        for(Iterator<Double> it = locks.iterator(); it.hasNext();) {
            double mz = it.next().doubleValue();
            lockmasses[count] = mz;
            tolerances[count] = mz*ppmcutoff/1000000;
            count++;
        }
        
        br.close();
    }
    public static void main(String args[]) throws Exception {
   

        // Command line processing
        Options opts = new Options();
        Option inputOpt = new Option
            ("i", "input", true, "Input file name - ms1filename");
        Option outOpt = new Option
            ("p", "ppmcutoff", true,  "ppm cut off for lock mass");
        Option modeOpt = new Option
            ("m", "mode", true,  "correction mode, 1 for avg, 2 for moving avg");
        Option windowsizeOpt = new Option
            ("w", "windowsize", true,  "moving average window size");
        Option maxinjtOpt = new Option
            ("t", "maxinjt", true,  "maximum ion injection time");
        //inputOpt.setRequired(true);
        //pValueOpt.setRequired(true);
        //hprdOpt.setRequired(false);
        //opts.addOption(dbNameOpt);
        opts.addOption(inputOpt);
        opts.addOption(outOpt);
        opts.addOption(modeOpt);
        opts.addOption(windowsizeOpt);
        opts.addOption(maxinjtOpt);

        BasicParser cliParser = new BasicParser();
        CommandLine cli = null;

        System.out.println(USAGE);
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

        ArrayList<String> ms1files = new ArrayList<String>(20);

        String inputFile = cli.getOptionValue("i");
        if(inputFile != null) {
            if(!inputFile.endsWith(".ms1")) { // not a ms1 file
                System.out.println(USAGE);
                System.err.println("Input file is not a ms1 file");
                System.exit(1);
                
            }
            String [] arr = inputFile.split(".ms1");
            ms1files.add(arr[0]);
        } else {
            ms1files = FileUtils.getFilePrefixes(".", ".ms1");
        }
        String ppm = cli.getOptionValue("p");
        if(ppm != null) {
            ppmcutoff = Double.parseDouble(ppm);
        }
        String corrmode = cli.getOptionValue("m");
        if(corrmode != null) {
            mode = Integer.parseInt(corrmode);
        }
        String windowsize = cli.getOptionValue("w");
        if( windowsize != null) {
            avgwindow = Integer.parseInt(windowsize);
        }
        String injt = cli.getOptionValue("t");
        if( injt != null) {
            maxinjectiontime = Double.parseDouble(injt);
        }
        readParameters("lockmass.params");
        try {
            //System.out.println("p value cutoff: " + minPValue);
            Runtime.getRuntime().exec("mkdir locked");
            for(Iterator<String> it = ms1files.iterator(); it.hasNext();) {
                String ms1file = it.next();
                System.out.print("Processing " + ms1file + ".ms1 and " +  ms1file + ".ms2 ...");
                outputSpectra(ms1file);
                System.out.println();
            }

        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
        }
    }
}
