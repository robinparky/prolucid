
/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2011</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Tao Xu 
 * @version 1.0
 */
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.io.*;
import java.util.*;
import java.io.*;

// this program can be used to removed unfragmented peaks from QTOF spectra in ms2 files
public class Ms2PrecursorRemover {
    private static final String USAGE = "\n\nUsage: java Ms2PrecursorRemover folder_with_ms2_files"; 
    private String ms2File;
    private String outPutFile;
    public static double massTolerance = 20.0/1000000;
    public static final double c13c12diff = 1.00033548;
    public static final int numpeaks = 7; // 2 smaller than prcmz and 5 greater than prcmz
    public static final String resultfolder = "xprec";
    

    public Ms2PrecursorRemover(String ms2File, String outfile) {
        this.ms2File = ms2File;
        this.outPutFile = outfile;
    }  

    public void outputMs2() throws IOException {

        SpectrumReader sr = new SpectrumReader(ms2File, "ms2");
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outPutFile)));
        outputHlines(sr.getHlines(), pw);
        for(Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList p = it.next(); 
            outputPeakList(pw, p);
        }
        sr.closeDataFile();
        pw.close();
    }
    //public void outputMs2() throws IOException {
        //String outfile = ms2File.substring(0, ms2File.lastIndexOf(".ms2")) + ".mgf";
        //String outfile = ms2File.substring(0, ms2File.lastIndexOf(".ms2")) + "_x_prec.ms2";
        //outputMs2(outfile);
     //   outputMs2();
    //}
    public static void main(String args[]) {
        try {
            String folder = args.length > 1? args[0] : ".";
            String xprecfolder = folder + "/" + resultfolder;
            File xprecdir = new File(xprecfolder);
            if(!xprecdir.exists()) {
                xprecdir.mkdir();
            } 
            for(Iterator<String> it = FileUtils.getFiles(folder, ".ms2").iterator(); it.hasNext();) {
                String ms2File = it.next();
                String outfile = ms2File.substring(0, ms2File.lastIndexOf(".ms2")) + "_x_prec.ms2";
                String outfileWithPath = folder + "/" + resultfolder + "/" + outfile;
                String ms2FileWithPath = folder + "/" + ms2File;
                System.out.println("I am trying remove precursor peaks in " + ms2File + ".");
                Ms2PrecursorRemover mpr = new Ms2PrecursorRemover(ms2FileWithPath, outfileWithPath);
                mpr.outputMs2();
            }
        } catch(Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
        }
     
    }

    public void outputHlines(Iterator<String> hlines, PrintWriter ps) {
        while(hlines.hasNext()) {
            ps.println(hlines.next());
        }
    }

    public void outputPeakList(PrintWriter ps, PeakList pl) throws IOException {
        pl.addIline("I\tPrecursorRemoved\tTrue");
        StringBuffer sb = new StringBuffer(5000);
        pl.getScanInfo(sb);
        
        double prcmz = pl.getPrecursorMass();
        double [] removelimits = new double[numpeaks*2]; // for 5 topic peaks
        int charge = 0;
        if(pl.getNumZlines() == 1) { // unambiguous charge state 
            for(Iterator<Zline> it = pl.getZlines(); it.hasNext();) {
                Zline z = it.next();
                charge = z.getChargeState();
                double isodiff = c13c12diff/charge;
                for(int i = -2; i < numpeaks -2; i++ ) {
                    double diff = (prcmz+(i*isodiff))*massTolerance; 
                    int index = (i+2)*2;
                    double isoshift = i*isodiff;
                    removelimits[index] = prcmz + isoshift - diff;
                    removelimits[index+1] = prcmz + isoshift + diff;
//System.out.println("index: " + index + "\tlimits i: " + i + "\t" + removelimits[index] + "\t" + removelimits[index+1] + "\tdiff: " + diff + "\tprcmz: " + prcmz + "\tcharge: " + charge);
                }
                break;
            }            
        }
        for(Iterator<Peak> peaks = pl.getPeaks(); peaks.hasNext();) {
            Peak p = peaks.next();
            if(!toBeRemoved(p, removelimits)) {
//System.out.println("not romoved\tprcmz: " + prcmz + "\tpeak mz: " + p.getM2z() + "\tintensity: " + p.getIntensity() + "\tprec charge: " + charge);
                sb.append(p.getM2z() + " " + p.getIntensity() + "\n");
            } else {
//System.out.println("romoved\tprcmz: " + prcmz + "\tpeak mz: " + p.getM2z() + "\tintensity: " + p.getIntensity() + "\tprec charge: " + charge);
            }
        }
        ps.print(sb.toString());
    }

    private boolean toBeRemoved(Peak p, double [] limits) {
        double mz = p.getM2z();
        int lastindex = limits.length - 1;
        if(mz >= limits[0] && mz <= limits[lastindex]) {
            for(int i = 0; i < lastindex-1; i+=2) {
                if(mz >= limits[i] && mz <= limits[i+1]) {
//System.out.println("true\tpeak mz: " + mz + "\ti: " + i + "\tlow limit: " + limits[i] + "\thigh lmit: " + limits[i+1]);
                    return true;
                } else {
//System.out.println("false\tpeak mz: " + mz + "\tlow limit: " + limits[i] + "\thigh lmit: " + limits[i+1]);
                }
            }
        }
        return false;
    }
}
