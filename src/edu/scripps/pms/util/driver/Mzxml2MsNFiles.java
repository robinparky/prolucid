/**
 * @file Mzxml2MsNFiles.java
 * This is the source file for edu.scripps.pms.util.io.Mzxml2MsNFiles
 * @author Tao Xu
 * @date $Date
 */
package edu.scripps.pms.util.driver;


import java.io.*;
import java.util.List;
import java.util.Stack;
import java.util.Set;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.io.MzxmlSpectrumReader;
//import edu.scripps.pms.util.MZXmlHandler;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
//import org.jdom.Namespace;
import org.jdom.input.SAXBuilder;


public class Mzxml2MsNFiles {

    public static final String USAGE = "USAGE: mzxml2msn foldername/mzxmlfilename mzxmlfileextension";
    public static final String version = "1.0.3";
    public static final String releaseDate = "08/08/2011";
    public static final int MAXMSLEVEL = 10;
    private static ArrayList<String> getFiles(String dir, String extension) throws IOException {
        
         ArrayList<String> files = new ArrayList<String>();
         File currentDir = new File(dir);
         if(currentDir.isFile()) {
             String [] arr = dir.split("\\." + extension);
             files.add(arr[0]); //dir is a mzxml file
             return files;
         } 
         for(String s : currentDir.list()) {
             if(s.endsWith(extension)) {
                 String [] arr = s.split("\\." + extension);
                 //files.add(dir + "/" + arr[0]);
                 files.add(dir + File.separator  + arr[0]);
             }
         }
         return files;    
    }

    private static void outputHLines(PrintStream ps, int mslevel, String mzxmlfile, String mzxmlextension) {
        StringBuffer sb = new StringBuffer(1000);
        sb.append("H\tExtractor\tMzxml2MsFiles\n");
        sb.append("H\tExtractorVersion\t" + version + "\treleased on " + releaseDate + "\n");
        sb.append("H\tComments\tMzxml2MsFiles written by Tao Xu in the Yates lab at The Scripps Research Institute\n");
        sb.append("H\tComments\tPlease send comments or bug reports to taoxu@scripps.edu\n");
        sb.append("H\tOriginalMzxmlFile\t" + mzxmlfile + "." + mzxmlextension + "\n");
        sb.append("H\tMsLevle\t" + mslevel + "\n");
        ps.print(sb.toString()); 
    }

    private static void outputSpectrum(PrintStream [] pss, int [] spectracounts, MzxmlPeakList mpl, String mzxmlfile, String mzxmlextension) throws Exception {
        //if(mpl == null) return;
        if(mpl == null || mpl.numPeaks() < 1) return; // ignore the spectrum that does not have any peaks
        int mslevel = mpl.getMsLevel();
//System.out.println("mslevel: " + mslevel);
        if(mslevel < pss.length) {
            PrintStream ps = pss[mslevel];
            if(ps == null) {
                 ps = new PrintStream(mzxmlfile + ".ms" + mslevel);
                 pss[mslevel] = ps;
                 outputHLines(ps, mslevel, mzxmlfile, mzxmlextension);
            }
             
//System.out.println(ps);
            ps.print(mpl.getSpectrumWithoutHlines()); 
            spectracounts[mslevel]++;

        }
    } 
    // for testing 
    public static void main(String args[]) throws Exception {
        System.out.println(USAGE);
        String fileextension = args[1];
        if(fileextension.startsWith(".")) fileextension = fileextension.substring(1);
        String folder = args[0]; // the folder can be a file name
        ArrayList<String> files = getFiles(folder, fileextension);
        System.out.println("Preparing to convert mzxml file/s in " + folder + " to ms files ...");


        for(String file : files) {
            String mzxmlfile = file + "." + fileextension;
            String ms2file = file + ".ms2";
            MzxmlSpectrumReader msr = new MzxmlSpectrumReader(mzxmlfile);
            System.out.println("Converting " + mzxmlfile + " to " + ms2file + " now. It may take a while, so please be patient ...");
            PrintStream [] pss = new PrintStream[MAXMSLEVEL];
            int [] spectracounts = new int[MAXMSLEVEL];
            //System.out.println("Number of scans: " + scan2Position.size());
            int mslevel = 0;
            Iterator<MzxmlPeakList> it = msr.getSpectra(mslevel);
            //for(it = msr.getSpectraWithChildren(); it.hasNext();) {
            while(it.hasNext()) {
                //msr.printMzxmlPeakList(it.next());
                MzxmlPeakList mpl = it.next();
                //ps.print(mpl.getSpectrumWithoutHlines()); 
                outputSpectrum(pss, spectracounts, mpl, file, fileextension); 
                //System.out.println("\nscanNumber: " + mpl.getLoscan() + "\tmsLevel: " + mpl.getMsLevel() + "\tnumChildren: " + mpl.getNumChildSpectrum());
                //msr.printMzxmlPeakList(it.next());      
        
            }
            for(int i = 0; i < MAXMSLEVEL; i++) {         
                PrintStream ps = pss[i];
                if(ps != null) ps.close();
                if(spectracounts[i] != 0) System.out.println("Number of ms" + i + " spectra extracted from " + file + "." + fileextension + ":\t" + spectracounts[i]);
       
            }
        }
    }

}
