/**
 * @file SpectrumPeakCounter.java
 * This is the source file for edu.scripps.pms.util.spectrum.SpectrumPeakCounter
 * @author Tao Xu
 * @date $Date: 2013/02/22 17:39:51 $
 */




import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
//import java.util.ListIterator;
import java.util.Collections;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.io.SpectrumReader;

public class SpectrumPeakCounter {

    private static String USAGE = "java SpectrumPeakCounter spectrum_file_name ms1_or_ms2";
    public static void main(String args[]) throws Exception {
        
        if(args.length < 2) {
            System.out.println(USAGE);
            System.exit(0);
        }
        SpectrumReader sr = new SpectrumReader(args[0], args[1]);
        int [] freq = new int[12];
        
        Iterator<PeakList> it = sr.getSpectra();
        int counter = 0;
        int numPeaks = 0;
        //boolean sortByIntensity = true;
        while (it.hasNext()) {
            PeakList list = it.next();
            int numpeaks = list.numPeaks();
            numPeaks += numpeaks;
            counter++;
            //if (counter > 20)
            //break;
            if(numpeaks < 1) {
                freq[0]++;    
            } else if(numpeaks < 6) {
                freq[1]++;
            } else if(numpeaks < 11) {
                freq[2]++;
            } else if(numpeaks < 16) {
                freq[3]++;
            } else if(numpeaks < 21) {
                freq[4]++;
            } else if(numpeaks < 26) {
                freq[5]++;
            } else if(numpeaks < 31) {
                freq[6]++;
            } else if(numpeaks < 51) {
                freq[7]++;
            } else if(numpeaks < 101) {
                freq[8]++;
            } else if(numpeaks < 151) {
                freq[9]++;
            } else if(numpeaks < 201) {
                freq[10]++;
            } else {
                freq[11]++;
            }


        }
        System.out.println("Total number of spectra processed: " + counter);
        System.out.println("Average number of peaks per list: " + numPeaks/counter);
        sr.closeDataFile();

        System.out.println("Finished");
        sr.closeDataFile();

        System.out.println("NumPeaks<=\t0\t5\t10\t15\t20\t25\t30\t50\t100\t150\t200\t>=200"); 
        StringBuffer sb = new StringBuffer();
        sb.append("Frequency\t");
        for(int i = 0; i < freq.length; i++) {
            sb.append(freq[i]);
            sb.append("\t");
        }
        System.out.println(sb.toString());
    }

}
