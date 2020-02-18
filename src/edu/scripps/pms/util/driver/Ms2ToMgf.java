
/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2005</p>
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

public class Ms2ToMgf {
    private static final String USAGE = "\n\nUsage: java Ms2ToMgf ms2File"; 
    private String ms2File;
    public Ms2ToMgf(String ms2File) {
        this.ms2File = ms2File;
    }  
    public void outputMgf(String outFileName) throws IOException {

        SpectrumReader sr = new SpectrumReader(ms2File, "ms2");
        PrintWriter outFile = new PrintWriter(new BufferedWriter(new FileWriter(outFileName)));
        for(Iterator<PeakList> it = sr.getSpectra(); it.hasNext();) {
            PeakList p = it.next(); 
            outputPeakList(outFile, p);
        }
        sr.closeDataFile();
        outFile.close();
    }
    public void outputMgf() throws IOException {
        String outfile = ms2File.substring(0, ms2File.lastIndexOf(".ms2")) + ".mgf";
        outputMgf(outfile);
    }
    public static void main(String args[]) {
        try {
            for(Iterator<String> it = FileUtils.getFiles(".", ".ms2").iterator(); it.hasNext();) {
                String ms2File = it.next();
                System.out.println("I am trying to convert " + ms2File + " to MGF file, please be patient.");
                Ms2ToMgf m2m = new Ms2ToMgf(ms2File);
                m2m.outputMgf();
            }
        } catch(Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
        }
     
    }

    public void outputPeakList(PrintWriter ps, PeakList pl) throws IOException {
        int lowScan = pl.getLoscan();
        int highScan = pl.getHiscan();

        ps.print("BEGIN IONS\n");
        ps.print("TITLE=" + lowScan + "." + highScan + "\n");
        ps.print("SCANS=" + lowScan + "\n");
        ps.print("PEPMASS=" + pl.getPrecursorMass() + "\n");
        if(pl.getNumZlines() == 1) { // unambiguous charge state 
            for(Iterator<Zline> it = pl.getZlines(); it.hasNext();) {
                Zline z = it.next();
                ps.print("CHARGE=" + z.getChargeState() + "+\n");
                break;
            }            
        }
        for(Iterator<Peak> peaks = pl.getPeaks(); peaks.hasNext();) {
            Peak p = peaks.next();
            ps.println(p.getM2z() + " " + p.getIntensity());
        }
        ps.print("END IONS\n\n");
    }
}
