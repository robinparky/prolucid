import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Collections;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.io.*;


// add precursor signal2noise level to peptide id compare result
// all the ms1 and ms2 files should be in the current working dir
// 
public class PrecursorSignal2NoisePPM {
    private static HashMap<String, double[]> file2signaltonoise = new HashMap();

    public void precursorIntensity(String fileName)
    {
	String line = null;
	String delimiter = "\t";
	String returnStr = "";

	try{
            BufferedReader reader = 
		    new BufferedReader(new FileReader(fileName));
	    //skip to the 3rd line
	    for(int i = 0; i < 2; i++)
	    {
	        System.out.println(reader.readLine());
                
	    }

     	    //looping through the file
	    while((line = reader.readLine()) != null)
	    {
	        String[] temp = line.split(delimiter);
		StringBuffer sb = new StringBuffer(40);
	        //for(int i = 1; i < temp.length; i+=12)
	        for(int i = 1; i < temp.length; i+=13)
	        {
	    	    returnStr = getPrecursorSignal2NoisePPM(temp[i]);
		    sb.append("\t").append(returnStr);
	    	}
		System.out.println(line + "\t" + sb.toString());
	    }
	    reader.close();
        } catch(FileNotFoundException ex){
	  System.out.println("Couldn't find" +fileName +
				" please pick it.");
	} catch(Exception ex){
	  System.out.println("Error reading file " +fileName);
          ex.printStackTrace();
	} 
    }


    private double [] getNoiseLevel(String filename) {
        double [] noiselevel = new double [100000];
        try{
            SpectrumReader sr = new SpectrumReader(filename+".ms1", ".ms1");
            Iterator<PeakList> it = sr.getSpectra();
            while(it.hasNext()) {
                PeakList pl = it.next();
                int scan = pl.getLoscan();
                noiselevel[scan] = pl.getNoiseLevel(10); // use the 10 percent of peaks with least intensity

            }
            sr.closeDataFile();
        } catch(IOException e) {
            e.printStackTrace();
        }

        return noiselevel;
    }

    public String getPrecursorSignal2NoisePPM(String str)
    {
        // 110525_HEK_DigDeAPr_75Percent_125ug_100ugmL_Step08.4279.4279.2
        // if str equal "X" or "" returh "X", otherwise, 
        
        if(str == null || "X".equals(str.trim()) || "".equals(str.trim())) {
            return "X";
        }

        String [] arr = str.split("\\.");
        if(arr.length < 2) return "";
        String filename = arr[0];
        int scannum = Integer.parseInt(arr[1]);

//System.out.println("Filename: " + filename + "\tscan number: " + scannum);
        double [] signalnoises = file2signaltonoise.get(filename);
        if(signalnoises == null) {
            signalnoises = new double[100000];

            double[] noiselevels = getNoiseLevel(filename);

            file2signaltonoise.put(filename, signalnoises);
            try{
                SpectrumReader sr = new SpectrumReader(filename+".ms2", ".ms2");
                Iterator<PeakList> it = sr.getSpectra();
                while(it.hasNext()) {
                    PeakList pl = it.next();
                    int scan = pl.getLoscan();
                    double prcint = pl.getPrecursorInt();

                    int prcscannumber =  pl.getPrecursorScanNumber();

                    signalnoises[scan] = prcint/noiselevels[prcscannumber];

                }
                sr.closeDataFile();
            } catch(IOException e) {
                e.printStackTrace();
            }
        }
        
//System.out.println("return: " + signalnoises[scannum]);
        return ""+signalnoises[scannum];
        
    }

    public static void main(String[] args)
    {
	PrecursorSignal2NoisePPM reader = new PrecursorSignal2NoisePPM();
        //String idcompareresultfile = "/home/mtlin/mtlin_on_data/precursorInstensity/ident_compare_PEPTIDE1703.txt";
        if(args.length < 1) {
            System.out.println("!!!Usage: addPrecursorSignal2NoiseToPeptideIdentCompare peptide_id_compare_result_file!!!\n" +

                               "All the ms1 and ms2 files should be linked to the current working directory");
            System.exit(0);
        }
        String idcompareresultfile = args[0];
	reader.precursorIntensity(idcompareresultfile);
    }
}

