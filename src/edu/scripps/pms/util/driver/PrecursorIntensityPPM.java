import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Collections;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.io.*;


// add precursor intensity to peptide id compare result
// all the ms2 files should be in the current working dir
// 
public class PrecursorIntensityPPM
{
    private static HashMap<String, double[]> file2intensities = new HashMap();
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
	        for(int i = 1; i < temp.length; i+=13)
	        {
	    	    returnStr = getPrecursorIntensityPPM(temp[i]);
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

    public String getPrecursorIntensityPPM(String str)
    {
        // 110525_HEK_DigDeAPr_75Percent_125ug_100ugmL_Step08.4279.4279.2
        // if str equal "X" or "" returh "X", otherwise, 
//System.out.println("filescanline: " + str);
        
        if(str == null || "X".equals(str.trim()) || "".equals(str.trim())) {
            return "X";
        }

        String [] arr = str.split("\\.");
        if(arr.length < 2) return "";
        String filename = arr[0];
        int scannum = Integer.parseInt(arr[1]);

//System.out.println("Filename: " + filename + "\tscan number: " + scannum);
        double [] intens = file2intensities.get(filename);
        if(intens == null) {
            intens = new double[100000];
            file2intensities.put(filename, intens);
            try{
                SpectrumReader sr = new SpectrumReader(filename+".ms2", ".ms2");
                Iterator<PeakList> it = sr.getSpectra();
                while(it.hasNext()) {
                    PeakList pl = it.next();
                    int scan = pl.getLoscan();
                    double prcint = pl.getPrecursorInt();
                    intens[scan] = prcint;

                }
                sr.closeDataFile();
            } catch(IOException e) {
                e.printStackTrace();
            }
        }
        
//System.out.println("return: " + intens[scannum]);
        return ""+intens[scannum];
        
    }

    public static void main(String[] args)
    {
	PrecursorIntensityPPM reader = new PrecursorIntensityPPM();
        //String idcompareresultfile = "/home/mtlin/mtlin_on_data/precursorInstensity/ident_compare_PEPTIDE1703.txt";
        if(args.length < 1) {
            System.out.println("!!!Usage: addPrecursorIntensityToPeptideIdentCompare peptide_id_compare_result_file!!!\n" +

                               "All the ms2 files should be linked to the current working directory");
            System.exit(0);
        }
        String idcompareresultfile = args[0];
	reader.precursorIntensity(idcompareresultfile);
    }
}

