import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import java.io.*;

import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.mspid.MassCalculator;
import edu.scripps.pms.mspid.MassSpecConstants;
import edu.scripps.pms.util.enzyme.Protease;

/**
 * <p>Program to six-frame translate a nucleotide sequence</p>
 */

// It now ignores frames with more than 2 stop condons
public class AddProteinLengthForAnita {



    /**
     * Call this to get usage info, program terminates after call.
     */
    public static void help() {
	System.out.println(
		"usage: java AddProteinLengthForAnita <DTArray> <fastafile>");
	System.exit( -1);
    }

    private static HashMap<String, Integer> mapGene2Length(String dtarrayfile, String fastaFileName) throws IOException {
        
	    BufferedReader br = new BufferedReader(new FileReader(dtarrayfile));
            String line = br.readLine();
            line = br.readLine();
            HashSet<String> genenames = new HashSet();
            HashMap<String, Integer> gene2length = new HashMap(); 
 
            while(line != null) {
                line = line.trim();
                if(!"".equals(line)) {
                    String genename = line.split("\t")[1];
System.out.println("Genename added: " + genename);
                    genenames.add(genename);
                } 
                line = br.readLine();
            } 
             
System.out.println("Number of genenames added: " + genenames.size());
            for (Iterator<Fasta> it = FastaReader.getFastas(fastaFileName); it.hasNext();) {
                Fasta f = it.next();
                String gene = f.getGeneName();
                if(genenames.contains(gene)) {
                    //System.out.println(acc + " is in list");
                    gene2length.put(gene, f.getLength());
                } else {

                    //System.out.println(acc + " not is in list");
                }
                
            }
            return gene2length;
    }

    public static void main(String[] args) throws Exception{

        
	BufferedReader br = null;
         

//	try {

            String dtarrayfile = args[0];
            String fastaFileName = args[1];

            HashMap<String, Integer> genename2length = mapGene2Length(dtarrayfile, fastaFileName);
//System.out.println("genename2length: " + genename2length);
System.out.println("Number of elements:  " + genename2length.size());

	    br = new BufferedReader(new FileReader(dtarrayfile));
            System.out.println(br.readLine());

            String line = br.readLine();
            System.out.println(line);
             
            line = br.readLine();
            while(line != null) {
//System.out.println("now processing before trim: " + line);
                line = line.trim();
System.out.println("now processing: " + line);
                if(!"".equals(line)) {
                    String [] arr = line.split("\t");
//System.out.println("number of element: " + arr.length);
                    String genename = arr[1];
System.out.println("genename:  "+ genename);
                    Integer len = genename2length.get(genename);
                    int length = 0;
                    if(len != null) {
                        length = len;
                    } else {
System.out.println("Length is not found for " + genename);
                    }
System.out.println("genename:  "+ genename + "\tlength" + length);
                    StringBuffer sb = new StringBuffer();
                    sb.append(arr[0]);
                    sb.append("\t");
                    sb.append(arr[1]);
                    sb.append("\t");
                    sb.append("" + length);
                    for(int i = 3; i < arr.length; i++) {
                        sb.append("\t");
                        sb.append(arr[i]);
                       
                    } 
                    System.out.println(sb); 
                } 
                line = br.readLine();
            } 
           
             
            

//	} catch(Exception e){
//            help();
            
//            e.printStackTrace(); 
//	} finally {
	    //tidy up
//	    if(br != null){
//		br.close();
//	    }
//	}
    }
}

