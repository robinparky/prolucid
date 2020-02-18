import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import java.io.*;

import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;

/**
 * <p>Program to six-frame translate a nucleotide sequence</p>
 */

// It now ignores frames with more than 2 stop condons
public class UniqueTmtOutput {



    /**
     * Call this to get usage info, program terminates after call.
     */
    public static void help() {
	System.out.println(
		"usage: java UniqueTmtOutput input_file");
	System.exit( -1);
    }

    private static double calcTotalIntensity(String out) {
        double totalinten = 0;
        if(out == null) return totalinten;
        String arr[] = out.split("\t");
        for(int i = 4; i < 10; i++) {
            if(arr[i] != null && !"".equals(arr[i])) {
                totalinten += Double.parseDouble(arr[i]);
            }
        }
 
        return totalinten;
    }


    public static void main(String[] args) throws Exception{

        
	BufferedReader br = null;
        HashMap<String, String> ac2output = new HashMap();

	try {

            String tmtfile = args[0];

	    br = new BufferedReader(new FileReader(tmtfile));
            String line = br.readLine();
             
            while(line != null) {
                line = line.trim();
                if(!"".equals(line)) {
                    String acc = line.split("\t")[1];
                   
                    String temp = ac2output.get(acc);
                    if(temp == null) {
                        ac2output.put(acc, line);
                    } else {
                        double intnew = calcTotalIntensity(line);
                        double intold = calcTotalIntensity(temp);
                        if(intnew > intold) {
                             ac2output.put(acc, line);
                        }
                    }
                } 
                line = br.readLine();
            } 

            Iterator<String> accset = ac2output.keySet().iterator();
            while(accset.hasNext()) {
                System.out.println(ac2output.get(accset.next()));
            }
	} catch(Exception e){
            help();
            e.printStackTrace(); 
	} finally {
	    //tidy up
	    if(br != null){
		br.close();
	    }
	}
    }
}

