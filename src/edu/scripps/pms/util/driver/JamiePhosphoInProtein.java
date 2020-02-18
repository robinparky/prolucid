import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import edu.scripps.pms.util.seq.Fasta;


/**
 * <p>Program to six-frame translate a nucleotide sequence</p>
 */

// It now ignores frames with more than 2 stop condons
public class JamiePhosphoInProtein {


    //public static final String phosphoresultfile = "/data/2/rpark/ip2_data/catclw/Jamie/2008_07_asyn_bluedot_c_2011_01_07_20_2904/search/upload2011_01_10_13_8781/phospho/phosphoresult.txt";
    //public static final String phosphoresultfile = "/data/3/taoxu/projects/catherine/jamie/phosphoresults.txt";
    //public static final String phosphoresultfile = "/data/3/taoxu/projects/catherine/jamie/all_phosphoresults_20110118.txt";
    public static final String phosphoresultfile = "/data/3/taoxu/projects/catherine/jamie/all_phosphoresults_20110119.txt";
    public static final String database = "/garibaldi/people-b/applications/yates/dbase/SGD_S-cerevisiae_na_12-16-2005_con_reversed.fasta";
    public static final String interestedacsfile = "/data/3/taoxu/projects/catherine/jamie/spblist.txt";

    /**
     * Call this to get usage info, program terminates after call.
     */
    public static void help() {
	System.out.println(
		"usage: java JamiePhosphoInProtein ");
	System.exit( -1);
    }

    public static int getPeptideLocation(String peptide, Fasta protein) {
        return protein.getSequence().indexOf(peptide);
    }

    public static void main(String[] args) throws Exception{

	BufferedReader br = null;
        HashSet<String> interestedacs = new HashSet();
        HashMap<String, ArrayList<PhosphoResult>> ac2results = new HashMap();

	try {
	    br = new BufferedReader(new FileReader(interestedacsfile));
            String line = br.readLine();
             
            while(line != null) {
                line = line.trim();
                if(!"".equals(line)) {
                    interestedacs.add(line);
                    ac2results.put(line, new ArrayList<PhosphoResult>());
                } 
                line = br.readLine();
            } 
            br.close();
            
	    br = new BufferedReader(new FileReader(phosphoresultfile));
            line = br.readLine();
            while(line != null) {
                line = line.trim();
                if(!"".equals(line)) {
                    for(Iterator<String> it = interestedacs.iterator(); it.hasNext();) { 
                        String ac = it.next();
                        if(line.indexOf(ac) > -1 ) { 
                            PhosphoResult pr = new PhosphoResult(line);
                            if(!pr.getProteinAcc().startsWith("Reverse")) {
                                ac2results.get(ac).add(new PhosphoResult(line));
                            }
                        }
                    }
                } 
                line = br.readLine();
            } 
            br.close();

            for(Iterator<String> it = interestedacs.iterator(); it.hasNext();) { 
                String ac = it.next();
                
                ArrayList<PhosphoResult> result = ac2results.get(ac);
                System.out.println("\n\n\nNumber of phospho ids for " + ac + " is " + result.size());
                //System.out.println("Sequence\tAscore\tdebunker score\tDTASelect confidence score\tXcorr\tDeltaCN\tScan\tFragmentation Method\tcharge\tSpectrum file\tProtein accession number\tDescription");
                System.out.println("Ascore Sequence\tOriginal Sequence\tAscore\tdebunker score\tDTASelect confidence score\tXcorr\tDeltaCN\tScan\tFragmentation Method\tcharge\tSpectrum file\tProtein accession number\tDescription");
                for(Iterator<PhosphoResult> it1 = result.iterator(); it1.hasNext();) {
                    System.out.println(it1.next().getLine());
                }
            }
            

	} finally {
	    //tidy up
	    if(br != null){
		br.close();
	    }
	}
    }
}

