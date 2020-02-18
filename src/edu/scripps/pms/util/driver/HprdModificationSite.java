
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

import java.util.*;
import java.io.*;
import edu.scripps.pms.util.seq.Fasta;

// this program is used to output HPRD PTM file.  It takes a DTASelect-mod.txt file as input 

public class HprdModificationSite {
    private static final String USAGE = "\n\nUsage: java HprdModificationSite DTASelect-mod.txtFile"; 
    
    // to store the first line of the input file

//    public static HashMap<String, String> acc2Protein = new HashMap<String, String>(2000); 
//    public static ArrayList<HprdModificationSites> file1Entries = new ArrayList<HprdModificationSites>();
 //   public static ArrayList<HprdModificationSites> file2Entries = new ArrayList<HprdModificationSites>();
    private String accession;
    private String description;
    
    private int site;
    private char residue;
    private double score;

    private String spectrumNumber;
    private String experimentType = "Mass Spectrometry";
    private String upstreamEnzyme = "Unknown";
    private String ptmType = "Phosphorylation";
    private String program = "SEQUEST";
    public String getHprdModification() {
        return accession + "\t" + ptmType + "\t" + site + "\t" + residue + "\t"
                + experimentType + "\t" + upstreamEnzyme + "\t" + score + "\t" + program 
                + "\t" + spectrumNumber;
    }
    public HprdModificationSite(String lLine, String mLine, String nLine) {
        
        String [] elements = lLine.split("\t");
        accession = Fasta.getAccession(elements[1]);
        description = elements[1];
        
        elements = mLine.split("\t");
        score = Double.parseDouble(elements[4]); 
        spectrumNumber = elements[5];
        
        elements = nLine.split("\t");
        site = Integer.parseInt(elements[1]);
        residue = elements[4].charAt(0);
        //System.out.println(line);
        //System.out.println("numElement: " + elements.length);
    
    }  
    public String getAccession() {
        return accession;
    }
    public boolean equals(HprdModificationSite h) {
        if(h == null) {
            return false;
        }
        
        return  (residue == h.residue) && site == h.site && accession.equals(h.accession);
    }

    
    private static void readModifications(String inputFile) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(new File(inputFile)));
        br.readLine(); // remove the header line
        br.readLine();
        br.readLine(); // remove the header line
        br.readLine();
 
        String line = null; 
        String lLine = null;
        String mLine = null; 
        String nLine = null; 
        int numPeptides = 0;
        HashSet<Integer> scans = new HashSet<Integer>();
        System.out.println("Sequence identifier\tPTM type\tSite\tResidue\tExperiment type\tUpstream enzyme\tPeptide score\tAlgorithms used\tSpectrum");
        while((line = br.readLine()) != null) {
            boolean isFirstSite = false;
            if(line.startsWith("L")) {
                lLine = line;
                 isFirstSite = true;
                line = br.readLine();
            } 

            while (!line.startsWith("M")) {
                line = br.readLine();
            }

            mLine = line;
            while(!line.startsWith("N")) {
                line = br.readLine();
            }

            nLine = line; 
            HprdModificationSite h = new HprdModificationSite(lLine, mLine, nLine);
            System.out.println(h.getHprdModification()); 
            if(isFirstSite) {
                numPeptides++;
                String [] elements = h.spectrumNumber.split("\\.");
//System.out.println("spectrumNumber: " + h.spectrumNumber + "\t" + elements);
                scans.add(new Integer(elements[1]));         
            }
        }
System.out.println("numPeptides: " + numPeptides + "\t" + "numOfScanNumbers: " + scans.size());
        br.close();
    }
    public static void main(String args[]) {
        try {
            readModifications(args[0]);
            
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
        }
     
    }

}
