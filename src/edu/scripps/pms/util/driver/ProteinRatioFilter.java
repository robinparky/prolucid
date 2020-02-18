
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

// this program is used to sort Lujian's spectra count file

public class ProteinRatioFilter {
    private static final String USAGE = "\n\nUsage: java ProteinRatioFilter dummyFileName file2 upRatioCutOff minPeptideCount minSpectrumCount"; 
    private static int minSpectrumCount = 0;
    private static int minPeptideCount = 0;
    private static double upRatioCutOff;
    private static double downRatioCutOff;
    
    // to store the first line of the input file
    private static String firstLine;
    private static String inputFile;

    private static PrintWriter outFile;
//    public static HashMap<String, String> acc2Protein = new HashMap<String, String>(2000); 
//    public static ArrayList<ProteinRatioFilter> file1Entries = new ArrayList<ProteinRatioFilter>();
 //   public static ArrayList<ProteinRatioFilter> file2Entries = new ArrayList<ProteinRatioFilter>();
    private String entryLine;
    private String [] elements;
    private String accession;
    
    public ProteinRatioFilter(String line) {
    
        entryLine = line;
        
        elements = entryLine.split("\t");
        String acc = Fasta.getAccession(elements[0]);
        int index = acc.indexOf("."); 
        if(index != -1) {
            acc = acc.substring(0, index);
        }
        accession = acc; 
        //System.out.println(line);
        //System.out.println("numElement: " + elements.length);
    
    }  
    public String getAccession() {
        return accession;
    }
    public int getNumElements() {
        return elements.length;
    }
    public String getEntryLine() {
        return entryLine;
    }
    public String getElement(int index) {
        return elements[index];
    }
    private static List<ProteinRatioFilter> readFile1(String inputFile) throws IOException {
        ArrayList<ProteinRatioFilter> file1Entries = new ArrayList<ProteinRatioFilter>();
        BufferedReader br = new BufferedReader(new FileReader(new File(inputFile)));
        String line = br.readLine(); // remove the header line
        while((line = br.readLine()) != null) {
            ProteinRatioFilter prf = new ProteinRatioFilter(line); 
            if(prf.getNumElements() < 17) {
                continue;
            }
            double ratio1 = Double.parseDouble(prf.getElement(6)); 
            double ratio2 = Double.parseDouble(prf.getElement(11)); 
            double ratio3 = Double.parseDouble(prf.getElement(16)); 
            int numRatioGood = 0;
            
            // need to consider ratio <= downRatioCutOff ? 
            numRatioGood += ratio1 > upRatioCutOff? 1 : 0; 
            numRatioGood += ratio2 > upRatioCutOff? 1 : 0; 
            numRatioGood += ratio3 > upRatioCutOff? 1 : 0; 
        
            if(numRatioGood >= 1) {
                file1Entries.add(prf);
            }
        }
        br.close();
        return file1Entries;
    }
    private static List<ProteinRatioFilter> readFile2(String inputFile) throws IOException {
        
        PrintWriter allAccs = new PrintWriter(new BufferedWriter(new FileWriter("allAccs.txt")));
        PrintWriter changedAccs = new PrintWriter(new BufferedWriter(new FileWriter("changedAccs.txt")));

        ArrayList<ProteinRatioFilter> file2Entries = new ArrayList<ProteinRatioFilter>();
        BufferedReader br = new BufferedReader(new FileReader(new File(inputFile)));
        String line = br.readLine(); // remove the header line
        while((line = br.readLine()) != null) {
            ProteinRatioFilter prf = new ProteinRatioFilter(line); 
            if(prf.getNumElements() < 9) {
               continue;
            }
            double ratio1 = Double.parseDouble(prf.getElement(10)); 
//System.out.println("ratio: " + ratio1);
            allAccs.println(prf.getAccession());
            // need to consider ratio1 <= downRatioCutOff? 
            boolean isGood = false;
            if(ratio1 >= upRatioCutOff) {
                if(Integer.parseInt(prf.getElement(2)) == 0) {
                    if(Integer.parseInt(prf.getElement(4)) > minPeptideCount && 
                             Integer.parseInt(prf.getElement(5)) > minSpectrumCount) {
                        isGood = true;
                    } 
                } else {
                    isGood = true;
//System.out.println("ratio: " + ratio1);
                }
            }
            if(ratio1 <= downRatioCutOff) {
                if(Integer.parseInt(prf.getElement(4)) == 0) {
                    if(Integer.parseInt(prf.getElement(2)) > minPeptideCount && 
                             Integer.parseInt(prf.getElement(3)) > minSpectrumCount) {
                        isGood = true;
                    } 
                } else {
//System.out.println("ratio: " + ratio1);
                    isGood = true;
                    
                }
            }
            if(isGood) {
                file2Entries.add(prf);
                changedAccs.print(prf.getAccession());
                if(ratio1 >= upRatioCutOff) {
                    changedAccs.println("\t" + "1");
                } else {
                    changedAccs.println("\t" + "-1");

                }
            }
        }
        br.close();
        changedAccs.close();
        allAccs.close();
        return file2Entries;
    }
    public static void outputResults(List<ProteinRatioFilter> file1Entries, List<ProteinRatioFilter> file2Entries) {
        HashSet<String> accsInFile1 = new HashSet<String>(1000);
        for(ProteinRatioFilter prf : file1Entries) {
            if(accsInFile1.contains(prf.getAccession())) {
                //System.out.println("Protein already included: " + prf.getAccession());
            }
            accsInFile1.add(prf.getAccession());
            
        } 
        System.out.println("Number of accessions in file1: " + accsInFile1.size());  
     
        for(ProteinRatioFilter prf : file2Entries) {
            String acc = prf.getAccession(); 
            if(accsInFile1.contains(acc)) {
                //System.out.println(acc + "\t" + prf.getEntryLine());
                System.out.println(prf.getEntryLine());
            }
        } 
    }
    public static void main(String args[]) {
        try {
            String file1 = args[0];
            String file2 = args[1];
            upRatioCutOff = Double.parseDouble(args[2]);
            downRatioCutOff = 1/upRatioCutOff; 
            minPeptideCount = Integer.parseInt(args[3]) - 1;
            minSpectrumCount = Integer.parseInt(args[4]) - 1;
          
     
            
//            List<ProteinRatioFilter> file1Entries = readFile1(file1);
 //           System.out.println("NumProteins in " + file1 + ": " + file1Entries.size());
            List<ProteinRatioFilter> file2Entries = readFile2(file2);
            System.out.println("NumProteins in " + file2 + ": " + file2Entries.size());
            
            outputResults(file2Entries, file2Entries);
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
        }
     
    }

}
