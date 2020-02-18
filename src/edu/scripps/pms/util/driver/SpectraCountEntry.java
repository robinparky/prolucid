
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

// this program is used to sort Lujian's spectra count file

public class SpectraCountEntry implements Comparable<SpectraCountEntry> {
    private static final String USAGE = "\n\nUsage: java SpectraCountEntry fileName"; 
    private String entryLine;
    private double [] countValues = new double[15];
    private double avgRatio;
    private int numNone0 = 0;
    
    // to store the first line of the input file
    private static String firstLine;
    private static String inputFile;
    private static PrintWriter outFile;
    private static ArrayList<SpectraCountEntry> allThree = new ArrayList<SpectraCountEntry>(1000);
    private static ArrayList<SpectraCountEntry> oneAnd2 = new ArrayList<SpectraCountEntry>(1000);
    private static ArrayList<SpectraCountEntry> oneAnd3 = new ArrayList<SpectraCountEntry>(1000);
    private static ArrayList<SpectraCountEntry> twoAnd3 = new ArrayList<SpectraCountEntry>(1000);
    private static ArrayList<SpectraCountEntry> one = new ArrayList<SpectraCountEntry>(1000);
    private static ArrayList<SpectraCountEntry> two = new ArrayList<SpectraCountEntry>(1000);
    private static ArrayList<SpectraCountEntry> three = new ArrayList<SpectraCountEntry>(1000);
    private static ArrayList<SpectraCountEntry> allZeros = new ArrayList<SpectraCountEntry>(1000);

    public int compareTo(SpectraCountEntry e1) {
        double diff = avgRatio - e1.avgRatio;
        if(diff > 0) {
            return -1;
        } else if (diff < 0) {
            return 1;
        } else {
            return 0;
        }
    }
    public SpectraCountEntry(String line) {
        entryLine = line;
        
        String [] elements = entryLine.split("\t");
        //System.out.println(line);
        //System.out.println("numElement: " + elements.length);
    
        for(int i = 3; i < elements.length; i++) {
            countValues[i] = Double.parseDouble(elements[i].trim());
        }   
        classifyEntry();
        entryLine += avgRatio;
    }  
    private void addCount(double ratio) {
        if(ratio != 0) {
            numNone0++;
            avgRatio += ratio;
        }
    }
    // return avg ratio
    private void classifyEntry() {

        addCount(countValues[6]);
        addCount(countValues[10]); 
        addCount(countValues[14]);

        if(countValues[6] != 0) {
            if(countValues[10] != 0) {
                if(countValues[14] != 0) {
                    allThree.add(this); 
                } else {
                    oneAnd2.add(this);
                }
            } else { // two is 0
                if(countValues[14] != 0) {
                    oneAnd3.add(this); 
                } else { // three is 0
                    one.add(this); // one only
                }
            }
        } else { // one is 0
            if(countValues[10] != 0) { // two is not 0, one is 0
                if(countValues[14] != 0) {
                    twoAnd3.add(this); 
                } else {
                    two.add(this);
                }
            } else { // two is 0
                
                if(countValues[14] != 0) {
                    three.add(this); 
                } else { // three is 0
                    allZeros.add(this); // all 0 
                }
            }
        }
        if(numNone0 > 0) {
            avgRatio = avgRatio/numNone0;
        }
    } 
    public static void outputResult(String outFileName) throws IOException {
        Collections.sort(one);
        Collections.sort(two);
        Collections.sort(three);
        Collections.sort(allThree);
        Collections.sort(oneAnd2);
        Collections.sort(oneAnd3);
        Collections.sort(twoAnd3);
        
        outFile.println(firstLine);  
        outputList(allThree); 
        outputList(oneAnd2); 
        outputList(oneAnd3); 
        outputList(twoAnd3); 
        outputList(one); 
        outputList(two); 
        outputList(three); 
        outputList(allZeros); 
       
        outFile.close();
    }
    private static void outputList(List<SpectraCountEntry> list) {
        Iterator<SpectraCountEntry> it = list.iterator();
        while(it.hasNext()) {
            outFile.println(it.next().entryLine);
        }    
    }
    public static void main(String args[]) {
        try {
            String inputFile = args[0];
            String outputFile = "sorted" + inputFile;
            BufferedReader br = new BufferedReader(new FileReader(new File(inputFile)));
            firstLine = br.readLine();
            firstLine += "avgRatio"; 
            String line = br.readLine();
            while(line != null) {
                new SpectraCountEntry(line);  
                line = br.readLine();
            } 
            br.close();
            outFile = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
            outputResult(outputFile);
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
        }
     
    }

}
