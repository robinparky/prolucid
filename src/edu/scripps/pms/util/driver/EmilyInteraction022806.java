
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

public class EmilyInteraction022806 {
    private static final String USAGE = "\n\nUsage: java EmilyInteraction022806 fileName"; 
    private String entryLine;
    
    // to store the first line of the input file
    private String firstLine;
    private String [] elements; // first line elements
    private ArrayList<HashSet<String>> fractionList;

    public static final int ACCPOSITION = 0;
    public static final int NUMOTHERELEMENTS = 5; // # of leading elements in each line
    public static final String DELIMITER = "\t"; // # of leading elements in each line
    
    
    public EmilyInteraction022806(String fileName) throws IOException {
         
        BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
        firstLine = br.readLine();
        elements = firstLine.split(DELIMITER);
        int numElements = elements.length;

        fractionList = new ArrayList<HashSet<String>>();

        for(int i = 0; i < numElements; i++) {
            fractionList.add(new HashSet<String>());
        }
        String line = br.readLine(); 
        while(line != null) {
            String [] e = line.split(DELIMITER);
            for(int i = NUMOTHERELEMENTS; i < numElements; i++) {
                String count = e[i].trim();
                if(!(count.equals("x") || count.equals("") || count.equals(" "))) {
                    fractionList.get(i).add(e[0]);
                } else {
//                    System.out.println(count); 
                }
            }     
            line = br.readLine();
        } 
        br.close();
    
    }  
    public int getNumElements() {
        return elements.length;
    } 
    public String index2FractionName(int index) {
        return elements[index];
    }
    public HashSet<String> index2FractionSet(int index) {
        HashSet<String> result = new HashSet<String>();
        result.addAll(fractionList.get(index)); 
        result.addAll(fractionList.get(index+1));
        result.addAll(fractionList.get(index+2));

        return result; 
    }
    public static void main(String args[]) {
        try {
            String inputFile1 = args[0];
            String inputFile2 = args[1];
            EmilyInteraction022806 method1 = new EmilyInteraction022806(inputFile1);
            EmilyInteraction022806 method2 = new EmilyInteraction022806(inputFile2);
            
            String input1 = inputFile1.substring(0, inputFile1.indexOf("."));
            String input2 = inputFile2.substring(0, inputFile2.indexOf("."));
            String outputFile = "interactionsIn" + input1 + "And" + input2 + ".txt";
            PrintWriter outFile = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
            int startIndex = NUMOTHERELEMENTS;
            int numElements1 = method1.getNumElements();
            int numElements2 = method2.getNumElements();

            for(int i = startIndex; i < numElements1; i += 3) {
                String fractionName1 = method1.index2FractionName(i);
                HashSet<String> set1 = method1.index2FractionSet(i);
                //System.out.println("i: " + i + "\t numOfElement in set1: " + set1.size());
                for(int j = startIndex; j < numElements2; j += 3) {
                    String fractionName2 = method2.index2FractionName(j);
                    HashSet<String> set2 = method2.index2FractionSet(j);
                    //System.out.println("\tj: " + j + "\t numOfElement in set2: " + set2.size());
                    Set<String> intersection = SetOperator.getIntersection(set1, set2);
                    if(intersection.size() > 1) {
                        StringBuffer sb = new StringBuffer();
                        String fractions = fractionName1 + "_" + fractionName2 + "\t" + intersection.size() + "\t";
                        sb.append(fractions);
                        Iterator it = intersection.iterator(); 
                        while(it.hasNext()) {
                            sb.append(it.next());
                            sb.append(" ");
                        }
                        outFile.println(sb);
                    }
                }
            }
            outFile.close();
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
        }
     
    }

}
