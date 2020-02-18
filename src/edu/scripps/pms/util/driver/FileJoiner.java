/**
 * @file FileJoiner.java
 * This is the source file for edu.scripps.pms.util.spectrum.FileJoiner
 * @author Tao Xu
 * @date $Date
 */



import edu.scripps.pms.util.io.*;

//import edu.scripps.pms.util.io.DTASelectFilterReader; 
import edu.scripps.pms.util.dtaselect.Peptide; 
import edu.scripps.pms.util.dtaselect.Protein; 
import edu.scripps.pms.util.io.SpectrumReader; 
import java.io.*;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

// comparing retention times of peptides identified in the sampe salt step in different MudPIT runs 
// to estimate the reproducibility of each salt step of LC/LC 
// this program expect a text file that contains all the directory names for comparison
// it also expect the number of salt steps the use wants to compare
public class FileJoiner {
    public static final String USAGE = "java FileJoiner fileName1 fileName2 [join_key_index]";
    private static String header = null;
    public static ArrayList<String> readFile(String file) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(file), 4096);
        
//System.out.println("Processing file: " + file);
        ArrayList<String> map = new ArrayList<String>(100000); 
        String line = null;
        while((line = br.readLine()) != null) {
           // System.out.println("key: " + key + "\tline: " + line);
            map.add(line);    
        }
//System.out.println("NumLines Added: " + numLines + "\tNumKeys: " + map.size());
        br.close();
        return map; 
    }

    public static int processFile(String file, int keyIndex, ArrayList<String> keys, HashMap<String, String> map) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(file), 4096);
        boolean isHeaderLine = true; 
//System.out.println("Processing file: " + file);
        String line = null;
        int numColumns = 0;
        while((line = br.readLine()) != null) {
//System.out.println("line: " + line);
            if(isHeaderLine) {
                numColumns = line.split("\t").length;
                if(header == null) {
                    header = line + "\t";
                } else {
                    header += line;
                }          
                isHeaderLine = false;
            } else {
                String [] arr = line.split("\t");
                String key = arr[keyIndex]; 
                //System.out.println("key: " + key + "\tline: " + line);
                map.put(key, line);  
                keys.add(key);
            }  
        }
//System.out.println("NumLines Added: " + numLines + "\tNumKeys: " + map.size());
        br.close();
        return numColumns;
    }
    public static void main(String args[]) throws Exception {
         
        try {
            HashMap<String, String> firstFileMap = new HashMap<String, String>(100000); 
            HashMap<String, String> secondFileMap = new HashMap<String, String>(100000); 
            ArrayList<String> firstFileKeys = new ArrayList<String>(100000);
            ArrayList<String> secondFileKeys = new ArrayList<String>(100000);

            int keyindex = args.length > 2? Integer.parseInt(args[2]) - 1 : 0;

            int firstFileNumColumns = processFile(args[0], keyindex, firstFileKeys, firstFileMap); 
            int secondFileNumColumns = processFile(args[1], keyindex, secondFileKeys, secondFileMap); 
            String firstFileTabs = "";
            for(int i = 0; i < firstFileNumColumns; i++) {
                firstFileTabs += "\t";
            } 
            String secondFileTabs = "";
            for(int i = 0; i < secondFileNumColumns; i++) {
                secondFileTabs += "\t";
            }
           
            System.out.println(header);
            // output all rows in file1 and rows in file2 that the keys were found in file1
            for(String key : firstFileKeys) {
                //System.out.print("key: " + key + "\t" + firstFileMap.get(key) + "\t");
                System.out.print(firstFileMap.get(key) + "\t");
                String linetoadd = secondFileMap.get(key);
                if(linetoadd != null) {
                    System.out.println(secondFileMap.get(key) + "\t0");
                } else {
                    System.out.println(secondFileTabs + "1");
                }
            }
            // output rows in file2 that the keys were exlusively found in file2
            for(String key : secondFileKeys) {
                String linekeep = firstFileMap.get(key);
                String lineadd = secondFileMap.get(key); 
                if(linekeep == null) {
                    System.out.println(firstFileTabs + lineadd + "\t2");
                }
            }          
 
          /* 
            for(String line : toKeep) {
                String [] arr = line.split("\t");
                String key = arr[3];
                //System.out.print("key: " + key + "\t" + firstFileMap.get(key) + "\t");
                System.out.print(firstFileMap.get(key) + "\t");
                String linetoadd = secondFileMap.get(key);
                if(linetoadd != null) {
                    System.out.println(secondFileMap.get(key));
                } else {
                    System.out.println();
                }
            }
            */
        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
            System.out.println(USAGE);
        }
    }
}
