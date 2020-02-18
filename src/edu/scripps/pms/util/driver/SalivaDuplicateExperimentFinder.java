/**
 * @file DuplicateEntryRemover.java
 * This is the source file for edu.scripps.pms.util.spectrum.DuplicateEntryRemover
 * @author Tao Xu
 * @date $Date: 2007/03/08 21:46:01 $
 */



import edu.scripps.pms.util.io.*;

//import edu.scripps.pms.util.io.DTASelectFilterReader; 
import edu.scripps.pms.util.dtaselect.Peptide; 
import edu.scripps.pms.util.dtaselect.Protein; 
import java.io.*;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import java.util.HashSet;
import edu.scripps.pms.util.spectrum.*;

// this is the class for Akira's  
public class SalivaDuplicateExperimentFinder {
    public static final String USAGE = "DuplicateEntryRemover inputFile outputFile";
    public static List<String> getFiles(String dir) {
        
         ArrayList<String> files = new ArrayList<String>();
         File currentDir = new File(dir);
         for(String s : currentDir.list()) {
             if(s.endsWith(".ms2")) {
                 files.add(s);
                 //System.out.println(s);
             }
         }
         return files;
    }
    public static void main(String args[]) throws Exception {
        try {
            String inputFile = args[0];
            BufferedReader br = new BufferedReader(new FileReader(inputFile));
            
            PrintStream ps = System.out; 
            HashSet<String> dirs = new HashSet<String>(5000);
            String line = null;
            int numLines = 0;
            while((line = br.readLine()) != null) {
                if(!line.equals("")) {
                    numLines++;
                    if(dirs.add(line)) {
         //               ps.println(line);
                    } else {
                        System.out.println("Duplicate Director: " + line);
                    }
                }
            }
            int numDirs = dirs.size();
            System.out.println("Number of duplicate entries: " + (numLines - numDirs));

            HashMap <String, String> ms2File2Dir = new HashMap<String, String>(10000);
            HashSet<String> dupLines = new HashSet<String>();
            for(String dir : dirs) {
                for(String ms2File : getFiles(dir)) {
                    String dirExist = ms2File2Dir.get(ms2File);
                    if(dirExist == null) {
                        ms2File2Dir.put(ms2File, dir);
                    } else {
                        dupLines.add("Duplicate ms2 file found in " + dir + "\t" + dirExist);
                    }
                }
            }
            for(String dup : dupLines) {
                System.out.println(dup); 
            } 

            br.close();
            ps.close();
        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
            System.out.println(USAGE);
        }
    }
}
