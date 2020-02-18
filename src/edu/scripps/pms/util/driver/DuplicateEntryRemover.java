/**
 * @file DuplicateEntryRemover.java
 * This is the source file for edu.scripps.pms.util.spectrum.DuplicateEntryRemover
 * @author Tao Xu
 * @date $Date: 2007/01/08 23:18:09 $
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
public class DuplicateEntryRemover {
    public static final String USAGE = "DuplicateEntryRemover inputFile outputFile";
    public static void main(String args[]) throws Exception {
        try {
            String inputFile = args[0];
            BufferedReader br = new BufferedReader(new FileReader(inputFile));
            
            PrintStream ps = System.out; 
            if(args.length > 1) {
                String outputFile = args[1];
                ps = new PrintStream(outputFile);
            }
            HashSet<String> dirs = new HashSet<String>(5000);
            String line = null;
            int numLines = 0;
            while((line = br.readLine()) != null) {
                if(!line.equals("")) {
                    numLines++;
                    if(dirs.add(line)) {
                        ps.println(line);
                    } else {
                        System.out.println("Duplicate entry: " + line);
                    }
                }
            }
            int numDirs = dirs.size();
            System.out.println("Number of duplicate entries: " + (numLines - numDirs));
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
