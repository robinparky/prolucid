/**
 * @file DtaselectTxtFilter.java
 * This is the source file for edu.scripps.pms.util.spectrum.DtaselectTxtFilter
 * @author Tao Xu
 * @date $Date: 2011/01/29 00:13:58 $
 */



import edu.scripps.pms.util.io.*;

//import edu.scripps.pms.util.io.DTASelectFilterReader; 
import edu.scripps.pms.util.dtaselect.Peptide; 
import edu.scripps.pms.util.dtaselect.Protein; 
import edu.scripps.pms.util.PmsUtil; 
import edu.scripps.pms.util.stat.StatCalc; 
import java.io.*;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import java.util.HashSet;
import edu.scripps.pms.util.spectrum.*;

// command line processor
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

// this is the class for calculating the z-score for delta mass 
public class DtaselectTxtFilter {

  
    public static void main(String [] args) throws Exception {
        String dtatxtfile = "/data/3/taoxu/projects/bingwen/salliva/smslmale/DTASelect.txt";
        String filestoretain = "/data/3/taoxu/projects/bingwen/salliva/smslmale/smslsqtfiles-male-only.txt";
        HashSet<String> fileswanted = new HashSet<String>(1000000);

        BufferedReader br = new BufferedReader(new FileReader(filestoretain));

        String line = null;
        while ((line = br.readLine()) != null) {
             fileswanted.add(line.trim());
        }
        //System.out.println();
        br.close();

        br = new BufferedReader(new FileReader(dtatxtfile));

        while ((line = br.readLine()) != null) {
             if(line.startsWith("D\t")) {

                 //System.out.println("line " + line);
                 String file = line.split("\t")[1].split("\\.")[0];
                 //System.out.println(file);
                 if(fileswanted.contains(file)) {
                     //System.out.println("we want to keep " + file);
                     System.out.println(line);
                 } else {

                     //System.out.println("we want to remove " + file);
                 }
                 
             } else {
                 System.out.println(line);
             }
        }
        br.close();

    }


}
