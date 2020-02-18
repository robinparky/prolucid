/**
 * @file Accession2Fastas.java
 * This is the source file for Accession2Fastas class
 * @author Tao Xu
 * @date $Date
 */


import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.io.FastaReader;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.HashMap;
import java.util.Set;
import java.io.FileInputStream;
import java.io.File;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.IOException;
import java.io.*;
import java.util.*;


public class FastaFilter {

    public static final String USAGE = "!!! Usage: fastafiler databaseFile to_include1 to_include2 ... !!!";
    // databaseName - the path and name of the protein database

    // assuming the mass tolerance is smaller than the mass of any AA residue
    public static void main(String [] args) throws Exception {
        try {
            String databaseName = args[0];
            ArrayList<String> toInclude = new ArrayList<String>();
            for(int i = 1; i < args.length; i++) {
                String s = args[i].trim();
                if(s.length() > 0) {
                    //System.out.println(i + " " + s + " added"); 
                    toInclude.add(s);
                }
            }

          
            FileInputStream fis = new FileInputStream(new File(databaseName));
            Iterator<Fasta> fastas = FastaReader.getFastas(fis);
            int i = 0;
            while(fastas.hasNext()) {
                Fasta f = fastas.next();
                String defline = f.getDefline();
                for(Iterator<String> it = toInclude.iterator(); it.hasNext();) {
                    if(defline.indexOf(it.next()) > -1 ) {
                        System.out.println(">" + f.getDefline());
                        System.out.println(f.getSequence());
                        break;
                    }
                }
            }
            fis.close();

        } catch (Exception e) {
            System.err.println(USAGE);
            e.printStackTrace(); 
        }
    }
}
