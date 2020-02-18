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


public class FastaAccessionSwitch {

    public static final String USAGE = "!!! Usage: fastaaccessionswitch databaseFile !!!";
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
            
            HashSet<String> outputseq = new HashSet(1000000);
 
            FileInputStream fis = new FileInputStream(new File(databaseName));
            Iterator<Fasta> fastas = FastaReader.getFastas(fis);
            int numredudant = 0;
            while(fastas.hasNext()) {
                Fasta f = fastas.next();
                String defline = f.getDefline();
                String seq = f.getSequence();
                if(outputseq.add(seq)) {
                    String [] arr = defline.split("\\s+");
                    System.out.println(">" +  arr[1] + "\t" + f.getDefline());
                    System.out.println(seq);
                } else {

                    numredudant++;
                }
            }
           // System.out.println("Number of redundant enrties removed: " + numredudant); 
            fis.close();

        } catch (Exception e) {
            System.err.println(USAGE);
            e.printStackTrace(); 
        }
    }
}
