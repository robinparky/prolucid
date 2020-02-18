/**
 * @file DecoyEntryRemover.java
 * This is the source file for DecoyEntryRemover class
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


public class DecoyEntryRemover {

    public static final String USAGE = "!!! Usage: decoyntryremover fasta_file Decoy_Entry_Prefix !!!";
    // databaseName - the path and name of the protein database

    // assuming the mass tolerance is smaller than the mass of any AA residue

    public static void removeDecoyEntries(String databaseName, String decoyentryprefix) throws IOException {

        FileInputStream fis = new FileInputStream(new File(databaseName));
        Iterator<Fasta> fastas = FastaReader.getFastas(fis);
        int i = 0;
        while(fastas.hasNext()) {
            Fasta f = fastas.next();
            String defline = f.getDefline();
            if(!defline.startsWith(decoyentryprefix)) {
                System.out.println(">" + f.getDefline());
                System.out.println(f.getSequence());
            } else {
                //System.out.println("Decoy entry " + f.getDefline());
            } 
        }
        fis.close();
    }
    public static void main(String [] args) throws Exception {
        try {
            String databaseName = args[0];
            String decoyentryprefix = args[1];
            removeDecoyEntries(databaseName, decoyentryprefix);
//System.out.println("Decoy prefix: " + decoyentryprefix)Accession2Fastas.java;
          

        } catch (Exception e) {
            System.err.println(USAGE);
            //e.printStackTrace(); 
        }
    }
}
