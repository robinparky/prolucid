/**
 * @file Accession2Fastas.java
 * This is the source file for Accession2Fastas class
 * @author Tao Xu
 * @date $Date
 */

import edu.scripps.pms.mspid.ProteinDatabase;

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


public class Accessions2Fastas {

    public static final String USAGE = "!!! Usage: accession2fastas accessionNumberFile databaseFile !!!";
    // databaseName - the path and name of the protein database

    private static HashSet<String> readAccessions(String accfile) throws IOException {
        
        HashSet<String> accessions = new HashSet<String>(1000000);
        BufferedReader br = new BufferedReader(new InputStreamReader( new FileInputStream(accfile)));    
        String ac = br.readLine();
        while ( ac != null) {
            ac = ac.trim();
            accessions.add(ac);
            ac = br.readLine(); 
        }
        br.close();
        return accessions;
    }
    // assuming the mass tolerance is smaller than the mass of any AA residue
    public static void main(String [] args) throws Exception {
        try {
            String accfile = args[0];
            String databaseName = args[1];
            HashSet accs = readAccessions(accfile);

          
            FileInputStream fis = new FileInputStream(new File(databaseName));
            Iterator<Fasta> fastas = FastaReader.getFastas(fis);
            int i = 0;
            while(fastas.hasNext()) {
                Fasta f = fastas.next();
                String myac = f.getAccession();
                String myacnoversion = f.getAccessionWithNoVersion();
                if(accs.contains(myac) || accs.contains(myacnoversion)) {
                    System.out.println(">" + f.getDefline());
                    System.out.println(f.getSequence());

                }
            }
            fis.close();

        } catch (Exception e) {
            System.err.println(USAGE);
            e.printStackTrace(); 
        }
    }
}
