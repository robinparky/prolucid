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


public class Fasta3PrimeUtrSplitter {

    public static final String USAGE = "!!! Usage: fasta3primeutrsplitter databaseFile [num_residues_for_3_prime_utr] !!!";
    // databaseName - the path and name of the protein database

    // assuming the mass tolerance is smaller than the mass of any AA residue
    public static void main(String [] args) throws Exception {
        try {
            String databaseName = args[0];
            int utrlength = 300;
            if(args.length > 0) {
                utrlength = Integer.parseInt(args[1]);
            }

          
            FileInputStream fis = new FileInputStream(new File(databaseName));
            Iterator<Fasta> fastas = FastaReader.getFastas(fis);
            int i = 0;
            int numtooshort = 0;
            while(fastas.hasNext()) {
                Fasta f = fastas.next();
                String defline = f.getDefline();
                String acc = f.getAccession();

                int pipeindex = acc.indexOf("|");
                if(pipeindex > 0) {
                    acc = acc.substring(0, pipeindex); 
                }


                String utracc = "3UTR_" + acc;
                String orgseq = f.getSequence();
                int orglen = orgseq.length();
                if(orglen > 300) {
                    String cdsseq = orgseq.substring(0, orglen - utrlength);
                    String utrseq = orgseq.substring(orglen - utrlength, orglen);
    //System.out.println(acc + "\toriginal length: " + orglen + "\tcds length: " + cdsseq.length() + "\tutr length: " + utrseq.length());
    System.out.println(">" +  acc + "\t" + f.getDefline());
    System.out.println(cdsseq);
    System.out.println(">" + utracc + "\t3' UTR of " + f.getDescription());
    System.out.println(utrseq);
                } else {
                    System.out.println(">" + acc + "\t" +  f.getDefline());
                    System.out.println(orgseq);

                    numtooshort++;
                }

            }
            fis.close();
            //System.out.println("Number too short: " + numtooshort);

        } catch (Exception e) {
            System.err.println(USAGE);
            e.printStackTrace(); 
        }
    }
}
