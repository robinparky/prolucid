
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

import edu.scripps.pms.pmsdb.*;

import ubic.db.*;
import java.sql.*;
import java.util.Date;
import java.util.*;
import java.io.*;

public class Accessions2Deflines {
    private static final String dbName = "pmsdb";
    public static final String USAGE = "\n\n--- java Accessions2Deflines accessionFileName ---\n";
    public static void main(String args[]) {
        try {
        Pmsdb pdb = new Pmsdb(dbName);
        String fileName = args[0];
        String outputFile = fileName+".out";
         
        PrintWriter output = new PrintWriter(outputFile);
        List<String> acs = readAccessions(fileName);
        for(String ac : acs) {
            String defline = pdb.accession2Defline(ac);
            String [] contents = defline.split(" ");
            output.print(ac + "\t");
            for(int i = 2; i < contents.length; i++) {
                output.print(contents[i] + " "); 
            }
            
            output.println();
        }
        output.close();
        System.out.println("Please find the output in " + outputFile + " file");
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
        }
    }
    
    private static List<String> readAccessions(String fileName) throws IOException {
        ArrayList<String> acs = new ArrayList<String>();
        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
        String line = null;
        while((line = br.readLine()) != null) {
            if(line.equals("") || line.indexOf(" ") == -1) {
                System.out.println(line);
                continue;
            }
            String [] contents = line.split(" ");
            acs.add(contents[0].trim()); 
        }
        br.close();
        return acs;
    }
}
