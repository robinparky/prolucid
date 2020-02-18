/**
 * @file JoinManyFiles.java
 * This is the source file for edu.scripps.pms.util.spectrum.JoinManyFiles
 * @author Tao Xu
 * @date $Date
 */



//import edu.scripps.pms.util.io.*;

//import edu.scripps.pms.util.io.DTASelectFilterReader; 
//import edu.scripps.pms.util.dtaselect.Peptide; 
//import edu.scripps.pms.util.dtaselect.Protein; 
//import edu.scripps.pms.util.io.SpectrumReader; 
import java.io.*;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

// this program expect a header line in all files to be joined
public class JoinManyFiles {
    public static final String USAGE = "java JoinManyFiles fileName1 fileName2 ... key_column_position\n" +
                                       "The columns must be tab delimited";
    private static String header = null;

    public static int processFile(String file, int keyIndex, ArrayList<String> keys, HashMap<String, String> map, StringBuffer header) throws IOException {
//System.out.println("Processing file " + file);
        BufferedReader br = new BufferedReader(new FileReader(file), 4096);
        boolean isHeaderLine = true; 
//System.out.println("Processing file: " + file);
        String line = null;
        int numColumns = 0;
        while((line = br.readLine()) != null) {
//System.out.println("line: " + line);
            if(isHeaderLine) {
                numColumns = line.split("\t").length;
                header.append(line);
                header.append("\t");
                isHeaderLine = false;
            } else {
                String [] arr = line.split("\t");
                if(arr.length < keyIndex || arr.length < 2) continue;
  //              System.out.println(line + "\t" + arr.length);
                String key = arr[keyIndex]; 
                map.put(key, line);  
                keys.add(key);
            }  
        }
//System.out.println("\tNumKeys: " + keys.size() + "\tmap size: " + map.size());
        br.close();
        return numColumns;
    }

    public static String joinFiles(ArrayList<String> files, int keyindex) throws IOException {
        StringBuffer sb = new StringBuffer(1000000);
        StringBuffer header = new StringBuffer(1000);

        HashSet<String> allkeys = new HashSet<String>(1000000);
        int numfiles = files.size();

        ArrayList<HashMap<String, String>> maps = new ArrayList<HashMap<String, String>>();
        ArrayList<ArrayList<String>> keys = new ArrayList<ArrayList<String>>();
        int numColumns [] = new int[numfiles];
        String emptyColumns [] = new String[numfiles];
        
        for(int i = 0; i < files.size(); i++) {
            HashMap<String, String> map = new HashMap<String, String>(100000);
            maps.add(map);
            ArrayList<String> key = new ArrayList<String>(100000);
            keys.add(key);

            //System.out.println("Number of keys: " + keys.size());
            numColumns[i] = processFile(files.get(i), keyindex, key, map, header); 
            allkeys.addAll(key);

            StringBuffer tabs  = new StringBuffer();
            for(int j = 0; j < numColumns[i]-1; j++) {
                tabs.append("\t");
            } 
            emptyColumns[i] = tabs.toString();  
        }

        header.append("sorting_index");
        sb.append(header.toString());
        sb.append("\n");
//System.out.println("Number of all keys: " + allkeys.size());
        for(Iterator<String> it = allkeys.iterator(); it.hasNext();) {
            String k = it.next();
//System.out.println("processing key: " + k);
            StringBuffer status = new StringBuffer();
            for(int i = 0; i < maps.size(); i++) {
                HashMap<String, String> map = maps.get(i);
                String s = map.get(k);
                if(s == null) { 
                    s = emptyColumns[i];
                    status.append("0_");
                } else {
                    status.append("1_");
                }
                sb.append(s); sb.append("\t");
            } 
            sb.append(status.toString());
            sb.append("\n");
        }

        return sb.toString();
    }
    public static void main(String args[]) throws Exception {
         
        try {

            int keyindex = Integer.parseInt(args[args.length-1]) - 1;
            ArrayList<String> files = new ArrayList<String>();
            for(int i = 0; i < args.length - 1; i++) {
                files.add(args[i]);
            }
            System.out.println(joinFiles(files, keyindex).toString()); 
        } catch (Exception e) {
            System.out.println(USAGE);
            e.printStackTrace();
            System.exit(1);
            System.out.println(USAGE);
        }
    }
}
