
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
import edu.scripps.pms.util.spectrum.*;
import java.util.*;
import java.io.*;

public class Randomizer {
    private String input;
    private ArrayList<String> toBeRandomized = new ArrayList<String>(1000000);
    private int numEntries;

    public Randomizer(String inputfile) throws IOException {
        input = inputfile;
        readInput();
        numEntries = toBeRandomized.size();
    }
    private void readInput() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(input));
        String s = br.readLine();
        while(s != null) {
            toBeRandomized.add(s);
            s = br.readLine();
        }
        br.close();
    }
    public void randomizedOutput() {
        Random rand = new Random(); 
        int max = numEntries;
        while (max > 0) {
            int index = rand.nextInt(max);
            String line = toBeRandomized.get(index);
            if(line != null) {
                System.out.println(line);
                toBeRandomized.remove(index);
                max--;
            } 
        }
    }
    public static void main(String [] args) {
        try {
            Randomizer r = new Randomizer(args[0]);
            r.randomizedOutput();
            
        } catch (Exception e) {
            e.printStackTrace();
        } 
    }
}
