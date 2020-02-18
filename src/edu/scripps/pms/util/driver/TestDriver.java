
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
import java.util.*;
import java.io.*;

public class TestDriver{
   // System.out.println("haha");
    //int i = 5;
    public static void main(String args[]) throws IOException {
      
        
        String f1 = "1.23456789";
        String f2 = "12345678901112";
        String f3 = "9999.123456789";
        System.out.println(f1 + "\t" + Float.parseFloat(f1));
        System.out.println(f2 + "\t" + Float.parseFloat(f2));
        System.out.println(f3 + "\t" + Float.parseFloat(f3));
        for(int i = 1; i < 100; i++) {
            double fact = 1;
            for(int j = i; j > 1; j--) {
                fact *= j;
            }
            System.out.println(i + "\t" + fact);
        }
        String s = "/data/1/root/";
        System.out.println("Index of / is: " + s.indexOf("/")); 
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        //System.err.println("Please input an integer for N:");
        //int n = Integer.parseInt(br.readLine());
        while(true) {
            System.err.println("Type exit or quit to exit, or enter to continure");
            String quit = br.readLine();
            if(quit == null || quit.startsWith("q") || quit.startsWith("e")) { 
               // System.exit(0);
               break;
            }
            System.err.print("Please input number of years: ");
            int numYears = Integer.parseInt(br.readLine());
            System.err.print("Please input original amount: ");
            double original = Double.parseDouble(br.readLine());
            System.err.print("Please input addition amount: ");
            double addition = Double.parseDouble(br.readLine());
            System.err.print("Please input rate: ");
            double rate = Double.parseDouble(br.readLine());
            double balance = original; 
            System.err.print("Please input inflation rate: ");
            double inflationRate = Double.parseDouble(br.readLine());
            System.out.println("Number of Years: " + numYears + "\tStart with: " + original + "\tAddition per year: " + addition + "\tAPR: " +  rate + "\tInflationRate: " + inflationRate);
            for(int i = 0; i < numYears; i++) {
                addition *= inflationRate;
                double gain = balance * (rate-1);
                balance = balance * rate + addition; 
                System.out.println("Year " + (i + 1) + ": " + balance + "\t" + balance/Math.pow(inflationRate, i+1) + "\t" + addition + "\t" + gain);
            } 
            for(int i = 0; i < numYears; i++) {
                balance /= inflationRate;
            }
            System.out.println("Current value: " + balance);
            HashMap<Integer, String> hm = new HashMap<Integer, String>();   
            hm.put(new Integer(8000000), "hahaha");
            System.out.println(hm.get(new Integer(8000000)));
        }
    }   

}

