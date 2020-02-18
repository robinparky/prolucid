
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
import java.lang.management.*;
import java.io.*;

public class TestDriver1{
   // System.out.println("haha");
    //int i = 5;
    public static void main(String args[]) {
        
        int numcpu = ManagementFactory.getOperatingSystemMXBean().getAvailableProcessors();
        System.out.println("number of processors: " + numcpu);

        MemoryMXBean mbean = ManagementFactory.getMemoryMXBean(); 
        System.out.println("heap size: " + mbean.getHeapMemoryUsage());

        String path = "/profuse/blahblah";
        int index = path.lastIndexOf("/");
        int index2 = path.lastIndexOf("\\/");
        
        System.out.println("index: " + index + "\tindex2: " + index2);
        String sub = path.substring(0, index);
        //String sub2 = path.substring(0, index2);
        String sub2 = path.substring(0, path.length());
        
        System.out.println("index: " + index + "\tsub: " + sub + "\tindex2: " + index2 + "\tsub2: " + sub2);


        String s1 = "grapefruit";
        String s2 = "grapefruit";


        String s3 = "grape"+"fruit"; 

        String s4 = new String("grapefruit");
        String s5 = new String("grapefruit");


       System.out.println("memory address of s1 to s5");
       System.out.println(mimicObjectToString(s1));
       System.out.println(mimicObjectToString(s2));
       System.out.println(mimicObjectToString(s3));
       System.out.println(mimicObjectToString(s4));
       System.out.println(mimicObjectToString(s5));

      
       s1 = "23";
       s2 = "0.3";
       s3 = ".3";
   
       System.out.println(s1.split("\\.")[0]);
       System.out.println(s2.split("\\.")[0]);
       System.out.println(s3.split("\\.")[0]);

       System.out.println(s1 + "\tindex of dot is: " + s1.indexOf("."));
       s1 = s1.split("\\.")[0];
       System.out.println(s1);

       System.out.println("Now checking file size");
       File f = new File("/lustre/people/applications/yates/dbase/EBI-IPI_human_3.24_12-01-2006_con_reversed.fasta");

       System.out.println("File size is: " + f.length());
    }   

    public static String mimicObjectToString(Object o)   
        {   
            //prevent a NullPointerException by returning null if o is null   
            String result = null;   
            if (o !=null)   
            {   
                result = o.getClass().getName() + "@" + Integer.toHexString(System.identityHashCode(o));   
            }   
            return  result;   
        }   

}

