
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

public class SetOperator {
    private static final String USAGE = "\n\nUsage: java SetOperator file1 file2"; 
    public static void main(String args[]) {
        try {
            String file1 = args[0];
            String file2 = args[1];
            List<String> list1 = getList(file1); 

            List<String> list2 = getList(file2);
            Set<String> union = getUnion(list1, list2);
            Set<String> intersection = getIntersection(list1, list2); 
            Set<String> list1Only = onlyInList1(list1, list2);
            Set<String> list2Only = onlyInList1(list2, list1);
            
            String header = "Number of proteins in the union of " + file1 + " and " + file2 + ": " + union.size();
            output(header, union);
            header = "Number of proteins in the intersection of " + file1 + " and " + file2 + ": " + intersection.size();
            output(header, intersection);
            header = "Number of proteins in " + file1 + " only: " + list1Only.size();
            output(header, list1Only);
            header = "Number of protins in " + file2 + " only: " + list2Only.size();
            output(header, list2Only);

        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
            System.exit(1);
        }
        
        
    }
    public static void output(String header, Set<String> result) {

        System.out.println(header);
        for(String s : result) {
            System.out.println(s);
        }
        System.out.println("\n\n");
    }
    public static Set<String> onlyInList1(List<String> list1, List<String> list2) {
        Set<String> result = new HashSet<String>(list1);
        result.removeAll(list2);
        return result;

    }
    public static Set<String> getIntersection(Set<String> set1, Set<String> set2) {
        Set<String> result = new HashSet<String>(set1);
        result.retainAll(set2);
        return result;
    }
    public static Set<String> getIntersection(List<String> list1, List<String> list2) {
        Set<String> result = new HashSet<String>(list1);
        result.retainAll(list2);
        return result;
    }
    public static Set<String> getUnion(Set<String> list1, Set<String> list2) {
        Set<String> result = new HashSet<String>(list1);
        result.addAll(list2);
        return result;
    }
    public static Set<String> getUnion(List<String> list1, List<String> list2) {
        Set<String> result = new HashSet<String>(list1);
        result.addAll(list2);
        return result;
    }
    private static List<String> getList(String fileName) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        List<String> list = new LinkedList<String>();
        String line = null;
        while((line = br.readLine()) != null) {
            //System.out.println(line.split("\\s")[0]);
            list.add(line.split("\\s")[0]);
        }

        return list;
    }
}
