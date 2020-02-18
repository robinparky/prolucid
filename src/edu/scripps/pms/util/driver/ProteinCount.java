
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
import edu.scripps.pms.util.io.*;
import java.util.*;
import java.io.*;

public class ProteinCount implements Comparable<ProteinCount> {
    private String accession;
    private double forwardCount = 0;
    private double reverseCount = 0;
    
    public ProteinCount(String accession) {
        this.accession = accession;
    }
    public void addCount(String acc, int rank) {
//System.out.println("acc: " + acc + "\taccession: " + accession);
        if(rank == 1) {
            if(acc.startsWith("Reverse_")) {
                reverseCount++;
            } else {
                forwardCount++;
            }
        } else {

            if(acc.startsWith("Reverse_")) {
                reverseCount += 0.5;
            } else {
                forwardCount += 0.5;
            }
        } 
 
    }
    public double getReverseCount() {
        return reverseCount;
    }
    public double getForwardCount() {
        return forwardCount;
    }
    public double getTotalCounts() {
        return reverseCount + forwardCount;
    }
    public String getAccession() {
        return accession;
    }
    public double getCountDifference() {
        return forwardCount - reverseCount;
    }
    public int compareTo(ProteinCount p) {
        double result = (p.getCountDifference()) - (getCountDifference());
        if(result > 0) {
            return 1;
        } else if (result < 0) {
            return -1;
        }
        return 0;
    }
}
