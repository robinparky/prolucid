/**
 * @file CamsiId.java
 * This is the source file for CamsiId
 * @author Tao Xu
 * @date $Date
 */



import edu.scripps.pms.util.io.*;

//import edu.scripps.pms.util.io.DTASelectFilterReader; 
import edu.scripps.pms.util.dtaselect.Peptide; 
import edu.scripps.pms.util.dtaselect.Protein; 
import edu.scripps.pms.util.PmsUtil; 
import edu.scripps.pms.util.stat.StatCalc; 
import edu.scripps.pms.util.sqt.SQTPeptide;
import edu.scripps.pms.util.sqt.MLine;


import java.io.*;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import java.util.HashSet;
import edu.scripps.pms.util.spectrum.*;

// command line processor
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

// this is the class for generate CAMSI repotring files. Modified from DeltaMassZScoreCalculator.java file
public class CamsiId implements Comparable<CamsiId> {
    private int scan;
    private HashSet<String> proteinids = new HashSet<String>(5);
    private String peptidesequence;
    private double score;
    private double fpr; // for false positive rate
    private String annotation;
    private boolean isModified = false; 
    private String contentline;

    public CamsiId(String line, int scoreindex) {
        contentline = line;
        String [] arr = line.split("\\t");
        scan = Integer.parseInt(arr[0]);
        String proteins = arr[1];
        
        String [] acs = null;
        if(proteins.indexOf("%") > -1) {
            acs = proteins.split("%");
        } else if(proteins.indexOf(";") > -1) {
            acs = proteins.split(";");
        } else if(proteins.indexOf(",") > -1) {
            acs = proteins.split(",");
        } else if(proteins.indexOf(" ") > -1) {
            acs = proteins.split(" ");
        } else {
            acs = proteins.split(":");
        } 
        for(int i = 0; i < acs.length; i++) {
            if(acs[i] != null && !acs[i].equals("")) {
                proteinids.add(acs[i].trim());
            }
        }
        peptidesequence = arr[2];
        score = Double.parseDouble(arr[scoreindex]);
        //fpr = Double.parseDouble(arr[4]);
        if(arr.length > 5) {
            annotation = arr[5];
        }
    }
    public String getContentLine() {
        return contentline;
    }
    public int compareTo(CamsiId id) {
        if(id.score < this.score) {
            return 1;
        } else if(id.score > this.score) {
            return -1;
        } else {
            return 0;
        }
    }
    public void isModified(boolean modified) {
        isModified = isModified || modified;
    }
    public CamsiId(int scan) {
        this.scan = scan;
    }   
    public Iterator<String> getProteinIds() {
         return proteinids.iterator();
    }
    public void addProteinId(String acc) {
        proteinids.add(acc);
    } 
    public void setPeptideSequence(String seq) {
        peptidesequence = seq;
    }
    public void setScore(double s) {
        score = s;
    }
    public void setFalsePositiveRate(double f) {
        fpr = f;
    }
    public int getScan() {
        return scan;
    }
    public double getFalsePositiveRate() {
        return fpr;
    }
    public double getScore() {
        return score;
    } 
    public String getPeptideSequence() {
        return peptidesequence;
    }
    public void resetProteinIds() {
        proteinids = new HashSet<String>(10);
    }
    public String output() {
        StringBuffer sb = new StringBuffer(200);
        sb.append(scan);
        sb.append("\t");
        Iterator<String> it = proteinids.iterator();
        if(it.hasNext()) {
            sb.append(it.next());
        }
        while(it.hasNext()) {
            sb.append(";" + it.next());
        }
        sb.append("\t");
     
        sb.append(peptidesequence + "\t");
        sb.append(score + "\t");
        sb.append(fpr);
        if(isModified) {
            sb.append("\tmodified");
        }
        if(annotation != null) {
            sb.append("\t" + annotation);
        }

        return sb.toString(); 
 
    }
}
