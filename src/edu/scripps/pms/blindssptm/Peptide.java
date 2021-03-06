/**
 * @file Peptide.java
 * This is the source file for edu.scripps.pms.util.seq.Peptide
 * @author Tao Xu
 * @date $Date
 */
package edu.scripps.pms.blindssptm;



import edu.scripps.pms.util.seq.Fasta;

public class Peptide {

    // description line of this Peptide sequence
    // the sequence string of this Peptide
    
    protected int start;
    protected int end;
    protected Fasta parent;
    //private int chargeState;

    public Peptide(Fasta parent, int start, int end) {
        this.parent = parent;
        this.start = start;
        this.end = end;
    }  

    public byte byteAt(int index) {
        return parent.byteAt(index + start);
    }
    public boolean equals(Peptide p) {

        if(p != null) {
            return parent == p.parent && start == p.start && end == p.end; 
        }
        return false;
    }
    public Fasta getParent() {
        return parent;
    }
    public String getExtendedSequence() {  

        // 4 more bytes for the heading and tailing aa and '.'
        int len = getLength() + 4;
        byte [] seq = new byte[len];
        if (start == 0) {
            seq[0] = '-';
        } else {
            seq[0] = parent.byteAt(start-1);
        }
        seq[1] = '.';

        int lastIndex = parent.getLength() - 1;
        if (end == lastIndex) {
            seq[len-1] = '-';
        } else {
            seq[len-1] = parent.byteAt(end+1);
        }
        seq[len-2] = '.';
        int index = 2;
        for (int i = start; i < end+1; i++) {
            seq[index++] = parent.byteAt(i);
        }

        return new String(seq);        
    }
    public String getSequence() {
        return new String(parent.getSequenceAsBytes(), start, getLength());
    }
    public int getStart() {
        return start;
    }
    public String getDefline() {
        return parent.getDefline();
    }
    public int getEnd() {
        return end;
    }

    public int getLength() {
        return end - start + 1;
    }
    public String getSequestLikeAccession() {
        return parent.getSequestLikeAccession();
    }
    public String getAccession() {
        return parent.getAccession();
    }
}
