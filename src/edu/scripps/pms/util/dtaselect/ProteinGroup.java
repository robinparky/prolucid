package edu.scripps.pms.util.dtaselect;

/**
 * @author Tao Xu 
 * @version $Id: ProteinGroup.java,v 1.2 2014/07/08 23:52:13 rpark Exp $
 */
import java.util.*;
import edu.scripps.pms.util.seq.Fasta;

public class ProteinGroup {
    private ArrayList<Protein> proteins = new ArrayList();
    private HashSet<String> proteinAccessions = new HashSet();
    private List<Peptide> peptides = new ArrayList(); 
 
    public void addProtein(Protein p) {
        proteins.add(p);
        if(p.getNumPeptides() > 0) {
            peptides = p.getPeptideList();
        } 
    }

    public Iterator<Protein> getProteins() {
        return proteins.iterator();
    }
        
    public Iterator<Peptide> getPeptides() {
        return peptides.iterator();
    }
    public int getNumPeptides() {
        return peptides.size();
    } 
    public int getNumProteins() {
        return proteins.size();
    } 
    public boolean containsProtein(String ac) {
        return proteinAccessions.contains(ac);
    }
    public Protein ac2Protein(String ac) {
        for(int i = 0; i < proteins.size(); i++) {
            Protein p = proteins.get(i);
            if(p.getAccession().equals(ac)) {
                return p;
            }
            
        }
        return null;
    }


}
