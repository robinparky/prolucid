/**
 * @file ModifiedPeptideHit.java
 * This is the source file for edu.scripps.pms.util.seq.ModifiedPeptideHit
 * @author Tao Xu
 * @date $Date: 2007/04/12 22:54:13 $
 */



package edu.scripps.pms.topdown;

import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.mspid.ProcessedPeakList;

import java.util.*;

public class ModifiedPeptideHit extends PeptideHit {
     
    protected Modifications mods;
    protected String extendedSeq = null;
    protected DiffMod [] diffMods;
    protected ArrayList<ModifiedPeptideHit> candidates;
    protected ArrayList<DiffModSite> sites = new ArrayList<DiffModSite>();
    protected int numPotentialSites = -1;
    protected double internalDiffModMassAdded = 0;
    protected int numInternalDiffModsAdded = 0;
    public static final boolean isModified = true;


    public ModifiedPeptideHit(Fasta parent, int start, int end, Modifications m) {
        super(parent, start, end);
        this.mods = m;
        diffMods = new DiffMod[end-start+1];
    }  
    public ModifiedPeptideHit(Peptide p, Modifications m) {
        super(p);
        mods = m;
        diffMods = new DiffMod[p.getEnd() - p.getStart() + 1];

    }
    public boolean isModified() {
        return isModified;
    }
    public void setDiffMod(int index, DiffMod diffmod) {
        diffMods[index] = diffmod;
        if(diffmod != null) {
            internalDiffModMassAdded += diffmod.getMassShift();
            numInternalDiffModsAdded++;
        }
    }
    public boolean isGoodCandidate() {
        return internalDiffModMassAdded == mods.getInternalDiffModsMass() &&
              numInternalDiffModsAdded == mods.getNumDiffMods();
    }
    public void setDiffMods(DiffMod [] dms) {
        this.diffMods = dms;
    }
    protected void generateCandidates() {
        candidates = new ArrayList<ModifiedPeptideHit>();
        calcModificationSites();
//System.out.println("numSites: " + numPotentialSites + "\t numDiffMods: " + mods.getNumDiffMods());
        if(getNumPotentialSites() < mods.getNumDiffMods()) {
            return;
        }
        candidates.add(this.copy()); 
        int numSites = sites.size(); 

//System.out.print("numPotentialSitesInPeptide: " + numPotentialSites +  "\t numDiffMods: " + mods.getNumDiffMods());
        // may have more efficient way if numSites == mods.getNumDiffMods()
        for(int i = 0; i < numSites; i++) {
           
            DiffModSite dms = sites.get(i); 
            processDiffModSite(dms);
        }
        getFinalCandidates();
    }
    private void getFinalCandidates() {
        
        ArrayList<ModifiedPeptideHit> temp = candidates;
        candidates = new ArrayList<ModifiedPeptideHit>();
        for(Iterator<ModifiedPeptideHit> it = temp.iterator(); it.hasNext();) {
            ModifiedPeptideHit m = it.next();
            if(m.isGoodCandidate()) {
                candidates.add(m);
            }
        } 
//System.out.println("\tNumTotal Candidates: " + temp.size() + "\tNumGoodCandidates: " + candidates.size());
    }
    private void processDiffModSite(DiffModSite dms) {
//System.out.println("numDiffModsForThisSite: " + dms.numDiffMods());
        int num = candidates.size();
        int position = dms.getPosition();
        for(int i = 0; i < num; i++) {
            ModifiedPeptideHit m = candidates.get(i);
            for(Iterator<DiffMod> it = dms.getDiffMods(); it.hasNext();) {
                DiffMod d = it.next();
                if(m.numInternalDiffModsAdded+1 <= mods.getNumDiffMods() &&
                          m.internalDiffModMassAdded+d.getMassShift() <= mods.getInternalDiffModsMass()) {
                    ModifiedPeptideHit newm = m.copy();
                    newm.setDiffMod(position, d);
                    candidates.add(newm); 
                }
            }
        }
    }   
    private ModifiedPeptideHit copy() {
        ModifiedPeptideHit newm = new ModifiedPeptideHit(peptide, mods);
        for(int i = 0; i < newm.diffMods.length; i++) {
            newm.setDiffMod(i, diffMods[i]);
        }              
        //newm.internalDiffModMassAdded = internalDiffModMassAdded;
        return newm;
    }
    protected void calcModificationSites() {
        if(numPotentialSites > -1) {
            return;
        }
        numPotentialSites = 0;
        int start = getStart();
        int end = getEnd();
        byte [] bytes = getParent().getSequenceAsBytes();
        for(int i = start; i <= end; i++) {
            if(mods.isModifiable(bytes[i])) {
                numPotentialSites++;
                sites.add(new DiffModSite(i - start, mods.residue2DiffMods(bytes[i])));
            }
        }
    }
    public int getNumPotentialSites() {
        if(numPotentialSites > -1) {
            calcModificationSites();
        }
        return numPotentialSites; 
    }    

    public Iterator<ModifiedPeptideHit> getAllModifiedPeptideHits() {
        generateCandidates(); 
        return candidates.iterator();
    }
    public double getTheorMass() {
        return ppl.getMassCalculator().getPrecursorMass(peptide.getSequence()) +  mods.getMassShift();
    }
    // override the function in super class
    public double getBStart() { // b == n term

        return ppl.getBStart() + mods.getDiffNTermMod();
    }    
    // override the function in super class
    public double getDbinwidthBStart() { // b == n term

        return (ppl.getBStart() + mods.getDiffNTermMod())*MassSpecConstants.DBINWIDTH;
    }    
    public double getYStart() { // y == c term
//System.out.println("in getYStart in ModifiedPeptideHit, " + mods.getDiffCTermMod());
        return ppl.getYStart() + mods.getDiffCTermMod();
    }    
    public double getDbinwidthYStart() { // y == c term
        return (ppl.getYStart() + mods.getDiffCTermMod())*MassSpecConstants.DBINWIDTH;
    }    
    public String getExtendedSequence() {
        if(extendedSeq != null) {
            return extendedSeq;
        }
        StringBuffer sb = new StringBuffer(getLength() + 15);
        // 4 more bytes for the heading and tailing aa and '.'
        int len = peptide.getLength() + 4;
        int start = peptide.getStart();
        int end = peptide.getEnd();
        if (start == 0) {
            sb.append('-');
        } else {
            sb.append((char)peptide.getParent().byteAt(start-1));
        }
        sb.append('.');

        double nTermMod = mods.getDiffNTermMod();
        double cTermMod = mods.getDiffCTermMod();
        if(nTermMod != 0) {
            sb.append("(" + nTermMod + ")"); 
        }

       // sb.append(getSequence()); // need to be modified for regular mod search
        for(int i = start; i <= end; i++) {
            sb.append((char)peptide.getParent().byteAt(i));
            DiffMod d = diffMods[i-start];
            if(d != null) {
                sb.append(d.getModInfo()); 
            }
        }
       
        if(cTermMod != 0) {
            sb.append("(" + cTermMod + ")"); 
        }
        sb.append('.');

        int lastIndex = peptide.getParent().getLength() - 1;
        if (end == lastIndex) {
            sb.append('-');
        } else {
            sb.append((char)peptide.getParent().byteAt(end+1));
        }

        extendedSeq = sb.toString(); 
        
        return extendedSeq;
    }


    public void getCharge3NumPeaksMatched() {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        byte [] seq = peptide.getParent().getSequenceAsBytes();
        //fragMasses = mc.getFragMasses(); no big difference

        double yMass = getYStart(); // overriden function
        double bMass = getBStart(); // overriden function
        int yIndex = peptide.getStart() + peptide.getEnd(); // index for y ion, i is used for b ion
        

        MassCalculator mc = ppl.getMassCalculator();
        boolean [] boolMasses = ppl.getBoolMasses();
        int start = peptide.getStart();
        int end = peptide.getEnd();
        for(int i = start; i < end; i++) {
            // for diffmods
            int bmi = i - start; // b diffmod index
            int ymi = yIndex - i - start; // y diffmod index  
            DiffMod bd = diffMods[bmi];
            DiffMod yd = diffMods[ymi]; 
            bMass += bd == null? 0 : bd.getMassShift();
            yMass += yd == null? 0 : yd.getMassShift();
            bMass += mc.getFragmentMass(seq[i]);
            yMass += mc.getFragmentMass(seq[yIndex-i]);

            int int2B = (int)(bMass*ppl.PRECISIONFACTOR + 0.5f);
            int int2Y = (int)(yMass*ppl.PRECISIONFACTOR + 0.5f);
            int int3B = (int) (0.5f + ppl.PRECISIONFACTOR*(bMass + MassSpecConstants.MASSH)/2f);
            int int3Y = (int)(0.5f + ppl.PRECISIONFACTOR*(yMass + MassSpecConstants.MASSH)/2f);
            // doubly charged fragment ions
            if(boolMasses[int3B]) numPeaksMatched++;
            if(boolMasses[int3Y]) numPeaksMatched++;

            // singly charged fragment ions
            if(boolMasses[int2B++]) numPeaksMatched++;
            if(boolMasses[int2Y++]) numPeaksMatched++;

            if(boolMasses[int2B]) numPeaksMatched++;
            if(boolMasses[int2Y]) numPeaksMatched++;
           // if(boolMasses[int2B++]) numPeaksMatched++;
           // if(boolMasses[int2Y++]) numPeaksMatched++;
            numTheroticPeaks += 6;
        }
        setNumPeaksMatched(ppl.getPTrue());

    }

    public void getCharge2NumPeaksMatched() {
        int numPeaksMatched = 0;
        int numTheroticPeaks = 0;
        byte [] seq = getParent().getSequenceAsBytes();
        boolean [] boolMasses = ppl.getBoolMasses();
        MassCalculator mc = ppl.getMassCalculator();
        //fragMasses = mc.getFragMasses(); no big difference
        double yMass = getYStart(); // from super class
        double bMass = getBStart();  // from super class

        int start = peptide.getStart();
        int end = peptide.getEnd();
        int yIndex = start + end; // index for y ion, i is used for b ion

        for(int i = start; i < end; i++) {
            // for diffmods
            int bmi = i - start; // b diffmod index
            int ymi = yIndex - i - start; // y diffmod index  
            DiffMod bd = diffMods[bmi];
            DiffMod yd = diffMods[ymi]; 
            bMass += bd == null? 0 : bd.getMassShift();
            yMass += yd == null? 0 : yd.getMassShift();

            bMass += mc.getFragmentMass(seq[i]);

            yMass += mc.getFragmentMass(seq[yIndex-i]);

            int tempB = (int)(bMass * ppl.PRECISIONFACTOR+ 0.5f);
            int tempY = (int)(yMass*ppl.PRECISIONFACTOR + 0.5f);
            //numPeaksMatched += boolMasses[tempB]? 1 : 0;
            //numPeaksMatched += boolMasses[tempY]? 1 : 0;
            if(boolMasses[tempB++]) numPeaksMatched++;
            if(boolMasses[tempY++]) numPeaksMatched++;

            if(boolMasses[tempB]) numPeaksMatched++;
            if(boolMasses[tempY]) numPeaksMatched++;

            numTheroticPeaks += 4;
        }
        setNumPeaksMatched(ppl.getPTrue());
    }

    public int [] getCharge2TheorMasses() {
        int theorMass[] = new int[(int)ppl.getPrecursorMass() + 100];
        //int theorMass[] = new int[(int)prcMass + 400];
        byte [] seq = getSequence().getBytes();
        double bMass = getDbinwidthBStart() + 0.5f;
        double yMass = getDbinwidthYStart() + 0.5f;

        int numAA = seq.length - 1;
        int yIndex = numAA; // index for y ion, i is used for b ion
        //int factor = ppl.PRECISIONFACTOR;
        for(int i = 0; i < numAA; i++) {
            int yi = yIndex -i;
            bMass += masses[seq[i]];
            yMass += masses[seq[yi]];

            // specific for diff mod
            DiffMod bm = diffMods[i];
            bMass += bm == null? 0 : bm.getDbinwidthMassShift();
            DiffMod ym = diffMods[yi];
            yMass += ym == null? 0 : ym.getDbinwidthMassShift();

            processBIon(theorMass, bMass);
            processYIon(theorMass, yMass);
        }
        return theorMass;
    }

    public boolean equals(ModifiedPeptideHit m) {
        return super.equals(m) && mods.equals(m.mods);
    }
}