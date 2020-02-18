/**
 * @file DiffModSite.java
 * This is the source file for edu.scripps.pms.blindptm.Modification
 * @author Tao Xu
 * @date $Date: 2007/01/20 00:52:13 $
 */

package edu.scripps.pms.blindptm;


import java.util.*;


public class DiffModSite {
//    private ArrayList<Double> differentialModifications = new ArrayList<Double>();
    private int position;
    private HashSet<DiffMod> mods;

    public DiffModSite(int pos, HashSet<DiffMod> mods) {
        position = pos;
        this.mods = mods;
    }
    public int numDiffMods() {
        return mods.size();
    }
    public Iterator<DiffMod> getDiffMods() {
        return mods.iterator();
    }
    public int getPosition() {
        return position;
    }
}



