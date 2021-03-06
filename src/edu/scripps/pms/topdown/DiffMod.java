/**
 * @file DiffMod.java
 * This is the source file for edu.scripps.pms.mspid.Modification
 * @author Tao Xu
 * @date $Date: 2007/04/12 22:49:56 $
 */

package edu.scripps.pms.topdown;


import java.util.*;


public class DiffMod {
    private double massShift;
    private char symbol;
    //private HashSet<Modification> mods = new HashSet<Modification>();
    int hashcode;
    private String info;
    private boolean [] modifiables = new boolean[256];
    private double dbinwidthMassShift;
    private double neutralLoss;
    private double dbinwidthNeutralLoss;
    
//    private ArrayList<Double> differentialModifications = new ArrayList<Double>();

    public DiffMod(double massShift, char symbol) {
        this.symbol = symbol;
        this.massShift = massShift;
        dbinwidthMassShift = massShift*MassSpecConstants.DBINWIDTH;
        hashcode = new Double(massShift).hashCode();
        info = "(" + massShift + ")";
    }
    public DiffMod(double massShift, double neutralLoss, char symbol) {
        this.symbol = symbol;
        this.massShift = massShift;
    
        dbinwidthMassShift = massShift*MassSpecConstants.DBINWIDTH;
        hashcode = new Double(massShift).hashCode();
        this.neutralLoss = neutralLoss;  
        dbinwidthNeutralLoss = neutralLoss*MassSpecConstants.DBINWIDTH;
        info = "(" + massShift + ")";
    }
    public void setModifiable(byte residue, boolean modifiable) {
        modifiables[residue] = modifiable;
    }
    public boolean isModifiable(byte residue) {
        return modifiables[residue];
    }
    
    public double getDbinwidthNeutralLoss() {
        return dbinwidthNeutralLoss;
    }
    public double getDbinwidthMassShift() {
        return dbinwidthMassShift;
    }
    public double getMassShift() {
        return massShift;
    }
    public char getSymbol() {
        return symbol;
    }
    public String getModInfo() {
        return info;
    }
    public boolean equals(DiffMod o) {
        return massShift == o.massShift;
    }
    public int hashCode() {
        return hashcode; 
    }
}



