package edu.scripps.pms.util.driver;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import java.io.*;

import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.mspid.MassCalculator;
import edu.scripps.pms.mspid.MassSpecConstants;
import edu.scripps.pms.util.PmsUtil;
import edu.scripps.pms.util.enzyme.Protease;

/**
 * <p>Program to six-frame translate a nucleotide sequence</p>
 */

// It now ignores frames with more than 2 stop condons
public class ProteinTrypticPeptideNumber {



    /**
     * Call this to get usage info, program terminates after call.
     */
    public static void help() {
	System.out.println(
		"usage: proteintrypticpeptides <fastafile>");
	System.exit( -1);
    }



    /* Loop through each residue in the database sequence for this
     * locus and accumulate molecular weight.  Use a simple algorithm
     * to estimate pI.
     */
    public static float [] CalculatepI(String Sequence) {
	//String         Sequence = this.TrimmedSequence();
	int            Length = Sequence.length();
	int            Looper;
	int            CountLys = 0;
	int            CountArg = 0;
	int            CountHis = 0;
	int            CountAsp = 0;
	int            CountGlu = 0;
	int            CountCys = 0;
	int            CountTyr = 0;
	float          CurrentPH = 7.0f;
	float          CurrentJump = 3.5f;
	float          CurrentCharge;
	float          LastCharge = 0;
	float          MWAccum = 0f;
	char           CurrentResidue;
	if (Length > 0) {
	    for (Looper = 0; Looper < Length; Looper++) {
		CurrentResidue = Sequence.charAt(Looper);
		switch (CurrentResidue) {
		case 'A':
		    MWAccum += 71.0;
		    break;
		case 'C':
		    MWAccum += 103.0;
		    CountCys ++;
		    break;
		case 'D':
		    MWAccum += 115.0;
		    CountAsp ++;
		    break;
		case 'E':
		    MWAccum += 129.0;
		    CountGlu ++;
		    break;
		case 'F':
		    MWAccum += 147.0;
		    break;
		case 'G':
		    MWAccum += 57.0;
		    break;
		case 'H':
		    MWAccum += 137.0;
		    CountHis ++;
		    break;
		case 'I':
		    MWAccum += 113.0;
		    break;
		case 'K':
		    MWAccum += 128.0;
		    CountLys ++;
		    break;
		case 'L':
		    MWAccum += 113.0;
		    break;
		case 'M':
		    MWAccum += 131.0;
		    break;
		case 'N':
		    MWAccum += 114.0;
		    break;
		case 'P':
		    MWAccum += 97.0;
		    break;
		case 'Q':
		    MWAccum += 128.0;
		    break;
		case 'R':
		    MWAccum += 156.0;
		    CountArg ++;
		    break;
		case 'S':
		    MWAccum += 87.0;
		    break;
		case 'T':
		    MWAccum += 101.0;
		    break;
		case 'V':
		    MWAccum += 99.0;
		    break;
		case 'W':
		    MWAccum += 186.0;
		    break;
		case 'Y':
		    MWAccum += 176.0;
		    CountTyr ++;
		    break;
		}
	    }
	    /* Use a bracketing strategy to identify the isoelectric
	     * point.  Calculate charge at pH of 7, and then move up
	     * 3.5 or down 3.5 as necessary.  Make each successive
	     * move up or down only half as large.  Keep going until
	     * two successive charges reported match to one place past
	     * the decimal point.
	     */
	    CurrentCharge = PmsUtil.ChargeAtPH(CurrentPH, CountLys,
					       CountArg, CountHis,
					       CountAsp, CountGlu,
					       CountCys, CountTyr);
	    while (PmsUtil.RoundTo(CurrentCharge,1) != PmsUtil.RoundTo(LastCharge,1)) {
		//		System.out.println("pH:\t" + new Float(CurrentPH).toString() 
		//              + "\tCharge\t" + new Float(CurrentCharge).toString());
		if (CurrentCharge > 0)
		    CurrentPH += CurrentJump;
		else
		    CurrentPH -= CurrentJump;
		CurrentJump /= 2;
		LastCharge = CurrentCharge;
		CurrentCharge = PmsUtil.ChargeAtPH(CurrentPH,
						   CountLys, CountArg,
						   CountHis, CountAsp,
						   CountGlu, CountCys,
						   CountTyr);
		if ( (CurrentPH > 14) || (CurrentPH < 0) ) {
		    System.out.println("pI can't be figured for " + Sequence);
		    System.exit(0);
		}
	    }
	}
        float [] result = new float[2];
        result[0] = CurrentPH;
        result[1] = MWAccum;
        return result;
	//return CurrentPH;

    }
    public static float CalculateKyteDoolittle(String Sequence) {
        int            Length = Sequence.length();
        int            Looper;
        float          MWAccum = 0f;
        char           CurrentResidue;
        if (Length > 0) {
            for (Looper = 0; Looper < Length; Looper++) {
                CurrentResidue = Sequence.charAt(Looper);
                switch (CurrentResidue) {
                case 'A':
                    MWAccum += 18.0;
                    break;
                case 'C':
                    MWAccum += 25.0;
                    break;
                case 'D':
                    MWAccum -= 35.0;
                    break;
                case 'E':
                    MWAccum -= 35.0;
                    break;
                case 'F':
                    MWAccum += 28.0;
                    break;
                case 'G':
                    MWAccum -= 4.0;
                    break;
                case 'H':
                    MWAccum -= 32.0;
                    break;
                case 'I':
                    MWAccum += 45.0;
                    break;
                case 'K':
                    MWAccum -= 39.0;
                    break;
                case 'L':
                    MWAccum += 38.0;
                    break;
                case 'M':
                    MWAccum += 19.0;
                    break;
                case 'N':
                    MWAccum -= 35.0;
                    break;
                case 'P':
                    MWAccum -= 16.0;
                    break;
                case 'Q':
                    MWAccum -= 35.0;
                    break;
                case 'R':
                    MWAccum -= 45.0;
                    break;
                case 'S':
                    MWAccum -= 8.0;
                    break;
                case 'T':
                    MWAccum -= 7.0;
                    break;
                case 'V':
                    MWAccum += 42.0;
                    break;
                case 'W':
                    MWAccum -= 9.0;
                    break;
                case 'Y':
                    MWAccum -= 13.0;
                    break;
                }
            }
        }
        return MWAccum/10.0f;
    }

    public static int getNumTrypticPeptides(Fasta f, MassCalculator mc, double highLimit, double lowLimit, int minlength, int maxmiscleav) {
        double [] precMasses = mc.getPrecMasses();
        int counts = 0;
        int proteinlength =  f.getLength();
        int lastIndex = f.getLength() - 1;
        byte [] seq = f.getSequenceAsBytes();
            // mass of H2O need to be added
            // because the mass of the residues does not count H2O
        double mass = MassSpecConstants.MASSH3O;
        double tempmass = 0;
        int templeft = 0;
        //int rightEnd = -1;
        int leftEnd = 0;


//System.out.println("Protein Length: " + (lastIndex + 1) + "\thighLimit: " + highLimit + "\tlowLimit: " + lowLimit);




        while(leftEnd < proteinlength) {

            int rightEnd = leftEnd;
            
            while(rightEnd < lastIndex) {
                rightEnd++;
                if(rightEnd == lastIndex || seq[rightEnd-1] == 'K' ||  seq[rightEnd-1] == 'R') { // tryptic end
                    if(rightEnd - leftEnd > minlength-1) {
                        double precmass = mc.getPrecursorMass(seq, leftEnd, rightEnd);
                        if(precmass >= lowLimit && precmass <= highLimit) {
//System.out.println("Protein Length: " + (lastIndex + 1) + "\tleftEnd: " + leftEnd + "\trightEnd: " + rightEnd +"\tmass: " + precmass + "\thighLimit: " + highLimit + "\tlowLimit: " + lowLimit);
                            int miscleav = 0;
                            for(int i = leftEnd; i < rightEnd; i++) {
                                if(seq[i] == 'K' || seq[i] == 'R') {
                                    miscleav++;
                                }
                            }
                            if(miscleav <= maxmiscleav) {
                                counts++;
                            }
                        }
                    }
                    break;
                }
            }
            leftEnd = rightEnd+1;

        }

        return counts;

    }
    //public static int getNumTrypticPeptides(Fasta f, Protease p, MassCalculator mc, double highLimit, double lowLimit, int minlength) {
    public static int getNumTrypticPeptides(Fasta f, MassCalculator mc, double highLimit, double lowLimit, int minlength) {
        double [] precMasses = mc.getPrecMasses();
        int counts = 0;
        int proteinlength =  f.getLength();
        int lastIndex = f.getLength() - 1;
        byte [] seq = f.getSequenceAsBytes();
            // mass of H2O need to be added
            // because the mass of the residues does not count H2O
        double mass = MassSpecConstants.MASSH3O;
        double tempmass = 0;
        int templeft = 0;
        //int rightEnd = -1;
        int leftEnd = 0;


//System.out.println("Protein Length: " + (lastIndex + 1) + "\thighLimit: " + highLimit + "\tlowLimit: " + lowLimit);




        while(leftEnd < proteinlength) {

            int rightEnd = leftEnd;
            
            while(rightEnd < lastIndex) {
                rightEnd++;
                if(rightEnd == lastIndex || seq[rightEnd-1] == 'K' ||  seq[rightEnd-1] == 'R') { // tryptic end
                    if(rightEnd - leftEnd > minlength-1) {
                        double precmass = mc.getPrecursorMass(seq, leftEnd, rightEnd);
                        if(precmass >= lowLimit && precmass <= highLimit) {
//System.out.println("Protein Length: " + (lastIndex + 1) + "\tleftEnd: " + leftEnd + "\trightEnd: " + rightEnd +"\tmass: " + precmass + "\thighLimit: " + highLimit + "\tlowLimit: " + lowLimit);

                            counts++;
                        }
                    }
                    break;
                }
            }
            leftEnd = rightEnd+1;

        }

        return counts;

    }

    public static void main(String[] args) throws Exception{

        
        HashSet<String> acset = new HashSet();
        ArrayList<String> aclist = new ArrayList();
        HashMap<String, Fasta> ac2protein = new HashMap();

        Protease p = new Protease("Trypsin");
        p.setType(true);

        p.addCleavageSite('R');
        p.addCleavageSite('K');

        MassCalculator mc = new MassCalculator();


	try {

            String fastaFileName = args[0];

           
             
            FileInputStream fis = new FileInputStream(new File(fastaFileName));
            for (Iterator<Fasta> it = FastaReader.getFastas(fastaFileName); it.hasNext();) {
                Fasta f = it.next();
                String acc = f.getAccession();
                ac2protein.put(acc, f);
                aclist.add(acc);
                
            }
            fis.close();
            System.out.println("Locus\tKD\tpI\tMolWt\tNumTrypticPeptides");
            for(Iterator<String> it = aclist.iterator(); it.hasNext();) { 
                String ac = it.next();
                Fasta f = ac2protein.get(ac);
                if(f != null) { 
                    double kd = CalculateKyteDoolittle(f.getSequence());
                    float [] piandmw = CalculatepI(f.getSequence());
                    double pi = piandmw[0];
                    double mw = piandmw[1];
                    //int numtrypticpepitdes = getNumTrypticPeptides(f, p, mc, 4800, 600, 6);
                    int numtrypticpepitdes = getNumTrypticPeptides(f, mc, 4800, 600, 6);
                    System.out.println(ac + "\t" + kd + "\t" + pi + "\t" + mw + "\t" + numtrypticpepitdes);
                    System.out.println(ac + "\t" + kd + "\t" + pi + "\t" + mw + "\t" + numtrypticpepitdes);
                } else {
        
                    System.out.println(ac + "not found");
                }
            }
            

	} catch(Exception e){
            help();
            e.printStackTrace(); 
	} finally {
	    //tidy up
	}
    }
}

