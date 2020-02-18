import java.io.BufferedReader;
import java.io.FileReader;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.RNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.SequenceTools;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.symbol.SymbolList;

/**
 * <p>Program to six-frame translate a nucleotide sequence</p>
 */

// It now ignores frames with more than 2 stop condons
public class ThreeFragmeTranslation {
    /**
     * Call this to get usage info, program terminates after call.
     */
    public static void help() {
	System.out.println(
		"usage: java ThreeFragmeTranslation <file> <format eg fasta> <dna|rna>");
	System.exit( -1);
    }

    public static String getDescription(Sequence seq) {
        String def = (String)seq.getAnnotation().asMap().get("description_line");
        String description = "";
        if(def != null) {
            int space = def.indexOf(" ");
//System.out.println("description: " + def);
//System.out.println("description: " + def.substring(space));
//System.out.println("sequence: " + seq.seqString());
        }
        return description; 
    }
    public static void main(String[] args) throws Exception{
	if (args.length != 3) {
	    help();
	}

	BufferedReader br = null;

	//file format (eg fasta)
	String format = args[1];

	//sequence type (eg dna)
	String alpha = args[2];

	try {
	    br = new BufferedReader(new FileReader(args[0]));

	    SequenceIterator seqi =
		(SequenceIterator)SeqIOTools.fileToBiojava(format, alpha, br);

	    //for each sequence
	    while(seqi.hasNext()){
		Sequence seq = seqi.nextSequence();
//getDescription(seq);
//System.out.println("name: " + seq.getName() + "\tAnnotation: " + seq.getAnnotation().asMap().get("description_line"));
    //if(true) 
    //continue;
                if(seq.length() < 10) continue;
		//for each frame
		for (int i = 0; i < 3; i++) {
		    SymbolList prot;
		    Sequence trans;

		    //take the reading frame
		    SymbolList syms = seq.subList(
			    i+1,
			    seq.length() - (seq.length() - i)%3);


		    //if it is DNA transcribe it to RNA
		    if(syms.getAlphabet() == DNATools.getDNA()){
			//before BJ1.4 use this method
			syms = RNATools.transcribe(syms);
			//after BJ1.4 use this method
			//syms = DNATools.toRNA(syms);
		    }

		    //output forward translation to STDOUT
		    prot = RNATools.translate(syms);
                    Annotation anno = new SimpleAnnotation();
                    String plusacc = seq.getName()+ "_Frame_Plus_" + i;
                    String minusacc = seq.getName()+ "_Frame_Minus_" + i;
                    String olddefline = "" + seq.getAnnotation().asMap().get("description_line");
                    String plusdefline = plusacc + " " + olddefline + " translation frame +" + i;
                    String minusdefline = minusacc + " " + olddefline + " translation frame -" + i;

                    anno.setProperty("description_line", plusdefline);
		    trans = SequenceTools.createSequence(prot, "", plusacc, anno);
                    String seqstring = trans.seqString();
                    int numstop = 0;
                    for(int j = 0; j < seqstring.length(); j++) {
                        if(seqstring.charAt(j) == '*') numstop++;
                    }
                    if(numstop > 1) {
                        //System.out.println("contains *, " + plusacc+ ": \t" + seqstring);
                    } else {

                        //System.out.println("no *, " + plusacc + ": \t" + seqstring);
		        SeqIOTools.writeFasta(System.out, trans);
                    }
			    //Annotation.EMPTY_ANNOTATION);
		    //SeqIOTools.writeFasta(System.out, trans);



		    //output reverse frame translation to STDOUT
		    syms = RNATools.reverseComplement(syms);
		    prot = RNATools.translate(syms);
                    Annotation anno1 = new SimpleAnnotation();
                    anno1.setProperty("description_line", minusdefline);
		    trans = SequenceTools.createSequence(prot, "", minusacc, anno1);
			 //   Annotation.EMPTY_ANNOTATION);
		    //SeqIOTools.writeFasta(System.out, trans);
		}
	    }
	}
	finally {
	    //tidy up
	    if(br != null){
		br.close();
	    }
	}
    }
}

