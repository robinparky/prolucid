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

public class OrfCreater {
    /**
     * Call this to get usage info, program terminates after call.
     */
    private static int minimumLength = 30;
    public static void help() {
	System.out.println(
		"usage: java OrfCreater <file> <format eg fasta> <dna|rna> <minimum_length\nDefault minimum length is 30");
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
	if (args.length < 3) {
	    help();
	}

	BufferedReader br = null;

	//file format (eg fasta)
	String format = args[1];

	//sequence type (eg dna)
	String alpha = args[2];

        int minlength = args.length < 4? 30 : Integer.parseInt(args[3]);
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
		    SymbolList syms = seq.subList(i+1, seq.length() - (seq.length() - i)%3);


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
                    String plusacc = seq.getName()+ "_FP_" + i; // F for Frame Plus
                    String minusacc = seq.getName()+ "_FM_" + i; // FM for Frame Minus

                    /*
                    String olddefline = "" + seq.getAnnotation().asMap().get("description_line");
                    String plusdefline = plusacc + " " + olddefline + " translation frame +" + i;
                    String minusdefline = minusacc + " " + olddefline + " translation frame -" + i;

                    anno.setProperty("description_line", plusdefline);
		    trans = SequenceTools.createSequence(prot, "", plusacc, anno);
			    //Annotation.EMPTY_ANNOTATION);
		    SeqIOTools.writeFasta(System.out, trans);
                    */
                     String seqstr = prot.seqString();
                     //String [] arr = seqstr.split("*");
                     outputOrfs(plusacc, seqstr, i); 
                     //System.out.println(seqstr);
                    
		    //output reverse frame translation to STDOUT
		    syms = RNATools.reverseComplement(syms);
		    prot = RNATools.translate(syms);
                    seqstr = prot.seqString();
                    outputOrfs(minusacc, seqstr, i);
                    /*
                    Annotation anno1 = new SimpleAnnotation();
                    anno1.setProperty("description_line", minusdefline);
		    trans = SequenceTools.createSequence(prot, "", minusacc, anno1);
			 //   Annotation.EMPTY_ANNOTATION);
		    SeqIOTools.writeFasta(System.out, trans);
                    */
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

    private static void outputOrfs(String accprefix, String seq, int frame) {
        String arr [] = seq.split("\\*");
        long pos = 0 + frame;
         
        for(int i = 0; i < arr.length; i++) {
             
            long end = pos + 3*arr[i].length();
            
            if(arr[i] != null && arr[i].length() >= minimumLength) { 
                String finalseq = arr[i].replaceAll("XX", "");
                System.out.println(">" + accprefix + "_" + pos + "_" + end + "\n" + finalseq);
            }
            pos = end + 3; // 3 is for the stop codon
        }
    }
}

