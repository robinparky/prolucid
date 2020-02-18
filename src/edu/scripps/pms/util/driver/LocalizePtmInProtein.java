

import java.util.Iterator;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.FileReader;
import java.io.FileInputStream;
import java.io.File;
import java.util.ArrayList;
import java.util.regex.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import edu.scripps.pms.mspid.ProteinDatabase;
import edu.scripps.pms.util.seq.Fasta;

/*
*LocalizePtmInProtein gets fastas, and peptides from proteins.
*It collects data on how frequent breaks or modifications occur in a position of a fasta sequence.
*/

// the first position is 1, not 0
public class LocalizePtmInProtein
{

	private ArrayList<String> filenames= new ArrayList<String>();
        private static String dbname = "";
        //private static String USAGE = "java LocalizePtmInProtein fasta_file_name phospho_result_file_names"; 
        private static String USAGE = "USAGE: localizeptminprotein fasta_file_name phospho_result_file_names"; 
   ///lustre/people/applications/yates/dbase/SGD_S-cerevisiae_na_06-30-2008_reversed.fasta     
   ///data/2/rpark/ip2_data/catclw/Jamie/2010_01_alpha_FT_c_2011_01_07_21_2922/search/upload2011_01_10_12_8770/phospho/phosphoresult.txt
           
	public static void main(String[] args) 
	{
		LocalizePtmInProtein pc= new LocalizePtmInProtein();
                if(args.length < 2) {
                    System.out.println(USAGE);
                    System.exit(0);
                }
        	try {
                        dbname = args[0];
			ArrayList<String> fileList= pc.setFileList(args);
			pc.analyzeProteins(fileList);
                }catch(IOException e){
                    System.out.println(USAGE);
                    e.printStackTrace();
		}
	}

	
	public ArrayList<String> setFileList(String [] args) throws IOException
	{

		for(int i = 1; i < args.length; i++)
		{
               		filenames.add(args[i]);
		}	
		
		return filenames;
	}	
        // not used
        private int [] getSeqAndAccessionIndex(String header) {
            int [] index = new int[5];
            String [] arr = header.split("\t");
            for(int i = 0; i < arr.length; i++) {
                String element = arr[i].trim();
                if(element.equals("Proteins")) {
                    index[0] = i;
                }else if(element.equals("Mod Sequence")) {
                    index[1] = i;
                }else if(element.equals("Sequence")) {
                    index[2] = i;
                } else if(element.equals("Localization Score")) {
                    index[3] = i;
                } else if(element.equals("Protein Descriptions")) {
                    index[4] = i;
                }
            } 

            return index;
        }


	public void analyzeProteins(ArrayList<String> filenames) throws IOException
	{
                ProteinDatabase pdb = new ProteinDatabase(dbname);

System.out.println("Corrected Peptide Sequence\tLocatlization Score\tProtein Accession\tModified Residue\tModified Residue Position\tProtein Description");
		for(int i = 0; i <filenames.size(); i++)
		{	
                        String filename = filenames.get(i);
                      
                        if(new File(filename).isDirectory()) {
                            filename = filename + "/DTASelect-filter.txt.phospho";
                        }

                        if(!new File(filename).isFile()) {
                            System.err.println(filename + " does not exist");
                            continue;
                        }
                        FileInputStream fis = new FileInputStream(filename);
                        BufferedReader br = new BufferedReader(new InputStreamReader(fis));                 

                        String line = br.readLine();
                        int [] posindex = getSeqAndAccessionIndex(line);
                        int accindex = posindex[0];
                        int seqindex = posindex[1]; // index of corrected sequence, - indicates there was no corrected sequence
                        int orgseqindex = posindex[2];
		 
                        int ascoreindex = posindex[3];
                        int descriptionindex = posindex[4];
	
                        line = br.readLine();
			while(line != null)
	       		{
                            String [] arr = line.split("\t");

                            line = br.readLine();	

                            if(arr.length < 6) { continue; }
                            String pSequence = arr[seqindex];
                            if("-".equals(pSequence)) { // no corrected seq, use original seq
                                pSequence = arr[orgseqindex];
                                int firstdot = pSequence.indexOf(".");
                                int lastdot = pSequence.lastIndexOf(".");
                                if(firstdot != -1 && lastdot != -1) { // get the sequence between the dots 
                                    String newpSequence = pSequence.substring(firstdot+1, lastdot);
//System.out.println(pSequence + "\t" + newpSequence);
                                    pSequence = newpSequence;
                                }
                            }
                            ArrayList<Integer> indexes= new ArrayList<Integer>();
                            StringBuffer sbSequence= new StringBuffer(pSequence);
                            
                            if(!pSequence.equals(""))
                            {
                                //Check for modifications

                                char c;
                                for(int index=0; index<sbSequence.length(); index++)
                                {
                                        c=sbSequence.charAt(index);
                                        if(c=='*'||c=='@'||c=='#')
                                        {
                                                indexes.add(index);
                                                sbSequence.deleteCharAt(index);
                                                index--;
                                        }
                                        else if(c=='(')
                                        {
                                                indexes.add(index);
                                                sbSequence.delete(index, sbSequence.indexOf(")")+1);
                                                index--;
                                        }
                                }

                            }
                            String seqnomod = sbSequence.toString(); 
                            String proteinaccs = arr[accindex];
                          
//System.out.println("before: " + proteinaccs + "\t");
                            proteinaccs = proteinaccs.substring(1, proteinaccs.length()-2);

//System.out.println("after: " + proteinaccs);
                            //String [] accs = arr[accindex].split(" , ");
                            String [] accs = proteinaccs.split(" , ");
                            String protdescs = arr[descriptionindex];
                            protdescs = protdescs.substring(1, protdescs.length()-2); 
                            String [] descs = protdescs.split(" , ");

                            String ascores = arr[ascoreindex];
                            ascores = ascores.substring(1, ascores.length()-1);                            
                            String [] locscores = ascores.split(", ");

                            for(int k = 0; k < accs.length; k++) {
                                String acc = accs[k];
                                Fasta f =  pdb.accession2Fasta(acc);
                                //String desc = descs[k];
                                String desc = f.getDescription();
                                //String protseq = pdb.accession2Fasta(acc).getSequence();  
                                String protseq = f.getSequence();  
                                int pos = protseq.indexOf(seqnomod);
                                int mlocindex = 0;
                                
                                for(Iterator<Integer> it = indexes.iterator(); it.hasNext();) {
                                    int peppos = it.next();
                                    int protpos = pos + peppos;
//System.out.println("mlocindex: " + mlocindex + "\t" + ascores);
                                     
                                    String mscore = ""; 
                                    if(locscores.length > 1) {
                                        mscore = locscores[mlocindex++];
                                    } else {
                                        mscore = locscores[0];
                                    } 
// the first position is 1, not 0
//System.out.println(acc + "\t" + pSequence + "\t" + peppos + "\t" + protpos + "\t" + pos + "\t" + seqnomod.charAt(peppos-1));
//System.out.println(pSequence + "\t" + locscores[mlocindex++] + "\t" + acc +  "\t" + seqnomod.charAt(peppos-1) + "\t"  + protpos + "\t" + desc);
System.out.println(pSequence + "\t" + mscore + "\t" + acc +  "\t" + seqnomod.charAt(peppos-1) + "\t"  + protpos + "\t" + desc);
//                                    System.out.println(acc + "\t" + seqnomod + "\t" +  "\t" + acc);
                                }

                            }
			}
			br.close();
		}
	}
}














