import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import java.io.*;
import edu.scripps.pms.util.seq.Fasta;
import edu.scripps.pms.util.io.FastaReader;
import edu.scripps.pms.util.io.DTASelectFilterReader;
import edu.scripps.pms.util.dtaselect.Peptide;
import edu.scripps.pms.util.dtaselect.Protein;
import edu.scripps.pms.util.spectrum.PeakList;
import edu.scripps.pms.util.io.SpectrumReader;


/**
 * <p>Program to six-frame translate a nucleotide sequence</p>
 */

// It now ignores frames with more than 2 stop condons
public class CopaSpectraRetriever {

    // mouse
    //public static final String database = "/lustre/people/applications/yates/dbase/UniProt_mouse_05-17-2011_reversed.fasta";
    //public static final String dtaselectfiles = "/data/8/taoxu_on_data8/projects/nhlbi/v1/MouseMitochondrialV1.txt";
    //public static String outputfile = "/data/8/taoxu_on_data8/projects/nhlbi/v1/MouseMitochondrialV1.copa";

    
    //public static final String database = "/lustre/people/applications/yates/dbase/UniProt_mouse_05-17-2011_reversed.fasta";
    //public static final String dtaselectfiles = "/data/8/taoxu_on_data8/projects/nhlbi/v1/MouseProteasome.txt";
    //public static String outputfile = "/data/8/taoxu_on_data8/projects/nhlbi/v1/MouseProteasome.copa";

    //public static final String database = "/lustre/people/applications/yates/dbase/UniProt_mouse_05-17-2011_reversed.fasta";
    //public static final String dtaselectfiles = "/data/8/taoxu_on_data8/projects/nhlbi/v1/Mouse_Nucleus.txt";
    //public static String outputfile = "/data/8/taoxu_on_data8/projects/nhlbi/v1/Mouse_Nucleus.copa";

     // Human 
//    public static final String database = "/lustre/people/applications/yates/dbase/UniProt_Human_05-17-2011_reversed.fasta";
  //  public static final String dtaselectfiles = "/data/8/taoxu_on_data8/projects/nhlbi/v1/HumanProteasome.txt";
  //  public static String outputfile = "/data/8/taoxu_on_data8/projects/nhlbi/v1/HumanProteasome.copa";

//    public static final String database = "/lustre/people/applications/yates/dbase/UniProt_Human_05-17-2011_reversed.fasta";
//    public static final String dtaselectfiles = "/data/8/taoxu_on_data8/projects/nhlbi/v1/HumanMitochondrialV1.txt";
//    public static String outputfile = "/data/8/taoxu_on_data8/projects/nhlbi/v1/HumanMitochondrialV1.copa";


    // C_elegan
    public static String database = "/lustre/people/applications/yates/dbase/wormbase_c-elegans_225_NoGeneDuplicates_05-10-2011_reversed.fasta";
    public static String dtaselectfiles = "/data/8/taoxu_on_data8/projects/nhlbi/v1/C_elegan_Mitochondria.txt";
    public static String outputfile = "/data/8/taoxu_on_data8/projects/nhlbi/v1/C_elegan_Mitochondria.copa.test";

    public static HashMap<String, Fasta> ac2Fasta = new HashMap();

    /**
     * Call this to get usage info, program terminates after call.
     */
    public static void help() {
	System.out.println("usage: java CopaSpectraRetriever file_contains_all_the_folders output_file_name fasta_file");
	System.exit(-1);
    }

    public static int getPeptideLocation(String peptide, Fasta protein) {
        return protein.getSequence().indexOf(peptide);
    }
    public static HashMap<String, Fasta> getProteinMap(String database) throws IOException {
        HashMap<String, Fasta> ac2fasta = new HashMap(1000000);
        FileInputStream fis = new FileInputStream(new File(database));
        Iterator<Fasta> fastas = FastaReader.getFastas(fis);
        while(fastas.hasNext()) {
            Fasta f = fastas.next();
            String myac = f.getAccession();
   //System.out.println("protien ac " + myac + "\tdescription " + f.getDescription());
            //if(myac != null && !myac.startsWith("Reverse")) {
            if(myac != null) {
                ac2fasta.put(myac, f);
            }
        }
        return ac2fasta;
    }
    public static void main(String[] args) throws Exception{

	BufferedReader br = null;
        if(args.length < 3) {
            help();
        } else {
            dtaselectfiles = args[0];
            outputfile = args[1];
            database = args[2];
        }
        //PrintStream ps = new PrintStream("/data/8/taoxu_on_data8/projects/nhlbi/test1.copa");
        PrintStream ps = new PrintStream(outputfile);
        HashMap<String, Peptide> seqcharge2BestPeptide = new HashMap(10000000);
        HashMap<String, String> seqcharge2BestCopaInfo = new HashMap(10000000);
        ac2Fasta =  getProteinMap(database);
        int copacount = 0;
	try {
	    br = new BufferedReader(new FileReader(dtaselectfiles));
            String line = null;
            line = br.readLine();
            int numpeptide = 0;
            int numspectraadded = 0;
            HashSet<String> uniquepepitdes = new HashSet(10000000);
            while(line != null) {
                line = line.trim();
                String folder = line;
                String dtafilter = folder + "/DTASelect-filter.txt";
                File tempdtafilterfile = new File(dtafilter);
                System.out.println("Now processing " + line);
                //if(!"".equals(line)) {
                if(tempdtafilterfile.exists()) {
                    DTASelectFilterReader reader = new DTASelectFilterReader(dtafilter);
                    HashMap<String, Peptide[]> filename2IdentifiedPeptides = new HashMap(1000);
                    ArrayList<String> proteinacs = new ArrayList(); // for protein group
                    for (Iterator<Protein> itr = reader.getProteins(); itr.hasNext();) {

                        Protein p = itr.next();
                         
                        String proteinid = p.getAccession();
//System.out.println("protien id: " + proteinid);
                        proteinacs.add(proteinid);
                        
                        for(Iterator<Peptide> pepItr=p.getPeptides(); pepItr.hasNext(); ) {
                            numpeptide++;
                            Peptide peptide = pepItr.next();
                            String pepseq = peptide.getSeqWithNoModification();
                            String seq = peptide.getSequence();
                            String filename = peptide.getFileName();
                            String scan = peptide.getLoScan();
                            String charge = peptide.getChargeState();
                            int scanInt = Integer.parseInt(scan);
                            
                            double spvalue = peptide.getSpScoreValue();
                            String seqcharge = seq + charge; 
                            uniquepepitdes.add(seqcharge);

                            Peptide bestpeptide = seqcharge2BestPeptide.get(seqcharge);
                            if(bestpeptide == null || spvalue > bestpeptide.getSpScoreValue()) {
                                numspectraadded++;
                                peptide.setProteinAccessions(proteinacs);
                                seqcharge2BestPeptide.put(seqcharge, peptide);

                                String filescan = filename + scan;
                            
                                Peptide [] peptides = filename2IdentifiedPeptides.get(filename);
                                if (peptides == null) {
                                    peptides = new Peptide[1000000];
                                    filename2IdentifiedPeptides.put(filename, peptides);
                                }
                                if(peptides[scanInt] == null || peptides[scanInt].getSpScoreValue() < peptide.getSpScoreValue()) {
                                    peptides[scanInt] = peptide;
                                } 

                            }
                            
                            //System.out.println("File name: " + filename + "\t" + scan);

                        }
                        if(p.getNumPeptides() != 0) {
                            proteinacs = new ArrayList();
                        }
                    }
                    System.out.println("finished reading DTASelect-filter.txt");
                    reader.close();
                    String [] arrs = folder.split("/");
                    //System.out.println("experiment_id: " + arrs[arrs.length - 1]);
                    String experimentId = arrs[arrs.length - 1];
                    if(!experimentId.startsWith("projects")) {
                        experimentId = arrs[arrs.length - 2]; 
                    }

                    for(Iterator<String> fileit = filename2IdentifiedPeptides.keySet().iterator(); fileit.hasNext();) {
                        String filename = fileit.next();
                        Peptide [] peptides = filename2IdentifiedPeptides.get(filename);
                        String ms2file = folder + "/" + filename + ".ms2";
System.out.println("Now reading " + filename);
                        File tempms2file = new File(ms2file);
                        if(tempms2file.exists()) {
                            SpectrumReader sr = new SpectrumReader(ms2file, "ms2");
                            Iterator<PeakList> plit = sr.getSpectraList().iterator();
                            //Iterator<PeakList> plit = sr.getSpectra();
    System.out.println("finished reading " + filename + "\tnow getting COPA info");
                            while(plit.hasNext()) {
                                PeakList pl = plit.next();
                                int scan = pl.getLoscan();
                                Peptide peptide = peptides[scan];
                                if(peptide != null) {
                                    
                                    String seqcharge = peptide.getSequence() + peptide.getChargeState(); 
                                    Peptide bestpeptide = seqcharge2BestPeptide.get(seqcharge);
                                    if(bestpeptide == peptide) {
                                        String copainfo = getCopaInfo(pl, peptide, ++copacount, experimentId);
                                        seqcharge2BestCopaInfo.put(seqcharge, copainfo);
                                    }
                                    //System.out.println("CopaInfo: " + copainfo);
                                    //ps.print(copainfo);
                                    //System.out.println("identified peptide: " + peptides[scan].getSequence() + " in " + filename + " scan " + scan );
                                }
                                
                            }
System.out.println("Number of COPA info added: " + seqcharge2BestCopaInfo.size() + " and copacount is: " + copacount);
                        } else {
                            System.out.println("MS2 file " + ms2file + " does not exist");
                        } 
                    }
                     
                } else {

                    System.out.println(dtafilter + " does not exist. Ignore");
                } 
                line = br.readLine();
            } 
            br.close();

            Iterator<String> seqcharges = seqcharge2BestCopaInfo.keySet().iterator();
            while(seqcharges.hasNext()) {
                String seqcharge = seqcharges.next();
                ps.print(seqcharge2BestCopaInfo.get(seqcharge));
            }
            ps.close();
            //System.out.println("Number of Peptide Identified: " + numpeptide + "\tNumber of unique peptides: " + uniquepepitdes.size() + "\tNumber from Map: " + seqcharge2BestPeptide.size());
            System.out.println("Number of Peptide Identified: " + numpeptide + "\tNumber of Spectra added: " + numspectraadded + "\tNumber of unique peptides: " + uniquepepitdes.size());
            

	} finally {
	    //tidy up
	    if(br != null){
		br.close();
	    }
	}
    }
    private static String getProteinInfo(ArrayList<String> acs) {
        Iterator<String> it = acs.iterator();
        StringBuffer sbac = new StringBuffer();
        StringBuffer sbdesc = new StringBuffer();
        StringBuffer sblen = new StringBuffer();
        sbac.append("UNIPROTIDS:::");
        sbdesc.append("DESCS:::");
        sblen.append("LEN:::");
        while(it.hasNext()) {
            String ac = it.next();
            Fasta f = ac2Fasta.get(ac);
//System.out.println("ac: " + ac + "\tand ac2Fasta size: " + ac2Fasta.size() + "\tfasta: " + f);
            ac = getUniProtAccession(ac);
            sbac.append(ac + ";");
            if(f != null) { 
                sbdesc.append(f.getDescription() + ";");
                sblen.append(f.getSequence().length() + ";");
            } else {
                System.out.println("Fasta is not found for " + ac);
            }
                        
        }
        return sbac + "|||" + sbdesc + "|||" + sblen + "|||";
    }
    private static String getReverseInfo(ArrayList<String> acs) {
        Iterator<String> it = acs.iterator();
        while(it.hasNext()) {
            String ac = it.next();
            if(!ac.startsWith("Reverse")) {
               return "NotReverseHit"; 
            }
        } 
        return "ReverseHit";
    }

    private static String getUniProtAccession(String longac) {
        String [] arr = longac.split("\\|");
        String ac = arr[0];
        if(arr.length > 1) {
            if(longac.startsWith("Reverse")) {
                ac = arr[0] + "|" + arr[1];
            } else {
                ac = arr[1];
            }
        }
        return ac;
    }
    private static String getCopaInfo(PeakList pl, Peptide peptide, int count, String experimentId) {
        StringBuffer sb = new StringBuffer(5000);
        sb.append("H|||");
        sb.append("PEPID:::" + count + "|||"); // dummy copa peptide id
        sb.append("SPECTRUMID:::" + count + "|||"); // dummy copa spectra id
        sb.append("MZ:::" + pl.getPrecursorMz() + "|||"); // 4 precursor mz
        sb.append("SEQ:::" + peptide.getSequence()); sb.append("|||"); // 5 peptide seq with modifications
        sb.append("MOD:::" + getDiffModInfo(peptide.getSequence())); sb.append("|||"); // 6 DiffModInfo
        sb.append("CHARGE:::" + peptide.getChargeState()); sb.append("|||"); // 7 precursor charge state 
         
        //sb.append("protein accessions"); sb.append("|"); // 8 protein id
        //sb.append("protein descriptions"); sb.append("|"); // 9 protein descriptions
        //sb.append("protein lengths"); sb.append("|"); // 10 protein lengths
if(peptide.getProteinAccessions() == null) {System.out.println("Protein info for " + peptide.getSequence() + " is " + peptide.getProteinAccessions());}
        sb.append(getProteinInfo(peptide.getProteinAccessions()));

        sb.append("EXPID:::" + experimentId); sb.append("|||"); // 11 protein lengths
       
          
        sb.append("SPECTRUMFILE:::" + peptide.getFileNameScanCharge()); sb.append("|||"); // 12 file name field that contains file, scan and charge
        sb.append("DTACONFIDENCE:::" + peptide.getConf()); sb.append("|||"); // 13 peptide id confidence 
        sb.append("XCORR:::" + peptide.getXCorr()); sb.append("|||"); // 14 peptide id XCorr 
        sb.append("DELTACN:::" + peptide.getDeltCN()); sb.append("|||"); // 15 peptide id DeltaCN 
        sb.append("ZSCORE:::" + peptide.getSpScore()); sb.append("|||"); // 16 peptide id ZScore 
        sb.append("PI:::" + peptide.getPi()); sb.append("|||"); // 17 peptide id pI
        sb.append("REVERSE:::" + getReverseInfo(peptide.getProteinAccessions())); // ReverseHit or NotReverseHit

        sb.append("|||\n"); // 17 peptide id pI
       
//System.out.println("COPA H line: " + sb);
        sb.append(pl.getFragmentIons()); 
        
        return sb.toString();
    }
    private static String getDiffModInfo(String seq) {
        return "DiffMod";
    }
}

