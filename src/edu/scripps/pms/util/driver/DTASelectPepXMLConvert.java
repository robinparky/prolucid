/**
 * @file DTASelectPepXMLConvert.java
 * This is the source file for DTASelectPepXMLConvert
 * @author Tao Xu
 * @date $Date
 */


import java.util.*;
import java.io.*;
import java.text.DecimalFormat;

import edu.scripps.pms.util.io.*;
import edu.scripps.pms.mspid.SearchParams;

import java.io.BufferedReader;
import java.util.*;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.FileReader;
import java.io.RandomAccessFile;

import edu.scripps.pms.util.dtaselect.Protein;
import edu.scripps.pms.util.dtaselect.Peptide;
import edu.scripps.pms.util.dtaselect.DiffModSite;


import org.jdom.*;
import org.jdom.output.*;
import org.jdom.input.*;

public class DTASelectPepXMLConvert
{
    public static final int SCANLENGTH = 6;
    private static final double [] monoAaMasses = new double[256];
    private static final double [] avgAaMasses = new double[256];
    private static final char [] symbols = new char[256];
    public static final DecimalFormat fourDigits = new DecimalFormat("0.0000");

   
    // populate the msmsrunsummary element with the information  
    private static double [] getMsmsRunElementFromSqt(Element msmsrunsummary, Element searchsummary, String sqtfile) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(sqtfile));
        double[] symbol2mass = new double[256];
        String line = br.readLine();
        ArrayList<Element> staticmods = new ArrayList<Element>();
        char diffmodsymbol = 35;
        while(line != null) {
            String temparr[] = line.split("\t");
            if(line.startsWith("S")) {
                break;
            } else if(line.startsWith("H\tStaticMod")) {
                String [] arr = line.split("\t")[2].split("="); 
                Element sm = new Element("aminoacid_modification");
                sm.setAttribute("aminoacid", arr[0]);
                sm.setAttribute("mass", fourDigits.format(Double.parseDouble(arr[1])));
                double massshift = Double.parseDouble(arr[1]) - monoAaMasses[arr[0].charAt(0)];
//System.out.println("massshift: " + massshift+"");
                sm.setAttribute("massdiff", fourDigits.format(massshift));
                sm.setAttribute("variable", "N");
                //sm.etAttribute("massdiff", "0"); need to be handled
                //searchsummary.addContent(sm); 
                staticmods.add(sm); 
            } else if(line.startsWith("H\tDiffMod")) {
                
                String [] arr = line.split("\t")[2].split("="); 
                String aas = arr[0];
                char symbol = aas.charAt(aas.length()-1);
                if(!Character.isLetter(symbol)) {
                    aas = aas.substring(0, aas.length()-1);
                } else {
                    symbol = diffmodsymbol++; // no symbol
                } 
                double massshift = Double.parseDouble(arr[1]);
                for(int i = 0; i < aas.length(); i++) {
                    Element dm = new Element("aminoacid_modification");
                    dm.setAttribute("aminoacid", ""+aas.charAt(i));
                    dm.setAttribute("massdiff", fourDigits.format(massshift));
                    dm.setAttribute("variable", "Y");
                    //sm.etAttribute("massdiff", "0"); need to be handled
                    double mass = monoAaMasses[aas.charAt(i)]; 
                    if(symbol != 0) {
                        dm.setAttribute("mass", fourDigits.format((mass + massshift)) );
                        dm.setAttribute("symbol", ""+symbol);
                        symbol2mass[symbol] = massshift;
                    }
                    searchsummary.addContent(dm); 
                }
            } else if(line.startsWith("H\tSQTGenerator\t")) {
                searchsummary.setAttribute("search_engine", temparr[2]); 
            } else if(line.startsWith("H\tPrecursorMasses")) {
                String masstype = "mono".equals(temparr[2])? "monoisotopic" : "avgisotopic";
                searchsummary.setAttribute("precursor_mass_type", masstype); 
            } else if(line.startsWith("H\tFragmentMasses")) {
                String masstype = "mono".equals(temparr[2])? "monoisotopic" : "avgisotopic";
                searchsummary.setAttribute("fragment_mass_type", masstype); 
            } else if(line.startsWith("H\tEnzymeName")) {
            } 
            line = br.readLine();
        }

        for(Iterator<Element> it = staticmods.iterator(); it.hasNext();) {

                searchsummary.addContent(it.next()); 
        }
        searchsummary.setAttribute("out_data_type", "out");
        searchsummary.setAttribute("out_data", ".tgz");
        searchsummary.setAttribute("search_id", "1");
        br.close();

        return symbol2mass;
    }

    public static void outputPepxml(String dtaselectfile, String sqtfile, 
                    String basename, String outputfile) throws IOException {
        
	Element rootEle = new Element("msms_pipeline_analysis");
	Element msmsrunsummary = new Element("msms_run_summary");
        
        msmsrunsummary.setAttribute("base_name", basename);
        Element searchsummary = new Element("search_summary");
        searchsummary.setAttribute("base_name", basename);
        
        Element searchdatabase = new Element("search_database");
        searchdatabase.setAttribute("local_path", "");
        searchdatabase.setAttribute("type", "AA");
        searchsummary.addContent(searchdatabase);

        msmsrunsummary.addContent(searchsummary);

        double [] symbol2Mass = getMsmsRunElementFromSqt(msmsrunsummary, searchsummary, sqtfile);
	//build jdom

        DTASelectFilterReader reader = new DTASelectFilterReader(dtaselectfile);

        Iterator<Protein> pitr = reader.getProteins();
        HashSet<Peptide> peptides = new HashSet<Peptide>(1000000);

        HashMap<Peptide, HashSet<Protein>> peptide2Proteins = new HashMap<Peptide, HashSet<Protein>>(1000000);

        ArrayList<Protein> aList = new ArrayList<Protein>();
        
        
	//HashSet set = new HashSet();

	int i=0;
        for (Iterator<Protein> itr = pitr; itr.hasNext(); )
        {
            Protein protein = itr.next();

            aList.add(protein);
            if(protein.getPeptideSize() > 0)
            {
		for(Iterator<Peptide> pepItr=protein.getPeptides(); pepItr.hasNext(); )
		{
		    Peptide peptide = pepItr.next();
                    peptides.add(peptide);
                    HashSet<Protein> proteins = peptide2Proteins.get(peptide);
                    if(proteins == null) {
                        proteins = new HashSet<Protein>();
                        peptide2Proteins.put(peptide, proteins);

                    }
                    proteins.addAll(aList);
                    
		    i++;
		}

                //aList.clear();
                aList = new ArrayList<Protein>();
            }
        }

        System.out.println("Number of peptide entries: " + i + "\tnumber of unique pepitdes: " + peptides.size()); 
        int index = 0;
        for(Iterator<Peptide> pepit = peptides.iterator(); pepit.hasNext();) {

            boolean isPhosphoPeptide = false;

            Peptide peptide = pepit.next();
            if("1".equals(peptide.getChargeState().trim())) continue; // Ascore cannot hanlde charge spectra, ignore
            HashSet<Protein> proteins = peptide2Proteins.get(peptide);
            Protein protein = proteins.iterator().next();
            //		if(protein.getLocus().equals("IPI00430839.1"))
            Element spectrumQueryEle = new Element("spectrum_query");
            Element searchResultEle = new Element("search_result");
            //searchResultEle.setAttribute("assumed_charge", peptide.getChargeState());
            Element searchHitEle = new Element("search_hit");
            searchHitEle.setAttribute("hit_rank", "1");
            String seq = peptide.getSequence();
            //searchHitEle.setAttribute("peptide", seq.substring(2, seq.length()-2));
            searchHitEle.setAttribute("peptide", getSeqWithNoModification(peptide));
            searchHitEle.setAttribute("peptide_prev_aa", seq.substring(0,1));
            searchHitEle.setAttribute("peptide_next_aa", seq.substring(seq.length()-1, seq.length()));
            searchHitEle.setAttribute("protein", protein.getLocus());
            searchHitEle.setAttribute("num_tot_proteins", proteins.size()+"");
            searchHitEle.setAttribute("calc_neutral_pep_mass", peptide.getCalcMHplus()+"");
            double massdiff = Double.parseDouble(peptide.getMhPlus()) - Double.parseDouble(peptide.getCalcMHplus());
            searchHitEle.setAttribute("massdiff", fourDigits.format(massdiff));
            searchHitEle.setAttribute("is_rejected", "0");

            Element score = new Element("search_score");
            score.setAttribute("name", "xcorr");
            score.setAttribute("value", peptide.getXCorr());
            searchHitEle.addContent(score);

            score = new Element("search_score");
            score.setAttribute("name", "deltacn");
            score.setAttribute("value", peptide.getDeltCN());
            searchHitEle.addContent(score);
            searchResultEle.addContent(searchHitEle);
            String scannum = padzeros2ScanNum(peptide.getLoScan());
            String chargestate = peptide.getChargeState();
            String filescancharge = peptide.getFileName() + "." + scannum + "." + scannum + "." + chargestate;
            spectrumQueryEle.setAttribute("spectrum", filescancharge);

            spectrumQueryEle.setAttribute("start_scan", peptide.getScanNum());
            spectrumQueryEle.setAttribute("end_scan", peptide.getScanNum());

            spectrumQueryEle.setAttribute("precursor_neutral_mass", peptide.getMhPlus());
            spectrumQueryEle.setAttribute("assumed_charge", peptide.getChargeState());
            spectrumQueryEle.setAttribute("index", ++i+"");

            ArrayList<DiffModSite> diffmods = getDiffModSites(peptide, symbol2Mass);
            if(diffmods != null && diffmods.size() > 0) {
                for(Iterator<DiffModSite> it = diffmods.iterator(); it.hasNext();) {
                    DiffModSite d = it.next();
                    isPhosphoPeptide = isPhosphoPeptide || isPhosphoSite(d); 
                    Element diffm = new Element("modification_info");
                    Element siteinfo = new Element("mod_aminoacid_mass");
                    siteinfo.setAttribute("position", d.getSite()+"" );
                    siteinfo.setAttribute("mass", fourDigits.format(d.getMassShift()+ monoAaMasses[d.getResidue()])+"" );
                    diffm.addContent(siteinfo);
                    searchHitEle.addContent(diffm);

                }
            }

            spectrumQueryEle.addContent(searchResultEle);		    
            msmsrunsummary.addContent(spectrumQueryEle);


        }


//	System.out.println(i);



        rootEle.addContent(msmsrunsummary);
	Document doc = new Document(rootEle);
        HashMap piMap = new HashMap( 2 );
        piMap.put( "type", "text/xsl" );
        piMap.put( "href", "/tools/bin/TPP/tpp-dshteynb/schema/pepXML_std.xsl" );
        ProcessingInstruction pi = new ProcessingInstruction( "xml-stylesheet", piMap );

        doc.getContent().add( 0, pi );
	//OutputStream os = new FileOutputStream(new File("converted_pepxml.xml")); //(filePath + "census_chro.xml");
	OutputStream os = new FileOutputStream(new File(outputfile)); //(filePath + "census_chro.xml");
	XMLOutputter outputter = new XMLOutputter();
	outputter.setFormat(Format.getPrettyFormat());
	outputter.output(doc, os);
	os.close();
        reader.close();

    }

    private static boolean isPhosphoSite(DiffModSite d) {
        
        double shift = d.getMassShift();
        char r = d.getResidue();
        return (shift > 79 && shift < 81) && (r == 'S' || r == 'T' || r == 'Y'); 

    }
    public static void main(String args[]) throws Exception
    {
	//Convert DTASelect to pepXML
	if(args.length<3)
	{
	    System.out.println("Usage: java DTASelectPepXMLConvert DTASelect-filter.txt sqt_file, base_name pepxml_file_name");
	    System.exit(0);
	}


        String dtaselectfile = args[0];
        String sqtfile = args[1];
        String basename = args[2];
        String outputfile = args[3];
    
        outputPepxml(dtaselectfile, sqtfile, basename, outputfile);

    }
    private static String padzeros2ScanNum(String scan) {
        StringBuffer sb = new StringBuffer(10);
        for(int i = 0; i < SCANLENGTH-scan.length(); i++) {
            sb.append('0');
        }
        sb.append(scan);
        return sb.toString();
    }

    public static boolean isModifiedPeptide(Peptide p) {
        String sequence = p.getSequence();
         return sequence.indexOf("(") != -1 || sequence.indexOf("*") != -1
                                           || sequence.indexOf("@") != -1
                                           || sequence.indexOf("#") != -1;

    }

    public static String getSeqWithNoModification(Peptide p) {

        StringBuffer tempseq= new StringBuffer(p.getMidSeq());
        char c;
        double mod;

        //modifications
        double[] modifications= new double[tempseq.length()];
        for(int k=0; k<tempseq.length(); k++)
        {
            c= tempseq.charAt(k);
            if(c=='(')
            {
                tempseq.delete(k, tempseq.indexOf(")")+1);
                k--;
            }
            else if(!Character.isLetter(c))
            {
                    //mod=configuration[c];
                tempseq.deleteCharAt(k);
                k--;
            }
        }

        return tempseq.toString();

    }
 
    public static ArrayList<DiffModSite> getDiffModSites(Peptide p, double [] symbol2mass) {
        ArrayList<DiffModSite> diffmods = new ArrayList<DiffModSite>();
        if (isModifiedPeptide(p)) {
            StringBuffer tempseq= new StringBuffer(p.getMidSeq());
            diffmods = new ArrayList<DiffModSite>();
            char c;
            double mod;

            //modifications
            double[] modifications= new double[tempseq.length()];
            for(int k=0; k<tempseq.length(); k++)
            {
                c= tempseq.charAt(k);
                if(c=='(')
                {
                        mod=Double.parseDouble(tempseq.substring(k+1, tempseq.indexOf(")")));
                        if(k==0)
                            modifications[k]+=mod;
                        else
                            modifications[k-1]+=mod;
                        tempseq.delete(k, tempseq.indexOf(")")+1);
                        k--;
                }
                else if(!Character.isLetter(c))
                {
                        //mod=configuration[c];
                        mod=symbol2mass[c];
                        if(k==0)
                            modifications[k]+=mod;
                        else
                            modifications[k-1]+=mod;
                        tempseq.deleteCharAt(k);
                        k--;
                }
            }

            //seqWithNoModification = tempseq.toString();
            for(int i = 0; i < modifications.length; i++) {
                if(modifications[i] != 0)
                    diffmods.add(new DiffModSite(i+1, modifications[i], tempseq.charAt(i)));
            }
        }
        return diffmods;
    }


        // initiatet the values of average and mono aa masses
    static {
        avgAaMasses['G'] =  57.05192f;   monoAaMasses['G'] =  57.0214636f;
        avgAaMasses['A'] =  71.07880f;   monoAaMasses['A'] =  71.0371136f;
        avgAaMasses['S'] =  87.07820f;   monoAaMasses['S'] =  87.0320282f;
        avgAaMasses['P'] =  97.11668f;   monoAaMasses['P'] =  97.0527636f;
        avgAaMasses['V'] =  99.13256f;   monoAaMasses['V'] =  99.0684136f;
        avgAaMasses['T'] = 101.10508f;   monoAaMasses['T'] = 101.0476782f;
        avgAaMasses['C'] = 103.13880f;   monoAaMasses['C'] = 103.0091854f;
        avgAaMasses['L'] = 113.15944f;   monoAaMasses['L'] = 113.0840636f;
        avgAaMasses['I'] = 113.15944f;   monoAaMasses['I'] = 113.0840636f;
        avgAaMasses['X'] = 113.15944f;   monoAaMasses['X'] = 113.0840636f;
        avgAaMasses['N'] = 114.10384f;   monoAaMasses['N'] = 114.0429272f;
        avgAaMasses['O'] = 114.14720f;   monoAaMasses['O'] = 114.0793126f;
        avgAaMasses['B'] = 114.59622f;   monoAaMasses['B'] = 114.5349350f;
        avgAaMasses['D'] = 115.08860f;   monoAaMasses['D'] = 115.0269428f;
        avgAaMasses['Q'] = 128.13072f;   monoAaMasses['Q'] = 128.0585772f;
        avgAaMasses['K'] = 128.17408f;   monoAaMasses['K'] = 128.0949626f;
        avgAaMasses['Z'] = 128.62310f;   monoAaMasses['Z'] = 128.5505850f;
        avgAaMasses['E'] = 129.11548f;   monoAaMasses['E'] = 129.0425928f;
        avgAaMasses['M'] = 131.19256f;   monoAaMasses['M'] = 131.0404854f;
        avgAaMasses['H'] = 137.14108f;   monoAaMasses['H'] = 137.0589116f;
        avgAaMasses['F'] = 147.17656f;   monoAaMasses['F'] = 147.0684136f;
        avgAaMasses['R'] = 156.18748f;   monoAaMasses['R'] = 156.1011106f;
        avgAaMasses['Y'] = 163.17596f;   monoAaMasses['Y'] = 163.0633282f;
        avgAaMasses['W'] = 186.21320f;   monoAaMasses['W'] = 186.0793126f;
    }

}

