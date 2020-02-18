
/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2009</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Tao Xu 
 * @version 1.0
 */

//import edu.scripps.pms.pmsdb.*;
import edu.scripps.pms.go.AnnotatedProtein;
import edu.scripps.pms.go.GoAssociationReader;
import edu.scripps.pms.go.GoAssociation;
import edu.scripps.pms.go.GoTerm;
import edu.scripps.pms.go.GeneOntologyExtReader;

import java.util.*;
import java.io.*;

public class AhnaGoAnnotation {
    public static final String USAGE = "\n\n--- java AhnaGoAnnotation minimum_occurence---\n";
    //public static final String resultfile = "/home/taoxu/taoxu_on_data/projects/ali/ahna/compare_peptide_count.txt";
    public static final String resultfile = "/home/taoxu/taoxu_on_data/projects/pms/pmsws/WEB-INF/src/edu/scripps/pms/idcompare/compare_peptide_count.txt";
    //public static final String annotatedresultfile = "/home/taoxu/taoxu_on_data/projects/ali/ahna/annotated_compare_peptide_count.txt";
    public static final String annotatedresultfile = "20101117_annotated_compare_peptide_count.txt";
    public static final String interestedtermfile = "/home/taoxu/taoxu_on_data/projects/ali/ahna/SkopGOterms.txt";
    //public static final String annotationfile = "/data/6/lliao/GO_related/gene_association.goa_mouse";
    public static final String humanAnnotationfile = "/home/taoxu/taoxu_on_data/projects/go/goa201003/gene_association.goa_human";
    public static final String mouseAnnotationfile = "/home/taoxu/taoxu_on_data/projects/go/goa201003/gene_association.goa_mouse";
    public static final String ratAnnotationfile = "/home/taoxu/taoxu_on_data/projects/go/goa201003/gene_association.goa_rat";
    public static final String  geneOntologyExtFile = "/home/taoxu/taoxu_on_data/projects/goa/20100326/gene_ontology_ext.obo";
    public static HashMap<String, AnnotatedProtein> ipi2AnnotatedProtein = new HashMap<String, AnnotatedProtein>(100000);
    public static HashMap<String, AnnotatedProtein> geneSymbol2AnnotatedProtein = new HashMap<String, AnnotatedProtein>(100000);
    public static final HashMap<String, Integer> goid2Counts = new HashMap<String, Integer>(10000);

    public static HashMap<String, GoTerm> goid2Term = null;
    public static HashMap<String, GoTerm> allgoid2Term = new HashMap<String, GoTerm>(100000);
    //public static final HashMap<String, GoTerm> cellularGoid2Term = new HashMap<String, GoTerm>(1000);
    //public static final HashMap<String, GoTerm> functionGoid2Term = new HashMap<String, GoTerm>(1000);
    //public static final HashMap<String, GoTerm> processGoid2Term = new HashMap<String, GoTerm>(1000);
 
    public static void main(String args[]) {
        try {
            //String fileName = args[0];
           /* 
            HashMap<String, GoTerm> goTermMap = getGoTermMap(gotermfile);
            HashMap<String, AnnotatedProtein> ipi2AnnotatedProteins = readAnnotations(annotationfile, goTermMap);
            System.out.println("gotermfile: " + gotermfile);
            System.out.println("associationfile: " + annotationfile);
            System.out.println("IPI\tCellular Component\tMolecular Function\tBiological Process");
            for(Iterator<String> it = ipi2AnnotatedProteins.keySet().iterator(); it.hasNext();) {
                System.out.println(ipi2AnnotatedProteins.get(it.next()).output());
            }
            */
            goid2Term = readInterestedTerms();
            ArrayList<String> annotationfiles = new ArrayList<String>();
            annotationfiles.add(humanAnnotationfile);
            annotationfiles.add(mouseAnnotationfile);
            annotationfiles.add(ratAnnotationfile);
            ipi2AnnotatedProtein = readAnnotations(annotationfiles);
          
            PrintStream ps = new PrintStream(annotatedresultfile); 
            BufferedReader br = new BufferedReader(new FileReader(resultfile));
            String line = br.readLine();
          
            int minoccurence = args.length == 0? 1 : Integer.parseInt(args[0]);

            Set<String> interestedgoids = goid2Term.keySet();
            while(line != null) { 
                if(line.startsWith("P\t")) {
 
                    String [] temparr = line.split("\t");


                    if(temparr.length > 5 && temparr[1].startsWith("contaminant") || temparr[1].startsWith("Rever")) {
                        line = br.readLine();
                        
                        continue;
                    }
                    String genesymbol = temparr[6].split("=")[2].split("\\s+")[0];
//System.out.println("symbol: " + genesymbol);
                    String acc = temparr[1];
                    int tempocc = Integer.parseInt(temparr[7]);


                    AnnotatedProtein ap = geneSymbol2AnnotatedProtein.get(genesymbol);
                    if(ap != null) {
                        
                        ps.println(line + "\t" + ap.output(interestedgoids));
                        if(tempocc >= minoccurence) {
                            count(ap);
                        }
                    } else {

                        ps.println(line + "\t" + genesymbol + "\t\t\t");
                    }
                } else if(line.startsWith("PL")) {
                    ps.println(line + "\t" + "Accession\tCellular_Component\tMolecular_Function\tBiological_Process");
                } else {
                    ps.println(line);
                } 
                line = br.readLine();
            }
            
            br.close();
            ps.close(); 

            StringBuffer cc = new StringBuffer(10000).append("Cellular Component ID\tTerm\tCounts\n");
            StringBuffer bp = new StringBuffer(10000).append("Biological Process ID\tTerm\tCounts\n");
            StringBuffer mf = new StringBuffer(10000).append("Molecular Function ID\tTerm\tCounts\n");
            for(Iterator<String> it = goid2Counts.keySet().iterator(); it.hasNext();) {
                String goid = it.next();
                GoTerm term = goid2Term.get(goid);
                if(term.isMolecularFunction()) addterm(mf, term);
                if(term.isCellularComponent()) addterm(cc, term);
                if(term.isBiologicalProcess()) addterm(bp, term);
            }
            
            System.out.println("Mininum Occurence is: " + minoccurence);  
            System.out.println(mf.toString() + "\n\n\n" + bp.toString() + "\n\n\n" + cc.toString());  
            
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
        }
    }
   
    private static void addterm(StringBuffer sb, GoTerm term) {
        String goid = term.getGoId(); 
        int count = goid2Counts.get(goid);
        sb.append(goid);
        sb.append("\t");
        sb.append(term.getGoName());
        sb.append("\t");
        sb.append(count);
        sb.append("\n");
    }
    private static void count(AnnotatedProtein ap) {
        
        count(ap.getMolucularFunctions()); 
        count(ap.getCellularComponents()); 
        count(ap.getBiologicalProcesses()); 
    }
    private static void count(HashSet<GoTerm> terms) {
        for(Iterator<GoTerm> it = terms.iterator(); it.hasNext();) {
            GoTerm gt = it.next();
            String goid = gt.getGoId();
            if(goid2Counts.keySet().contains(goid)) {
                int counts = goid2Counts.get(goid).intValue() + 1; 
 //System.out.println("counts: " + counts);
                goid2Counts.put(goid, new Integer(counts));
            }
        }
    }
    private static HashMap<String, GoTerm> readInterestedTerms() throws IOException {
        HashMap<String, GoTerm> goid2term = new HashMap<String, GoTerm>(1000);
        BufferedReader br = new BufferedReader(new FileReader(interestedtermfile));
        String line = br.readLine();
        while(line != null) {
            //String arr[] = line.trim().split("\\s+"); // Separated by "whitespace"
            String arr[] = line.trim().split("\t"); // Separated by tab
            if(arr.length == 3) {
                String goid = arr[1];
                if("biological_process".equals(arr[2].trim())) {
                    goid2term.put(arr[1], new GoTerm('b', arr[1], arr[0]));
                } else if("cellular_component".equals(arr[2].trim())) {
                    goid2term.put(arr[1], new GoTerm('c', arr[1], arr[0]));
                } else if("molecular_function".equals(arr[2].trim())) {
                    goid2term.put(arr[1], new GoTerm('m', arr[1], arr[0]));
                } 
                goid2Counts.put(goid, new Integer(0));
            }
             
            line = br.readLine();
        }
        br.close();
        return goid2term;
        //System.out.println("b: " + processGoid2Term.size() + "\tc: " + cellularGoid2Term.size() + "\tm: " + functionGoid2Term.size());
    }

    private static HashMap<String, AnnotatedProtein> readAnnotations(ArrayList<String> annotationfiles) throws IOException  {

        allgoid2Term = readExtGeneOntology(geneOntologyExtFile);            
        HashMap<String, AnnotatedProtein> ipi2AnnotatedProteins = new HashMap<String, AnnotatedProtein>(100000);
        int counter = 0;
        //boolean sortByIntensity = true;
        int unknowngoid = 0;
        HashSet<String> unknow = new HashSet<String>(50000);
        for(Iterator<String> fileit = annotationfiles.iterator(); fileit.hasNext(); ) {
            String file = fileit.next();
            GoAssociationReader gar = new GoAssociationReader(file);
            Iterator<GoAssociation> it = gar.getGoAssociations();
            while (it.hasNext()) {
                GoAssociation g = it.next();
                //swissProtAccs.add(g.getDbObjectId());
                //ipiAccs.add(g.getDbObjectSynonym());
                                
                 
 
                String dos = g.getDbObjectSynonym();
                String [] arr = dos.split("\\|");
    //System.out.println(dos + "\t" + arr.length);
                String goid = g.getGoId().trim();
                //GoTerm gt = gotermmap.get(goid);
                GoTerm gt = allgoid2Term.get(goid);
                if(gt == null) {
                    unknowngoid++; 
                    //System.out.println("Unknow goid: " + goid); 
                    unknow.add(goid); 
            //        continue;
                }

                for(String s : arr) {
                    if(s != null && s.startsWith("IPI")) {
                        
                        AnnotatedProtein ap = ipi2AnnotatedProteins.get(s);
                        if(ap == null) {
                            ap = new AnnotatedProtein(s);
                            ipi2AnnotatedProteins.put(s, ap);
                            
                        }
                       
                        //ap.addAnnotation(gt); // need to add parent terms too
                        addAnnotations(ap, gt, allgoid2Term);
                    }
                    
                }           
    
                String genesymbol = g.getDbObjectSymbol().trim();
                AnnotatedProtein ap = geneSymbol2AnnotatedProtein.get(genesymbol);
                if(ap == null) {
                    ap = new AnnotatedProtein(genesymbol);
                    geneSymbol2AnnotatedProtein.put(genesymbol, ap);
                    
                }
                //ap.addAnnotation(gt);
                addAnnotations(ap, gt, allgoid2Term);
                counter++;
            }
        }
        //System.out.println("Total number of associations processed: " + counter);
        //System.out.println("Total number of unknown annotations: " + unknowngoid);
        //System.out.println("Total number of unknown goids: " + unknow.size());
        return ipi2AnnotatedProteins;
    } 
    
    private static void addAnnotations(AnnotatedProtein ap, GoTerm gt, HashMap<String, GoTerm> allgoid2Term) {
        if(ap == null || gt == null) return; 
        //System.out.prinlnt("Add Annotations for protein " + ap + "\tfor go term " + gt.getGoId());
        ap.addAnnotation(gt);
        for(Iterator<String> it = gt.getIsAs(); it.hasNext();) {
            String goid = it.next();
            GoTerm goterm = allgoid2Term.get(goid);
            if(goterm != null) {
                ap.addAnnotation(goterm);
                addAnnotations(ap, goterm, allgoid2Term);
            } else {
                 System.out.println("Unknown GoId: " + goid);
            }
        }
    } 

    private static HashMap<String, GoTerm> readExtGeneOntology(String extgofile) throws IOException {

        HashMap<String, GoTerm> goid2term = new HashMap<String, GoTerm>(1000000);
        //GeneOntologyExtReader goer = new GeneOntologyExtReader("/home/taoxu/taoxu_on_data/projects/goa/20100326/gene_ontology_ext.obo");
        GeneOntologyExtReader goer = new GeneOntologyExtReader(extgofile);
        ArrayList<GoTerm> terms = goer.getExtGoTerms();
        Iterator<GoTerm> it = terms.iterator();
        while(it.hasNext()) {
            GoTerm goterm = it.next();
            goid2term.put(goterm.getGoId(), goterm);
            
        }
        return goid2term;
    } 
}
