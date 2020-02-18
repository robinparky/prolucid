
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

import edu.scripps.pms.pmsdb.*;
import edu.scripps.pms.go.AnnotatedProtein;
import edu.scripps.pms.go.GoAssociationReader;
import edu.scripps.pms.go.GoAssociation;
import edu.scripps.pms.go.GoTerm;

import java.util.*;
import java.io.*;

public class GetIpiGoAnnotations {
    public static final String USAGE = "\n\n--- java GetIpiGoAnnotations accessionFileName ---\n";
    public static final String gotermfile = "/data/3/taoxu/projects/go/200904/go_200904-assocdb-tables/term.txt";
    //public static final String annotationfile = "/data/6/lliao/GO_related/gene_association.goa_mouse";
    public static final String annotationfile = "/home/taoxu/taoxu_on_data/projects/go/200904/gene_association.goa_mouse";

    public static void main(String args[]) {
        try {
            //String fileName = args[0];
            
            HashMap<String, GoTerm> goTermMap = getGoTermMap(gotermfile);
            HashMap<String, AnnotatedProtein> ipi2AnnotatedProteins = readAnnotations(annotationfile, goTermMap);
            System.out.println("gotermfile: " + gotermfile);
            System.out.println("associationfile: " + annotationfile);
            System.out.println("IPI\tCellular Component\tMolecular Function\tBiological Process");
            for(Iterator<String> it = ipi2AnnotatedProteins.keySet().iterator(); it.hasNext();) {
                System.out.println(ipi2AnnotatedProteins.get(it.next()).output());
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println(USAGE);
        }
    }
    
    private static HashMap<String, AnnotatedProtein> readAnnotations(String file, HashMap<String, GoTerm> gotermmap) throws IOException  {

        HashMap<String, AnnotatedProtein> ipi2AnnotatedProteins = new HashMap<String, AnnotatedProtein>(100000);
        GoAssociationReader gar = new GoAssociationReader(file);
        Iterator<GoAssociation> it = gar.getGoAssociations();
        int counter = 0;
        //boolean sortByIntensity = true;
        int unknowngoid = 0;
        HashSet<String> unknow = new HashSet<String>(50000);
        while (it.hasNext()) {
            GoAssociation g = it.next();
            //swissProtAccs.add(g.getDbObjectId());
            //ipiAccs.add(g.getDbObjectSynonym());
            String dos = g.getDbObjectSynonym();
            String [] arr = dos.split("\\|");
//System.out.println(dos + "\t" + arr.length);
            String goid = g.getGoId().trim();
            GoTerm gt = gotermmap.get(goid);
            if(gt == null) {
                unknowngoid++; 
                System.out.println("Unknow goid: " + goid); 
                unknow.add(goid); 
                continue;
            }

            for(String s : arr) {
                if(s != null && s.startsWith("IPI")) {
                    
                    AnnotatedProtein ap = ipi2AnnotatedProteins.get(s);
                    if(ap == null) {
                        ap = new AnnotatedProtein(s);
                        ipi2AnnotatedProteins.put(s, ap);
                        
                    }
                    ap.addAnnotation(gt);
                    
                }
                
            }           
 
            counter++;
         }
        System.out.println("Total number of associations processed: " + counter);
        System.out.println("Total number of unknown annotations: " + unknowngoid);
        System.out.println("Total number of unknown goids: " + unknow.size());
        return ipi2AnnotatedProteins;
    } 
     
    private static HashMap<String, GoTerm> getGoTermMap(String fileName) throws IOException {
        HashMap<String, GoTerm> goTermMap = new HashMap<String, GoTerm>(1000000);
        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
        String line = null;
HashSet<String> goids = new HashSet<String>(1000000);
        int numterms = 0;
        while((line = br.readLine()) != null) {
            if(line.equals("") || line.indexOf("\t") == -1) {
                continue;
            }
            
            GoTerm gt = new GoTerm(line); 
            goTermMap.put(gt.getGoId().trim(), gt);
        }
        br.close();
        return goTermMap;
    }
}
