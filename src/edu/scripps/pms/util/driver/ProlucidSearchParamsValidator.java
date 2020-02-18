/**
 * @file ProlucidSearchParamsValidator.java
 * This is the source file for edu.scripps.pms.util.driver.ProlucidSearchParamsValidator
 * @author Tao Xu
 * @date $Date: 2007/09/10 20:25:41 $
 */
package edu.scripps.pms.util.driver;
import edu.scripps.pms.mspid.SearchParams;

public class ProlucidSearchParamsValidator {

    public static void main(String args[]) throws Exception {

        try {
            String paramfile = "search.xml";
            if(args.length > 0) {
                paramfile = args[0];
            }
            SearchParams sp = new SearchParams(paramfile);
            
        } catch(Exception e) {
            e.printStackTrace();
            System.err.println("\n\n!!! Something is wrong in your ProLuCID search parameter file !!!");
            System.err.println("!!! If you cannot figure out what the problem is, please talk to Tao !!!");
            System.exit(1); 
        }
        
        System.err.println("\n\nEverething seems fine in your ProLuCID search parameter file.  Please go ahead to submit your jobs.");
        System.exit(0);

    }
}


