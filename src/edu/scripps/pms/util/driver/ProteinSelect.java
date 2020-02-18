import edu.scripps.pms.util.io.SQTParser;
import java.io.BufferedReader;
import java.util.*;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.FileReader;
import java.io.RandomAccessFile;
import edu.scripps.pms.util.sqt.SQTPeptide;
import edu.scripps.pms.util.sqt.MLine;
import edu.scripps.pms.util.seq.Fasta;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Tao Xu 
 * @version $Id
 */
public class ProteinSelect
{
    private String fileName;
    private BufferedReader br;
    private String lastLine = "";
    private boolean isValidPeptide=true;
    private static HashMap<String, ProteinCount> acc2ProteinCount = new HashMap<String, ProteinCount>();    


    public static String parseAccession(String locus) {
        if(locus == null || !locus.startsWith("Reverse_")) {
            return Fasta.getAccession(locus);
        }
        
        String acc = null;
        if(locus.indexOf("|") > 0) {
            acc = locus.substring(locus.indexOf("_")+1, locus.indexOf("|")); 
        } else {
            acc = locus.substring(locus.indexOf("_")+1); 
        }
        if(locus.startsWith("Reverse_")) {
            //System.out.println(locus + "\t" + acc);
        }
        return Fasta.getAccession(acc); 
    }
    public static void main(String args[]) throws IOException
    {
         ProteinSelect reader = new ProteinSelect(args[0]);
         SQTPeptide peptide;
         MLine mLine;
         String temp;
         for(Iterator<SQTPeptide> itr = reader.getSQTPeptide(); itr.hasNext(); )
         {
             peptide = itr.next();

             double prcMass = Double.parseDouble(peptide.getCalMZ().trim());
             //System.out.print("prcMass: " + prcMass + "\t");
             for(Iterator<MLine> mItr = peptide.getMLine(); mItr.hasNext(); )
             {
                 mLine = mItr.next();
                 double calMz = Double.parseDouble(mLine.getCalMZ().trim()); 
                 int rank = Integer.parseInt(mLine.getXcorrRank().trim());
                 if (rank < 2) {
                     //System.out.println("clacMass: " + calMz);
                     for(Iterator<String> lit = mLine.getLLine(); lit.hasNext();) {
                         String locus = lit.next();
                         String acc = parseAccession(locus);
//System.out.println("locus: " + locus + "\taccession: " + acc);
                         ProteinCount pc = acc2ProteinCount.get(acc);
                         if(pc == null) {
                             pc = new ProteinCount(acc);
                             acc2ProteinCount.put(acc, pc);
                         }
                         pc.addCount(locus, rank);
                     }
                    // break;
                 }
                 
             }
         }
         Collection<ProteinCount> proteinCounts = acc2ProteinCount.values();
         
         List<ProteinCount> counts = new ArrayList<ProteinCount>(proteinCounts.size());
         for(ProteinCount p : proteinCounts) {
             counts.add(p);  
         }
      
         Collections.sort(counts);
         double sigma = calcSigma(counts);
         for(ProteinCount p : counts) {
             System.out.println(p.getAccession() + "\t" + p.getCountDifference() + "\t" + p.getForwardCount() + "\t" + p.getReverseCount());
         }
    }
    
    public static double calcSigma(List<ProteinCount> proteinCounts) {
        double minimum = 100;
         
        for(ProteinCount p : proteinCounts) {
            if(p.getCountDifference() < minimum) {
                minimum = p.getCountDifference();
            }    
        }
        System.out.println("minimum: " + minimum);
        double max = Math.abs(minimum);
        System.out.println("max: " + max);
        double sum = 0;
        int numForMean = 0;
        int n = proteinCounts.size();
        System.out.println("numProteins: " + n);
        for(ProteinCount p : proteinCounts) {
            double countDiff = p.getCountDifference();
            if(countDiff <= max) { 
                sum += p.getCountDifference();
                numForMean++;
            }
        }
        System.out.println("numForMean: " + numForMean);
        double mean = sum/(0.0+numForMean);
        System.out.println("mean: " + mean);
        double sSquare = 0; 
        for(ProteinCount p : proteinCounts) {
            double diff = p.getCountDifference()- mean;
             
            double countDiff = p.getCountDifference();
            if(countDiff <= max) {
                sSquare += diff*diff;
            }
        }
System.out.println("sSQUARE: " + sSquare);
        sSquare /= (numForMean-1);
System.out.println("sSQUARE: " + sSquare);
        double s = Math.sqrt(sSquare);
        System.out.println("s : " + s); 
        return s;
    }
    public ProteinSelect(String fileName) throws IOException
    {
        this.fileName = fileName;
//        this.reader = new FileReader(fileName);
        br = new BufferedReader(new FileReader(fileName));

        readHeader();
    }


    private void readHeader() throws IOException
    {

        while( (lastLine = br.readLine()).startsWith("H\t") );  //read head line
        while( (lastLine = br.readLine()).equals("") ); //read empty line
    }

    public Iterator<SQTPeptide> getSQTPeptide() throws IOException
    {
        return new Iterator<SQTPeptide> ()
        {
            private SQTPeptide peptide;
            private MLine mLine;

            public boolean hasNext()
            {
                return lastLine != null;
            }

            public SQTPeptide next()
            {
                try 
		{
                    peptide = getPeptide();
		    while(!isValidPeptide)
		    {
			isValidPeptide = true;
			peptide = getPeptide();
		    }
                }
                catch (IOException e) {
                    e.printStackTrace();
                }

                return peptide;
            }

            public void remove()
            {
                throw new UnsupportedOperationException("Not supported");
            }

            private SQTPeptide getPeptide() throws IOException
            {

//                String[] strArr = lastLine.split("\t");

                peptide = new SQTPeptide(lastLine);


                java.util.ArrayList<MLine> list = new java.util.ArrayList<MLine>();
                
                while ( ( (lastLine = br.readLine()) != null && !lastLine.equals("")) &&
                        (lastLine.startsWith("M\t") || lastLine.startsWith("L\t")) )
                {
                    //sqt file has wrong deltaCN value.  So, grep the right value from it.                
                    boolean isDeltaCN=false;
                    double deltaCN=0;                
                    double temp;

                    list.clear();
                    
                    while(lastLine!=null && lastLine.startsWith("M\t"))
                    {
			//System.out.println(lastLine);
			/*  If M line is messed up in sqt file, discard it */
			String[] arr = lastLine.split("\t");
			if(arr.length<10)
			{
				isValidPeptide = false;
				break;
			}

                        mLine = new MLine(arr);

                        if(!isDeltaCN)
                        {
                            temp = Double.parseDouble( mLine.getDeltCN() );                        
                        
                            if(temp>0)
                            {
                                deltaCN=temp;                            
                                isDeltaCN=true;
                            }
                        }                                                
                        
                        while ( (lastLine = br.readLine()) != null &&  lastLine.startsWith("L\t") )
                        {
                            mLine.addLLine(lastLine);
                        }
                        
                        list.add(mLine);                        
                    }
                    
                    for(int i=0;i<list.size();i++)
                    {
                        list.get(i).setDeltCN(deltaCN);
                        peptide.addMLine(list.get(i));
                    }
                        
                    //peptide.addMLine(mLine);

                    if(lastLine==null || lastLine.startsWith("S\t"))
                       break;
                }

                return peptide;
            }
        };
    }

    public void close() throws IOException
    {
        br.close();
    }

}
