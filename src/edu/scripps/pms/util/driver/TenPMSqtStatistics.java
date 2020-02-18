import edu.scripps.pms.util.io.SQTParser;
import java.io.BufferedReader;
import java.util.Iterator;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.FileReader;
import java.io.RandomAccessFile;
import edu.scripps.pms.util.sqt.SQTPeptide;
import edu.scripps.pms.util.sqt.MLine;

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
 * @version $Id: TenPMSqtStatistics.java,v 1.3 2007/03/08 21:46:01 taoxu Exp $
 */

import edu.scripps.pms.mspid.MassCalculator;

public class TenPMSqtStatistics
{
    private String fileName;
    private BufferedReader br;
    private String lastLine = "";
    private boolean isValidPeptide=true;
    private static String USAGE = "Usage: java TenPMSqtStatistics sqtfile mono/avg(boolean, true for avg)";
    public static void main(String args[]) throws IOException
    {
         TenPMSqtStatistics reader = new TenPMSqtStatistics( args[0] );
         final int numShift = 50;
         final int freqLen = 100;
         boolean isAvg = false;
         try {
             isAvg = Boolean.parseBoolean(args[1]);
         } catch(Exception e) {
             System.err.println(USAGE);
             System.exit(1);
         }
         int [] frequence = new int[freqLen];
         int [] badFreq = new int[freqLen];
         SQTPeptide peptide;
         MLine mLine;
         String temp;

         int numPassed = 0;
         int numRealHitsPassed = 0;
         int numDiffTooGreat = 0;
         int numDiffTooGreatBad = 0;
         int totalRealHits = 0;
         int numBadPassed = 0;
         int numSpectra = 0;
         for(Iterator<SQTPeptide> itr = reader.getSQTPeptide(); itr.hasNext(); )
         {
             numSpectra++;
             peptide = itr.next();
             int chargeState = Integer.parseInt(peptide.getChargeState());
             double prcMass = Double.parseDouble(peptide.getCalMZ().trim());
             
             for(Iterator<MLine> mItr = peptide.getMLine(); mItr.hasNext(); )
             {
                 mLine = mItr.next();
                 // String lLine = mLine.getLLine().next();
                 double calMz = Double.parseDouble(mLine.getCalMZ().trim()); 
                 int rank = Integer.parseInt(mLine.getXcorrRank().trim());
                 double deltaCN = Double.parseDouble(mLine.getDeltCN());
                 if (rank == 1 && isPassed(mLine, chargeState)) {
                     numPassed++;

                     if(isAvg) {
                         String origseq = mLine.getSequence().trim();
                         String seq = origseq.substring(2, origseq.length()-2);
                         // replace the calculated mass with mono mass
                         calMz = MassCalculator.getMonoMass(seq, chargeState);  
                   
//System.out.println("calMz: " + calMz  + "\tnewCalcMz: " + newCalMz + "\tprcMass: " + prcMass + "\tavgMass: " + MassCalculator.getAvgMass(seq, chargeState));
//System.out.println("torigSeq: " + origseq + "\tseq: " + seq);
                     }
                     double diff = calMz - prcMass;
                     boolean isTooBig = diff > 0.4;
                     diff *= 10;
                     int intDiff = (int)Math.round(diff)+numShift;
                     
                     if(isRealHit(mLine)) {
                         numRealHitsPassed++; 
                         frequence[intDiff]++;
                         if(isTooBig) {
                             numDiffTooGreat++;
                         }
                     } else { // passed but not real hits
                         badFreq[intDiff]++;
                         if(isTooBig) {
                             numDiffTooGreatBad++;
                         }
                     } 
                 }
                 if( rank == 1 && isRealHit(mLine)) {
                     totalRealHits++;
                 }

                 break;
                 
             }
         }
         for (int i = 0; i < frequence.length; i++) {
             System.out.println((i-numShift)/10.0 + "\t" + frequence[i] + "\t" + badFreq[i]);
        }
        System.out.println("NumSpectra: " + numSpectra + "\tNumPassed: " + numPassed + "\ttotalRealHits: " + totalRealHits + "\tNumRealHitPassed: " + numRealHitsPassed + "\tNum m2z diff too great: " + numDiffTooGreat + "\tnumDiffTooGreatBad: " + numDiffTooGreatBad);
    }

    // return true if the peptide hit is a contaminant or 10PM peptide
    public static boolean isRealHit(MLine m) {

        boolean isRealHit = false;
        for(Iterator<String> locus = m.getLLine(); locus.hasNext();) {
            String l = locus.next(); // get a locus
            if(l.startsWith("10PM") || l.startsWith("contaminant")) {
                isRealHit = true;
                break;
            }
        }
        return isRealHit;
    } 
    // return true is this peptide hit pass the deltCN and Xcorr theshholds 
    public static boolean isPassed(MLine m, int chargeState) {
        boolean isTrueHit = false;
        double deltaCN = Double.parseDouble(m.getDeltCN());
        double xcorr = Double.parseDouble(m.getXcorr());
        if(deltaCN >= 0.08) {
            switch (chargeState) {
                case 1: isTrueHit = (xcorr >= 1.8)? true : false; break;
                case 2: isTrueHit = (xcorr >= 2.5)? true : false; break; 
                case 3: isTrueHit = (xcorr >= 3.5)? true : false; break; 
            }
        }

        return isTrueHit;
    }
    
    public TenPMSqtStatistics(String fileName) throws IOException
    {
        this.fileName = fileName;
//        this.reader = new FileReader(fileName);
        init();
    }

    private void init() throws IOException
    {
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
                    double temp = 0;

                    list.clear();
                   
                    // read M lines 
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
                            try {
                                temp = Double.parseDouble(mLine.getDeltCN().trim());                        
                            } catch(Exception e) {
                                System.out.println("lastLine: " + lastLine);
                                e.printStackTrace();
                            }
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

                    if(lastLine==null || lastLine.startsWith("S\t")) {
                        break;
                    }
                }
                if(lastLine != null && lastLine.equals("")) {
                    lastLine = br.readLine();
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
