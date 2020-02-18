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
 * @version $Id: SqtStatistics.java,v 1.5 2013/02/22 17:39:51 taoxu Exp $
 */
public class SqtStatistics
{
    private String fileName;
    private BufferedReader br;
    private String lastLine = "";
    private boolean isValidPeptide=true;

    public static void main(String args[]) throws Exception
    {
         SqtStatistics reader = new SqtStatistics( args[0] );
         int [] frequence = new int[100];
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
//System.out.println("calMZ: " + mLine.getCalMZ() + "\tprcMass: " + peptide.getCalMZ());
                 if (rank == 1) {
//System.out.println("clacMass: " + calMz);
                     double diff = calMz - prcMass;
                     if(diff > 0.5) {
                         System.out.println(mLine.getSequence());
                     }
                     diff *= 10;
                     int intDiff = (int)Math.round(diff)+40;
try{
                     frequence[intDiff]++;
} catch(Exception e) {
System.err.println("Problem found in M Line: " + mLine.getMLine() + "\t in S Line: " + peptide.getSLine());
e.printStackTrace();
    throw e;
} 
                     break;
                 }
                 
             }
         }
         for (int i = 0; i < frequence.length; i++) {
             System.out.println((i-40)/10.0 + "\t" + frequence[i]);
        }
    }

    public SqtStatistics(String fileName) throws IOException
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
        while( lastLine != null && lastLine.equals("")) {
            lastLine = br.readLine();
        } //read empty line
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
