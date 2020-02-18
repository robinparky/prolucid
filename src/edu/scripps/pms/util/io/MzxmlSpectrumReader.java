/**
 * @file MzxmlSpectrumReader.java
 * This is the source file for edu.scripps.pms.util.spectrum.MzxmlSpectrumReader
 * @author Tao Xu
 * @date $Date
 */



package edu.scripps.pms.util.io;

import java.io.*;
import java.util.List;
import java.util.Stack;
import java.util.Set;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import edu.scripps.pms.util.spectrum.*;
//import edu.scripps.pms.mspid.MassSpecConstants;
import edu.scripps.pms.util.MZXmlHandler;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
//import org.jdom.Namespace;
import org.jdom.input.SAXBuilder;


public class MzxmlSpectrumReader {

    //public static final String MZXML = "mzxml";
    // delimiter of m2z and intensity
    public static final String MZINTENSITYDELIMITER = " ";
    private RandomAccessFile raf; // reads from the mzxml file
    // the format of ms data file, e.g., ms1, ms2, dta
    private String msFileName;
    private HashMap<Integer, Long> scanNum2Position = null;
    private ArrayList<Integer> scanNums = null;
    private long endIndex;
    private int msLevel = 0;
    private long fileLength;
    public static final String ENDSCANTAG = "</scan";
    public static final String BEGINSCANTAG = "<scan";
    public static final String ENDMSRUNTAG = "</msRun";
    public static final float MASSPROTON = 1.007276466f;
    
    // this reader will read only spectra of  specified ms level 
    public MzxmlSpectrumReader (String msFileName, int msLevel) throws IOException {
        this.msLevel = msLevel;
        init(msFileName);
    }
   
    
    public MzxmlSpectrumReader (String msFileName) throws IOException {
        init(msFileName);
    }

    private void init(String msFileName) throws IOException {
        this.msFileName = msFileName;
        raf = new RandomAccessFile(msFileName, "r"); 

        //System.out.println("Begin reading Indeces");
        getIndices(); 
        //System.out.println("Finish reading Indeces");
    }

    public MzxmlPeakList scanNum2PeakList(Integer scanNum, boolean withChildren) 
                                                throws IOException, JDOMException, Exception {
        return getPeakList(scan2Position(scanNum), withChildren);
    } 
    public MzxmlPeakList scanNum2PeakList(int scanNum, boolean withChildren) 
                                          throws IOException, JDOMException, Exception {
        return getPeakList(scan2Position(scanNum), withChildren);
    } 

    // this function has to read through the spectrum file to count the number of spectra,
    // so don't use it unless you have to.
    // subclass may need to override this function
    public int getNumSpectra() throws IOException { 
        return scanNum2Position.size();
    }
    private Stack<Integer> getScanNumStack() {
        Stack<Integer> s = new Stack<Integer>();
        for(int i = scanNums.size()-1; i >=0; i--) {
//System.out.print("ScanNum: " + scanNums.get(i).intValue());
            s.push(scanNums.get(i)); 
        }
        return s;
    }
 
    // return -1 if the scan number is invalid
    private long scan2Position(Integer scan) {
         Long index =  scanNum2Position.get(scan);
        if(index != null) {
            return index.longValue();
        } else {
            return -1;
        }
    }

    // get the index of the start position of the given scan number
    // return -1 if the scan number is invalid
    private long scan2Position(int scan) {

        Long index =  scanNum2Position.get(new Integer(scan));
        if(index != null) {
            return index.longValue();
        } else {
            return -1;
        }
    }
    // return the position where the index starts in the xml file.  return -1 if not found
    private void getIndices() throws IOException {
        fileLength = raf.length();
        endIndex = fileLength - 1;
        raf.seek(fileLength - 500);
        String line = "";
        long indexStartPosition = -1;
        // get the start position of indices 
        while((line = raf.readLine()) != null) {
            if(line.trim().startsWith("<indexOffset")) {
                String indexElement = line.split(">")[1];
                int pos = indexElement.indexOf("<");
                if(pos > 0) {
                    indexElement = indexElement.substring(0, pos);
                }
                 indexStartPosition = Long.parseLong(indexElement); 
            }
        }
         
        scanNum2Position = new HashMap<Integer, Long>(100000);
        if(indexStartPosition < 0) {
            return;
        }
        // get all scan numbers and indices
        raf.seek(indexStartPosition);
        while((line = raf.readLine()) != null) {
            line = line.trim();    
            if(line.startsWith("<offset id=")) {
               
                String [] arr = line.split(">");
                String positionStr = arr[1].substring(0, arr[1].indexOf("<"));
                if(positionStr != null && !positionStr.equals("")) {
                    long pos = Long.parseLong(positionStr);
                    //System.out.println("in mzXML reader: " + line + "\tarr[0]: " + arr[0]);
 
                     // the following are index from Aginent mzXML file:w

                     //<index name="scan" >
                     //    <offset id="1" >1003</offset>
                     //    <offset id="2" >1322</offset>
                     //    <offset id="3" >1640</offset>

                    int scan = Integer.parseInt(arr[0].substring(arr[0].indexOf("\"")+1, arr[0].lastIndexOf("\"")));
                    scanNum2Position.put(new Integer(scan), new Long(pos));
                }
            }
        }
        Set<Integer> scans = scanNum2Position.keySet(); 
        scanNums = new ArrayList<Integer>(scans.size());
        for(Integer i : scans) {
            scanNums.add(i);
        }
        Collections.sort(scanNums);

    }
    public ArrayList<Integer> getScanNums() {
        return scanNums;
    } 

    // for testing 
    public static void main(String args[]) throws Exception {
        MzxmlSpectrumReader msr = new MzxmlSpectrumReader(args[0]);
        //System.out.println("Number of scans: " + scan2Position.size());
        int numSpectra = 0;
        int mslevel = Integer.parseInt(args[1]);
        Iterator<MzxmlPeakList> it = null; 

        if(mslevel > 0) {
            it = msr.getSpectra(mslevel);
        } else {
            it = msr.getSpectra();
        }
/*
        while(it.hasNext()) {
            MzxmlPeakList mpl = it.next();
            numSpectra++;
            //System.out.println("ScanNum: " + mpl.getLoscan() + "\tMsLevle: " + mpl.getMsLevel());
        }

        System.out.println("numSpectra: " + numSpectra);
        System.out.println("Start to retrieve peaklist without children");
        for(Integer scan : msr.getScanNums()) {
            MzxmlPeakList peaks = msr.scanNum2PeakList(scan, false); 
            System.out.println("scan: " + peaks.getLoscan() + "\tNumChildren: " + peaks.getNumChildSpectrum()); 
        } 

        System.out.println("Start to retrieve peaklist with children");
        for(Integer scan : msr.getScanNums()) {
            MzxmlPeakList peaks = msr.scanNum2PeakList(scan, true); 
            
            System.out.println("scan: " + peaks.getLoscan() + "\tNumChildren: " + peaks.getNumChildSpectrum()); 
            //System.out.println(peaks.getSpectrumWithoutHlines());
        } 
*/
        int i = 0;
        //for(it = msr.getSpectraWithChildren(); it.hasNext();) {
        while(it.hasNext()) {
            //msr.printMzxmlPeakList(it.next());
            MzxmlPeakList mpl = it.next();
            System.out.print("in mzxmlspectrumReader, spectra number " + ++i); 
            System.out.println("\nscanNumber: " + mpl.getLoscan() + "\tmsLevel: " + mpl.getMsLevel() + "\tnumChildren: " + mpl.getNumChildSpectrum());
            //msr.printMzxmlPeakList(it.next());      
       
        }         
    }
    // for testing
    private void printMzxmlPeakList(MzxmlPeakList mpl) {
        System.out.println("scanNumber: " + mpl.getLoscan() + "\tmsLevel: " + mpl.getMsLevel() + "\tnumChildren: " + mpl.getNumChildSpectrum());
        for(Iterator<Peak> it = mpl.getPeaks(); it.hasNext();) {
            Peak p = it.next();
            System.out.println("mass: " + p.getM2z() + "\tIntensity: " + p.getIntensity()); 
        }
        for(Iterator<MzxmlPeakList> it = mpl.getChildren(); it.hasNext();) {
             printMzxmlPeakList(it.next());
        }
    }
    public Iterator <MzxmlPeakList> getSpectraWithChildren() throws IOException, JDOMException {
           
        return new Iterator<MzxmlPeakList>() {
            String lastLine = "";
            MzxmlPeakList nextPeakList = null;
            BufferedReader br = new BufferedReader(new FileReader(msFileName), 40960);
            { readFirstPeakList(); };
            private void readFirstPeakList() throws IOException, JDOMException {
                while((lastLine = br.readLine()) != null && !lastLine.trim().startsWith("<scan")){};
//System.out.println(lastLine);
                nextPeakList = retrievePeakList();
                
            }
            private MzxmlPeakList retrievePeakList() throws IOException, JDOMException {
                int numOpenTag = 1;
                int numCloseTag = 0; 
                StringBuffer sb = new StringBuffer(150000);
                if(!lastLine.trim().startsWith(BEGINSCANTAG)) {
                    return null;
                } 
                sb.append(lastLine); // read first open tag line 
                while(numOpenTag != numCloseTag && (lastLine = br.readLine()) != null) {
                    String trimmedline = lastLine.trim();
                    if(trimmedline.startsWith("<scan")) {
                       numOpenTag++; 
                    } else if(trimmedline.startsWith("</scan") || trimmedline.endsWith("/scan>")) {
                        numCloseTag++;
                    }
             
                    sb.append(lastLine);

                }
                int t = sb.toString().indexOf("-int\"\""); 
                if (t > 0) {
                    //System.out.println("here");
                    //remove extra quote
                    sb.deleteCharAt(t+5);
                }
                lastLine = br.readLine();
                InputStream is = new StringBufferInputStream(sb.toString()); 
//System.out.println("Sb.tostring(): " + sb.toString());
                Document doc = new SAXBuilder().build(is);
                return readPeakList(doc.getRootElement(), true); // decode
               
            }
            public boolean hasNext() {
                if(nextPeakList == null) { 
                    try {
                        br.close(); // close the inputstream
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
                return  nextPeakList != null; 
            }

            public MzxmlPeakList next() {

                MzxmlPeakList peaks = null;
                try {
                    peaks = nextPeakList;
                    nextPeakList = retrievePeakList();
                } catch (IOException e) {
                    e.printStackTrace();
                } catch (JDOMException je) {
                    je.printStackTrace();
                }
                return peaks;
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported");
            }
        };
    }
   
    // return an Iterator of all MzxmlPeakList for all spectra, regardless of the ms level 
    public Iterator <MzxmlPeakList> getSpectra() throws IOException, JDOMException, Exception {
        return getSpectra(0);
    }

    public Iterator <MzxmlPeakList> getScanNums(final int msLevel) throws IOException, JDOMException, Exception {
           
        return new Iterator<MzxmlPeakList>() {
            String lastLine = "";
            private int mslevel = msLevel;
            MzxmlPeakList nextPeakList = null;
            BufferedReader br = new BufferedReader(new FileReader(msFileName), 40960);
            { readFirstPeakList(mslevel); };
            private void readFirstPeakList(int mslevel) throws IOException, JDOMException, Exception {
                while((lastLine = br.readLine()) != null && !lastLine.trim().startsWith("<scan")){};
//System.out.println(lastLine);
                nextPeakList = retrievePeakList(mslevel);
                
            }
            private MzxmlPeakList retrievePeakList(int mslevel) throws IOException, JDOMException, Exception {

                StringBuffer sb = new StringBuffer(40960);
                if(lastLine == null || !lastLine.trim().startsWith(BEGINSCANTAG)) {
                    return null;
                }
                sb.append(lastLine);

		System.out.println(lastLine);
                while((lastLine = br.readLine()) != null && !lastLine.trim().startsWith(BEGINSCANTAG) 
                         && !lastLine.trim().startsWith(ENDSCANTAG)){
                    sb.append(lastLine);
                }
                //sb.append(ENDSCANTAG+">");
                String trimmed = sb.toString().trim();
                if(trimmed.endsWith("/mzXML>")) {
                    int index = trimmed.lastIndexOf("/scan>");
                    sb = new StringBuffer(trimmed.substring(0, index+6));
                }
                //if(!trimmed.endsWith("/scan>")) {
                if(!sb.toString().endsWith("/scan>")) {
                    sb.append("</scan>");
                }
                
//System.out.println("Sb.tostring(): " + sb.toString());
                // remove unused endscantags
                while(lastLine != null && !lastLine.trim().startsWith(BEGINSCANTAG)) {
                    lastLine = br.readLine(); 
                }     

                InputStream is = new StringBufferInputStream(sb.toString()); 
//System.out.println("Sb.tostring(): " + sb.toString());
                Document doc = new SAXBuilder().build(is);
                MzxmlPeakList peaklist = readPeakList(doc.getRootElement(), false);
                if(mslevel != 0 && peaklist.getMsLevel() != mslevel) {
                    return retrievePeakList(mslevel);
                } else {
                    //MZXmlHandler.decode32(peaklist.getEncodedM2zAndIntensities(), peaklist);
                    MZXmlHandler.decode(peaklist);
                    return peaklist;
                }
               
            }
            public boolean hasNext() {
                if(nextPeakList == null) { 
                    try {
                        br.close(); // close the inputstream
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
                return  nextPeakList != null; 
            }

            public MzxmlPeakList next() {

                MzxmlPeakList peaks = null;
                try {
                    peaks = nextPeakList;
                    nextPeakList = retrievePeakList(mslevel);
                } catch (IOException e) {
                    e.printStackTrace();
                } catch (JDOMException je) {
                    je.printStackTrace();
                } catch (Exception je) {
                    je.printStackTrace();
                }
                return peaks;
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported");
            }
        };
    }
    // return a Iteractor of MzxmlPeakList for spectra of given ms level, like 0, 1, 2, 3, and etc.
    // msLevel 0 for all spectra, 1 for ms1, 2 for ms2 and 3 for ms3 
    // child spectrum are not included 
    public Iterator <MzxmlPeakList> getSpectra(final int msLevel) throws IOException, JDOMException, Exception {
           
        return new Iterator<MzxmlPeakList>() {
            String lastLine = "";
            private int mslevel = msLevel;
            MzxmlPeakList nextPeakList = null;
            BufferedReader br = new BufferedReader(new FileReader(msFileName), 40960);
            { readFirstPeakList(mslevel); };
            private void readFirstPeakList(int mslevel) throws IOException, JDOMException, Exception {
                while((lastLine = br.readLine()) != null && !lastLine.trim().startsWith("<scan")){};
//System.out.println(lastLine);
                nextPeakList = retrievePeakList(mslevel);
                
            }
            private MzxmlPeakList retrievePeakList(int mslevel) throws IOException, JDOMException, Exception {

                StringBuffer sb = new StringBuffer(40960);
                if(lastLine == null || !lastLine.trim().startsWith(BEGINSCANTAG)) {
                    return null;
                }
                sb.append(lastLine);
                while((lastLine = br.readLine()) != null && !lastLine.trim().startsWith(BEGINSCANTAG) 
                         && !lastLine.trim().startsWith(ENDSCANTAG)){
                    sb.append(lastLine);
                }
                //sb.append(ENDSCANTAG+">");
                String trimmed = sb.toString().trim();
                if(trimmed.endsWith("/mzXML>")) {
                    int index = trimmed.lastIndexOf("/scan>");
                    sb = new StringBuffer(trimmed.substring(0, index+6));
                }
                //if(!trimmed.endsWith("/scan>")) {
                if(!sb.toString().endsWith("/scan>")) {
                    sb.append("</scan>");
                }
                
                //System.out.println("Sb.tostring(): " + sb.toString());
                int t = sb.toString().indexOf("-int\"\""); 
                if (t > 0) {
                    //System.out.println("here");
                    //remove extra quote
                    sb.deleteCharAt(t+5);
                }

                //System.out.println("Sb.tostring(): " + sb.toString());
                // remove unused endscantags
                while(lastLine != null && !lastLine.trim().startsWith(BEGINSCANTAG)) {
                    lastLine = br.readLine(); 
                }     
//System.out.println(sb);
                InputStream is = new StringBufferInputStream(sb.toString()); 
//System.out.println("Sb.tostring(): " + sb.toString());
                Document doc = new SAXBuilder().build(is);
                MzxmlPeakList peaklist = readPeakList(doc.getRootElement(), false);
                if(mslevel != 0 && peaklist.getMsLevel() != mslevel) {
                    return retrievePeakList(mslevel);
                } else {
                    //MZXmlHandler.decode32(peaklist.getEncodedM2zAndIntensities(), peaklist);
                    MZXmlHandler.decode(peaklist);
                    return peaklist;
                }
               
            }
            public boolean hasNext() {
                if(nextPeakList == null) { 
                    try {
                        br.close(); // close the inputstream
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
                return  nextPeakList != null; 
            }

            public MzxmlPeakList next() {

                MzxmlPeakList peaks = null;
                try {
                    peaks = nextPeakList;
                    nextPeakList = retrievePeakList(mslevel);
                } catch (IOException e) {
                    e.printStackTrace();
                } catch (JDOMException je) {
                    je.printStackTrace();
                } catch (Exception je) {
                    je.printStackTrace();
                }
                return peaks;
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported");
            }
        };
    }
    private MzxmlPeakList getPeakList(long start, boolean withChildren) throws IOException, JDOMException, Exception {
        raf.seek(start);       
        
        int numOpenTag = 1;
        int numCloseTag = 0; 
 //       System.out.println("\n\n\nStart reading one spectrum \n\n");
        String line = "";
        StringBuffer sb = new StringBuffer(150000);
        //    System.out.println(raf.readLine() + "\n\n");
        sb.append(raf.readLine()); // read first open tag line 
        while(numOpenTag != numCloseTag && (line = raf.readLine()) != null) {
            //line = line.trim();
            if(line.trim().startsWith("<scan")) {
               if(!withChildren) {
                    String trimmed = sb.toString().trim();
                    if(trimmed.endsWith("/mzXML>")) {
                        int index = trimmed.lastIndexOf("/scan>");
                        sb = new StringBuffer(trimmed.substring(0, index+6));
                        break;
                    }
                    if(!sb.toString().endsWith("/scan>")) {
                        sb.append("</scan>");
                    }
            
                   //sb.append("</scan>");
                   break;
               }
               numOpenTag++; 
//System.out.println("open tag found");
            } else if(line.trim().startsWith("</scan")) {
                numCloseTag++;
//System.out.println("close tag found");
            }
             
            sb.append(line);

        }
  //      System.out.println("\n\n\nFinished reading one spectrum \n\n");
//System.out.println(sb.toString());
        int t = sb.toString().indexOf("-int\"\""); 
        if (t > 0) {
            //System.out.println("here");
            //remove extra quote
            sb.deleteCharAt(t+5);
        }

        InputStream is = new StringBufferInputStream(sb.toString()); 
        
        Document doc = new SAXBuilder().build(is);

        //System.out.println(root.getAttributeValue("num"));
        MzxmlPeakList peaklist = readPeakList(doc.getRootElement(), true); // decode
        //MZXmlHandler.decode32(peaklist.getEncodedM2zAndIntensities(), peaklist);
        MZXmlHandler.decode(peaklist);
        is.close();
         
        return peaklist; 
    }
    //protected MzxmlPeakList readPeakList(Document doc) {

    protected MzxmlPeakList readPeakList(Element root, boolean decode) {
        MzxmlPeakList peakList = new MzxmlPeakList();
        //Element root = doc.getRootElement();
        String num = root.getAttributeValue("num");
        String msLevel = root.getAttributeValue("msLevel");
        String peaksCount = root.getAttributeValue("peaksCount");
        String retentionTime = root.getAttributeValue("retentionTime");
        String basePeakMz = root.getAttributeValue("basePeakMz");
        String basePeakIntensity = root.getAttributeValue("basePeakIntensity");

        int scannum = num == null? 0 : Integer.parseInt(num);
        peakList.setLoscan(scannum);
        peakList.setHiscan(scannum);

        int mslevel = msLevel == null? 0 : Integer.parseInt(msLevel); 
        peakList.setMsLevel(mslevel);
        if(retentionTime != null) {
            if(retentionTime.startsWith("PT")) {
                retentionTime = retentionTime.substring(2, retentionTime.length()-1);
            }
        }
        double retTime = retentionTime == null? 0 : Double.parseDouble(retentionTime);
        peakList.setRetentionTime(retTime);
 
        double baseMz = basePeakMz == null? 0 : Double.parseDouble(basePeakMz);
        double baseIntensity = basePeakIntensity == null? 0 : Double.parseDouble(basePeakIntensity);
        peakList.setBasePeakM2z(baseMz);
        peakList.setBasePeakIntensity(baseIntensity);

        Element peaks = root.getChild("peaks");
        int precision = Integer.parseInt(peaks.getAttributeValue("precision"));
        String encodedMzAndIntensity = peaks.getTextTrim();
        peakList.setEncodedM2zAndIntensities(encodedMzAndIntensity, precision); 

        Element precursor = root.getChild("precursorMz");

        String precursorCharge = null;

        if(precursor != null) {
            String precursorScanNum = precursor.getAttributeValue("precursorScanNum");
            int precScan = precursorScanNum == null? 0 : Integer.parseInt(precursorScanNum);
            peakList.setPrecursorScan(precScan);
            float prcMz = Float.parseFloat(precursor.getTextTrim());
            peakList.setPrecursorMass(prcMz);
            peakList.setPrecursorIntensity(Double.parseDouble(precursor.getTextTrim()));
            precursorCharge = precursor.getAttributeValue("precursorCharge");
            if(precursorCharge == null) {
                peakList.addZline(new Zline(1, prcMz)); 
                peakList.addZline(new Zline(2, prcMz*2-MASSPROTON)); 
                peakList.addZline(new Zline(3, prcMz*3-2*MASSPROTON)); 
                //peakList.addZline(new Zline(4, prcMz*4-3*MASSPROTON)); 
            } else {
                int precCharge = Integer.parseInt(precursorCharge);
                peakList.addZline(new Zline(precCharge, prcMz*precCharge-(precCharge-1)*MASSPROTON)); 

            }
        }
        List<Element> childrenScans = root.getChildren("scan"); 
        if(childrenScans != null) {
            Iterator<Element> children = childrenScans.iterator();
            while(children.hasNext()) {
                peakList.addChildSpectrum(readPeakList(children.next(), decode)); 
            }
        }
        if(decode) {
            MZXmlHandler.decode32(peakList.getEncodedM2zAndIntensities(), peakList);
        }

        if(precursorCharge == null) {
            // need to determine if it is charge 1 or 2 and 3 
            peakList.reAssignCharge123();     
        }
        return peakList;
    }
    public void closeDataFile() throws IOException {
        raf.close();
    }

}
