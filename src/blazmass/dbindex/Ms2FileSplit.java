/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package blazmass.dbindex;

import edu.scripps.pms.util.io.FileUtils;
import java.awt.font.NumericShaper;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.Range;

/**
 *
 * @author rampuria
 */
public class Ms2FileSplit {
    private int max = 6000000;
    
    
    public static void main(String [] args){
        //String file ="/home/rampuria/data/newprolucid/20141101_HeLa_1ug_BEH60_140min_35ms_IC_DE5_5e3_1.ms2";
        String file = "/data/2/rpark/ip2_data/cbamberg/Fibulin_3/wt_2_mono_2016_07_11_17_83427/search/projects2016_07_12_10_100050/temp/test.ms2";
        Ms2FileSplit m = new Ms2FileSplit();
        m.splitfile(file,1000,6000000);
        
    }
    
    public void splitfile(String ms2file,int splitnum,int maxrange){
        this.max = maxrange;
       
        File ms2File = new File(ms2file);
        String ms2name = FilenameUtils.removeExtension(ms2File.getName());
        int totalprecmass = this.max;
        List<Integer> range = new ArrayList<>();
        range.add(0);
        for(int i=0;i<totalprecmass;){
            i=i+totalprecmass/splitnum;
            range.add(i);
        }
        BufferedReader br = null;
       try{
           File foldercheck = new File(ms2File.getParent()+File.separator+"splitms2");
           if(foldercheck.exists() && foldercheck.isDirectory()){
                for(File f : foldercheck.listFiles()){
                    f.delete();
                }
                foldercheck.delete();
           }
           new File(ms2File.getParent()+File.separator+"splitms2").mkdir();
           br = new BufferedReader(new FileReader(new File(ms2file)));
           String  eachLine = br.readLine();
           List<BufferedWriter> bout = new ArrayList<>();
           for(int i =0;i<splitnum;i++){
               bout.add(new BufferedWriter(new FileWriter(ms2File.getParent()+File.separator+"splitms2"+File.separator+ms2name+"_"+i+".ms2")));
           }
           for(int i =0;i<splitnum;i++){
               bout.get(i).write("H\tRANGE\t"+(range.get(i)-10000)+"\t"+(range.get(i+1)+10000)+"\n");
           }
           while(eachLine != null){
               if(eachLine.startsWith("S\t")){
                   StringBuffer slinebuffer = new StringBuffer();
                   slinebuffer.append(eachLine+"\n");
                   StringBuffer iLinebuffer = new StringBuffer();
                   StringBuffer massbuffer = new StringBuffer();
                   List<String> zline = new ArrayList<>();
                   String [] words = new String[2];
                   List<Double> precmass = new ArrayList<>();
                   //String [] words = eachLine.split("\t");
                   while((eachLine=br.readLine()) != null && !(eachLine.startsWith("S\t"))){
                       
                       if(eachLine.startsWith("Z\t")){
                          // zlineBuffer.append(eachLine+"\n");
                           zline.add(eachLine+"\n");
                           words = eachLine.split("\t"); 
                           precmass.add(Double.parseDouble(words[2]));
                       }
                       else if(eachLine.startsWith("I\t")){
                           iLinebuffer.append(eachLine+"\n");
                       }
                       else{
                           massbuffer.append(eachLine+"\n");
                       }
                   }
                   int bufferindex =0;
                   int zlineindex=0;
                   for(int k=0;k<precmass.size();k++){
                   //(int)((Double.parseDouble(words[3])+0.0005)*1000);
                      
                   int massint = (int)((precmass.get(k)+0.0005)*1000);
                   for(int i=0;i<range.size()-1;i++){
                      Range<Integer> r = Range.between(range.get(i), range.get(i+1));
                      if(r.contains(massint)){
                          bufferindex=i;
                          zlineindex=k;
                          break;
                      }
                      else{
                          continue;
                      }
                   }
                   
                   bout.get(bufferindex).write(slinebuffer.toString()+iLinebuffer.toString()+zline.get(zlineindex)+massbuffer);
                   }
                  /* StringBuffer sbtotal = buffer[bufferindex];
                   if(sbtotal == null){
                       sbtotal=sb;
                       buffer[bufferindex]=sbtotal;
                   }else{
                        buffer[bufferindex]=sbtotal.append(sb); 
                   }*/
                   continue;
               }
               
               eachLine=br.readLine();
                       
           }
           for(BufferedWriter out:bout){
               out.close();
           }
           /*
           deleting null files
           */
           File folder = new File(ms2File.getParent()+File.separator+"splitms2"+File.separator);
           if(folder.exists()) {
            if(folder.isDirectory()) {
          String [] allfiles = folder.list();
          for(String f :allfiles){
              File check = new File(folder+File.separator+f);
              long length = check.length();
              if(length < 100){
                 check.delete();
              }
              
          }
            }
           }
           
           /*new File(ms2File.getParent()+File.separator+"splitms2").mkdir();
           for(int i=0;i<buffer.length;i++){
               if(buffer[i] == null){
                   continue;
               }
               else{
                   File f = new File(ms2File.getParent()+File.separator+"splitms2"+File.separator+ms2name+"_"+i+".ms2");
                   BufferedWriter out = new BufferedWriter(new FileWriter(f));
                   out.write(buffer[i].toString());
                   out.close();
               }
           }*/
           
           
       }catch (Exception e){
           e.printStackTrace();
       }
    }
    
}
