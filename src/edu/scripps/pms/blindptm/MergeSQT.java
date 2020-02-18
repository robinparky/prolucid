package edu.scripps.pms.blindptm;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintStream;
import java.util.ArrayList;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author rampuria
 */
public class MergeSQT {
    
    public static void main(String [] args){
        String foldername ="/home/rampuria/prolucid/completed/";
        MergeSQT s = new MergeSQT();
        s.mergeSQT(foldername);
    }
    
    public void mergeSQT(String foldername){
        File outfile = new File(foldername+File.separator+"mergedfile.sqt");
                if(outfile.exists()){
                    outfile.delete();
                }
        ArrayList<String> sqtfiles = new ArrayList<String>();
        File myfile = new File(foldername);
        if(myfile.exists()) {
            if(myfile.isDirectory()) {
                String [] allfiles = myfile.list();
                for(int i = 0; i < allfiles.length; i++) {
                    String file = allfiles[i];
                    if(file != null) {
                        if(file.endsWith(".sqt")) {
                            String ms2filename = foldername + File.separator + file;
                            sqtfiles.add(ms2filename);
                            //System.out.println("adding " + ms2filename);
                        }
                    }
                }
            } else {
                sqtfiles.add(foldername);
            }        
        }

        BufferedReader br = null;
        boolean isheader = false;
        try{   
        BufferedWriter out = new BufferedWriter(new FileWriter(foldername+File.separator+"mergedfile.sqt"));
        StringBuffer hb =new StringBuffer();
        for(String sqtfile:sqtfiles){
            File f = new File(sqtfile);
                br = new BufferedReader(new FileReader(f));
                String eachLine = null;
                if(hb.length() != 0){
                        isheader=true;
                    }
                while((eachLine = br.readLine()) != null){
                    if(eachLine.startsWith("H\t") && isheader ==false){
                        hb.append(eachLine);
                        out.write(eachLine+"\n");
                    }
                    else if(!eachLine.startsWith("H\t")){
                        //sb.append(eachLine).append("\n");
                        out.write(eachLine+"\n");
                    }
                }
                br.close();
            }
        out.close();
                
            }catch(Exception e){
                e.printStackTrace();
            }
            
        
        
    }
    
    
    
}
