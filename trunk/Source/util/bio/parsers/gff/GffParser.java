package util.bio.parsers.gff;

import java.util.*;
import java.io.*;
import javax.swing.*;

/**
 * Parses GFF files extracting the annotation according to GFF(2).
 * http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml
 * <seqname>tab<source>tab<feature>tab<start>tab<stop>tab<score>tab
	<strand>tab<frame>tab[space delimited key=value; key=value; attributes]
 * Uses some hashes to count the sources, features, and attributes
 * Generates an array of GffFeature objects, see to learn about what's
 *  used to parse the Attributes.
 *
 * This parser needs to be tested on each new GFF file, the GFF format is so
	flexible as to cause all kinds of problems.
 *
 * Works with flybase, and the drosophila release 3.1 gff annotation.
 */
public class GffParser {
    
    private String fileName;
    private String[] allFeatures;  //contains all types of features found in file
    private String[] allSources;   //ditto
    private String[] allAttributes; //ditto
    private String comments;
    private GffFeature[] gffLineObjs;
    
    //private GffFeature[] features;  //contains an array of GffFeatures
    
    public GffParser(File file, int startNt, int stopNt) {
    //read in GFF file one line at a time
        int counter =0;
        try{
            BufferedReader in = new BufferedReader(new FileReader (file));  
            String line;
            StringBuffer sb = new StringBuffer();
            ArrayList al = new ArrayList();
                
            while ((line = in.readLine()) != null){
                line = line.trim(); //kill whitespace on ends
                if (line.length() == 0) continue; //skip blank lines
                if (line.startsWith("#")) {  //save comments
                    sb.append(line+"\n");
                }   
                else{   
                    GffFeature gf = new GffFeature(line);
                    if (gf.isGffGood()==false) break;
                    //check if the feature is inbetween the start and stop
                    if(gf.getStart()>=startNt && gf.getEnd()<=stopNt ) al.add(gf); 
					//System.out.println(gf);
                }
                counter++;
                if (counter > 20000000){
					JOptionPane.showMessageDialog(
						null,"Sorry, your gff file is too large and has been truncated!\n",null,JOptionPane.WARNING_MESSAGE);
                	 break; } //limit GFF upload
            }
            //save comments and convert arraylist to array
            comments = sb.toString();
            gffLineObjs = new GffFeature[al.size()];
            al.toArray(gffLineObjs);
        }
        catch (IOException e){
            e.printStackTrace();
        }
        /*
        System.out.println("Total number of lines processed: "+counter);
        int orig = gffLineObjs.length;
        System.out.println("Number of GFFline objects created: "+orig);
        System.out.println("Summary listing of all the features, sources, and attributes found in the GFF file.");
        quantifyItems();
        */
    }
    public GffFeature[] getGffFeatures(){return gffLineObjs;}
    
    public void quantifyItems(){
        LinkedHashSet features = new LinkedHashSet(100);
        LinkedHashSet sources = new LinkedHashSet(10);
        LinkedHashMap attributes = new LinkedHashMap(100);
        
        int len = gffLineObjs.length;
        for (int i=0; i<len; i++){
            features.add(gffLineObjs[i].getFeature());
            sources.add(gffLineObjs[i].getSource());
            attributes.putAll(gffLineObjs[i].getAttsHash());
        }
        System.out.println("Features: "+features);
        System.out.println("Sources: "+sources);
        System.out.println("Attributes: "+attributes);
        
    }

    //public static void main(String[] args) {
    //    GffParser gf = new GffParser("/Users/nix/Desktop/GIClass/Generic-Genome-Browser-1.54/sample_data/yeast_data.gff", 0, 1000000000);
    //} 
    
}
