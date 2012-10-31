package util.bio.parsers;

import java.io.*;
import java.util.*;

/**This class will read in combine quality files associated with raw sequencing results, skips blank lines, parses off the > character.
 * Quals and names are index matched, ie qual[3] and names[3] refer to the same thing.
 */

public class QualityFileParser {
    
    private String[] quals; //space delimited
    private String[] names;
    private String fileName;
    private HashMap qualNames;
    
    //Constructor
    public QualityFileParser(String file) {
        //Read in file line by line, skipping blank lines, only saving > and number lines
        fileName = file;
        StringBuffer sb = new StringBuffer();
        ArrayList namesAL = new ArrayList();
        ArrayList qualsAL = new ArrayList();
        String line;
        boolean inside = false;
        int counter = 0;
        
        try {
            BufferedReader in = new BufferedReader(new FileReader(file));
            while ((line = in.readLine()) !=null) {
                line = line.trim();                     //kill whitespace and test if exists
                if (line.length() == 0) continue;       //skip blank lines
                
                if (line.startsWith(">")){              //if a comment line
                    if (inside) {
                        qualsAL.add((new String(sb)).trim());  //need to clip off extra space at stop of text
                        sb = new StringBuffer();
                    }
                    namesAL.add(line.substring(1).trim());      //add comment to names[]
                    counter++;
                }
                else {
                    inside = true;
                    sb.append(line+" ");
                }
            }
            in.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        
        qualsAL.add((new String(sb)).trim());           //add last seq
        
        //convert ALs to arrays
        names = new String[namesAL.size()];
        namesAL.toArray(names);
        
        quals = new String[qualsAL.size()];
        qualsAL.toArray(quals);   
    }
    
    //primary methods
    public void printMFP() {
        for (int i=0; i<names.length; i++){
            System.out.println("Name: "+names[i]);
            System.out.println("Quals: "+quals[i]);
            System.out.println("");
        }
    }
    public void printReport() {
        System.out.println("***************************************");
        System.out.println("File Name              : "+ fileName);
        System.out.println("Number Quality Reports : "+quals.length);
        System.out.println("");
    }
    
    /**Returns a hashmap of the fasta text: sequence*/
    public HashMap getQualNames(){
    	if (qualNames == null){
    		qualNames = new HashMap (quals.length);
    		for (int i=0; i< quals.length; i++) qualNames.put(names[i], quals[i]);
    	}
    	return qualNames;
    }
    
    //getter methods
    public String[] getQuals(){return quals;}
    public String[] getNames(){return names;}
    public String getFileName(){return fileName;}
    
    /*public static void main (String[] args){
        QualityFileParser m = new QualityFileParser("/Users/nix/FlyBioIn/SELEX/TestData/NOB-001.fasta.screen.qual");
        m.printMFP();
        m.printReport();
    }*/
    
}
