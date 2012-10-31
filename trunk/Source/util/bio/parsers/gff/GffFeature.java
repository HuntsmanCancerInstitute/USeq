package util.bio.parsers.gff;

import java.util.*;
import javax.swing.*;

/**
A container for each of the values in a GFF 2 line.
<seqname>tab<source>tab<feature>tab<start>tab<stop>tab<score>tab
<strand>tab<frame>tab[space delimited key=value; key=value; attributes]
Each param is split into a String[] on \t
Attributes are further processed splitting on ;\s into another String[] each
of theses Strings are then split into key value pairs on an = or \s and stuffed
into a LinkedHashMap
Note Stop = End +1, needed to convert to real coordinates!
*/
public class GffFeature {
    
    private String gffLine;     //incoming line from GFF file
    private String seqName;     //all the following are the GFF(2) labels
    private String source;      //  to each tab delimited catagory
    private String feature;
    private int start;          
    private int end;			
    private double score;          
    private String strand;
    private String frame;
    private String attributes;
    private LinkedHashMap attributesHash;
    private boolean gffGood = true;
    
    public GffFeature(String unParsedGffLine){
        //make field assignments
        gffLine = unParsedGffLine;
        String[] items = gffLine.split("\\t");
        //check length
        if (items.length!=9){
			JOptionPane.showMessageDialog(
				null,
				"Sorry, Your GFF file appears to be incomplete.  Too few arguments coming off this gff line:\n   "+
				gffLine+"\n Check the file, your annotation is incomplete.",null,JOptionPane.WARNING_MESSAGE);
        	gffGood = false;
        }
        else {
        	seqName = items[0];
        	source = items[1];
        	feature = items[2];
        	strand = items[6];
        	if (items[5].startsWith(".")) score =0;
        	else score = GffFeature.parseStringToDouble(items[5]);
        	frame = items[7];
        	attributes = items[8];
        	start = (int)GffFeature.parseStringToDouble(items[3]); //start is always less than stop, independent of strand/ orientation
        	end = (int)GffFeature.parseStringToDouble(items[4]);
       
        	//create HashMap of attributes
        	String[] atts = attributes.split(";\\s*"); //split on ; and zero or more spaces
        	attributesHash = new LinkedHashMap();
        	for (int i=0; i<atts.length; i++){
            	String[] kv = atts[i].split("=|\\s"); //splitting on = or a space
            	if (kv.length>1) attributesHash.put(kv[0], kv[1]);
            	else if (kv.length==1) attributesHash.put(kv[0], "");
        	}
        }
    }
    public String toString(){
        StringBuffer ggf = new StringBuffer(
          "SeqName   : "+seqName+
        "\nSource    : "+source+
        "\nFeature   : "+feature+
        "\nStart     : "+start+
        "\nEnd       : "+end+
        "\nScore     : "+score+
        "\nStrand    : "+strand+
        "\nFrame     : "+frame+
        "\nAttributes: "+attributes+"\n");
        
        //pull hash
       Set enteries = attributesHash.entrySet();
       Iterator iter = enteries.iterator();
        while (iter.hasNext()){
            Map.Entry entry = (Map.Entry)iter.next();
            ggf.append("\tkey: "+entry.getKey()+"\tvalue: "+entry.getValue()+"\n");
        }
        return ggf.toString();
    }

    public static double parseStringToDouble(String text){
        try { return Double.parseDouble(text);}
        catch(NumberFormatException n) {return 0;
        }
    }
    
    public String getGffLine(){return gffLine;}
    public String getSeqName(){return seqName;}
    public String getSource(){return source;}
    public String getFeature(){return feature;}
    public int getStart(){return start;}
    public int getEnd(){return end;}
    public double getScore(){return score;}
    public String getStrand(){return strand;}
    public LinkedHashMap getAttsHash(){return attributesHash;}
    public String getAttsString(){return attributes;}
	public boolean isGffGood() {
		return gffGood;
	}

}
