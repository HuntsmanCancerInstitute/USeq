package util.bio.parsers;

import java.io.*;
import java.util.*;

import util.gen.IO;

/**
 * Extracts sequences and headers out of multi fasta files.
Will read in combine fasta files, skips blank lines, does not filter actual seq
parses off the > character. Seqs and names are index matched, ie seqs[3] and names[3] refer 
to the same thing.*/

public class MultiFastaParser {

	private String[] seqs;
	private ArrayList seqsAL = new ArrayList();
	private String[] names;
	private ArrayList namesAL = new ArrayList();
	private Fasta[] fastas;
	private HashMap namesSeqs;
	private HashMap seqsNames;
	private boolean fastaFound = false;
	private int minLength = 0;
	private boolean justPrint = false;

	//Constructors
	public MultiFastaParser(File file){ 
		seqsAL.clear();
		namesAL.clear();
		parseIt(file);
	}
	public MultiFastaParser(File file, int minLength, boolean justPrint){ 
		seqsAL.clear();
		namesAL.clear();
		this.minLength = minLength;
		this.justPrint = justPrint;
		parseIt(file);
	}
	public MultiFastaParser(String file) {
		seqsAL.clear();
		namesAL.clear();
		parseIt(new File(file));
	}
	public MultiFastaParser(){}

	public void resetFields(){
		seqs = null;
		names = null;
		fastaFound = false;
		seqsAL.clear();
		namesAL.clear();
	}

	public MultiFastaParser(File[] fastaFiles){
		for (int i=0; i< fastaFiles.length; i++){
			parseIt(fastaFiles[i]);
		}
	}

	//Methods
	
	/**Every call adds to the list of sequences and names.*/
	public void parseIt (File file){   
		StringBuffer sb = new StringBuffer();
		String seqName = null;
		String line;
		boolean inside = false;
		int counter = 0;


		try {
			BufferedReader in = IO.fetchBufferedReader(file);
			while ((line = in.readLine()) !=null) {
				line = line.trim();                     //kill whitespace and test if exists
				if (line.length() == 0) continue;       //skip blank lines
				if (line.startsWith(">")){          //if a comment line
					if (inside) {
						if (minLength == 0 || sb.length() >= minLength){
							if (justPrint){
								System.out.println(">"+seqName+"\n"+sb);
							}
							else {
								if (sb.length()>0) seqsAL.add(new String(sb));
								if (seqName != null) namesAL.add(seqName);      
							}
						}
						sb = new StringBuffer();
					}
					//set new seqName;
					seqName = line.substring(1).trim();
					counter++;
				}   
				else {
					inside = true;
					sb.append(line);
				}   
			}
			in.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}

		//add last?
		if (minLength == 0 || sb.length() >= minLength){
			if (justPrint){
				System.out.println(">"+seqName+"\n"+sb);
			}
			else {
				if (sb.length()>0) seqsAL.add(new String(sb));
				if (seqName != null) namesAL.add(seqName);
			}  		
		}

		//check
		if (namesAL.size()>0  && seqsAL.size() == namesAL.size()) fastaFound=true;
		sb = null;
		//null fields due to new addition
		namesSeqs = null;
		seqsNames = null;
		names = null;
		seqs = null;
		fastas = null;

	}

	public void printFasta() {
		int num = seqsAL.size();
		for (int i=0; i<num; i++){
			System.out.println(">"+namesAL.get(i));
			System.out.println(seqsAL.get(i));
		}
	}

	//getter methods
	/**Returns a hashmap of the fasta text: sequence*/
	public HashMap getNamesSeqs(){
		if (namesSeqs == null){
			int num = seqsAL.size();
			namesSeqs = new HashMap (num);
			for (int i=0; i< num; i++) namesSeqs.put(namesAL.get(i), seqsAL.get(i));
		}
		return namesSeqs;
	}
	/**Returns a hashmap of the fasta sequence:text*/
	public HashMap getSeqsNames(){
		if (seqsNames == null){
			int num = seqsAL.size();
			seqsNames = new HashMap (seqs.length);
			for (int i=0; i< num; i++) seqsNames.put(seqsAL.get(i),namesAL.get(i));
		}
		return seqsNames;
	}

	public String[] getSeqs(){
		if (seqs == null){
			seqs = new String[seqsAL.size()];
			seqsAL.toArray(seqs);
		}
		return seqs;
	}
	public String[] getNames(){
		if (names == null){
			names = new String[namesAL.size()];
			namesAL.toArray(names);
		}
		return names;
	}
	public Fasta[] getFastas(){
		if (fastas == null){
			int num = namesAL.size();
			fastas = new Fasta[num];
			for (int i=0; i< num; i++){
				fastas[i] = new Fasta ((String)namesAL.get(i), (String)seqsAL.get(i));
			}
		}
		return fastas;
	}
	public int getNumReads(){return seqsAL.size();}

	public static void main (String[] args){
        MultiFastaParser m = new MultiFastaParser(args[0]);
        System.out.println("Found seq? "+m.isFastaFound());
        System.out.println("Length "+m.getSeqs()[0].length());
        System.out.println("Names[0] "+m.getNames()[0]);
        System.out.println("Seqs[0] "+m.getSeqs()[0]);
        File parent = new File (args[0]).getParentFile();
        //print em?
        for (int i=0; i<m.getSeqs().length; i++){
        	String name = m.getNames()[i].replace(".1", "");
        	File f = new File (parent, name+".fasta");
        	IO.writeString(">"+name+"\n"+m.getSeqs()[i], f);
        }
    }
	public boolean isFastaFound() {
		return fastaFound;
	}
	public int getMinLength() {
		return minLength;
	}
	public void setMinLength(int minLength) {
		this.minLength = minLength;
	}
	public boolean isJustPrint() {
		return justPrint;
	}
	public void setJustPrint(boolean justPrint) {
		this.justPrint = justPrint;
	}
	public ArrayList getNamesAL() {
		return namesAL;
	}
	public void setNamesAL(ArrayList namesAL) {
		this.namesAL = namesAL;
	}
	public ArrayList getSeqsAL() {
		return seqsAL;
	}
	public void setSeqsAL(ArrayList seqsAL) {
		this.seqsAL = seqsAL;
	}
	public void setFastaFound(boolean fastaFound) {
		this.fastaFound = fastaFound;
	}

}
