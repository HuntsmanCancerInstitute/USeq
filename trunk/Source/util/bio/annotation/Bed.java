package util.bio.annotation;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import edu.utah.seq.useq.data.Region;
import edu.utah.seq.useq.data.RegionScoreText;


import util.gen.IO;

public class Bed extends Coordinate implements Serializable{

	//fields
	private String name;
	private double score; 
	private char strand;

	//constructor
	public Bed (String chromosome, int start, int stop, String name, double score, char strand){
		super(chromosome, start, stop);
		this.name = name;
		this.score = score;
		this.strand = strand;
	}
	
	public Bed (String chromosome, char strand, RegionScoreText nsss){
		super(chromosome, nsss.getStart(), nsss.getStop());
		this.name = nsss.getText();
		this.score = nsss.getScore();
		this.strand = strand;
	}


	/**Parses a tab delimited bed file: chrom, start, stop, text, score, strand. Only the first three are required.
	 * @param subStart - bases to be subtracted from region starts
	 * @param subEnd - bases to be subtracted from region ends
	 * */
	public static Bed[] parseFile(File bedFile, int subStart, int subEnd){
		Bed[] bed =null;
		String line = null;
		try{
			BufferedReader in = IO.fetchBufferedReader(bedFile);
			String[] tokens;
			ArrayList<Bed> al = new ArrayList<Bed>();
			//chrom, start, stop, text, score, strand
			//0       1       2     3     4      5
			int counter = 0;
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() ==0 || line.startsWith("#")) continue;
				tokens = line.split("\\s+");
				if (tokens.length < 3) continue;
				String name = ""+counter++;
				double score = 0;
				char strand = '.';
				if (tokens.length >= 6) {
					name = tokens[3];
					score = Double.parseDouble(tokens[4]);
					strand = tokens[5].charAt(0);
				}
				else if (tokens.length == 5){
					name = tokens[3];
					score = Double.parseDouble(tokens[4]);
				}
				else if (tokens.length == 4){
					name = tokens[3];
				}
				al.add(new Bed(tokens[0], Integer.parseInt(tokens[1])-subStart, Integer.parseInt(tokens[2])- subEnd, name, score, strand));
			}
			bed = new Bed[al.size()];
			al.toArray(bed);
		}catch (Exception e){
			e.printStackTrace();
			System.out.println("Bad line? -> "+line);
		}
		return bed;
	}

	/**Splits a bed file by chromosome and strand writing the data in binary form into the temp directory.
	 * Assumes a 6 column tab delimited bed file of chrom, start, stop, text, score, strand.
	 * Unsorted!*/
	public static HashMap<String,File> splitBedFileByChromStrandToTempDir(File bedFile, File tempDirectory){
		HashMap<String,File> splitData = new HashMap<String,File>();
		Pattern tab = Pattern.compile("\\t");
		String line = null;
		
		try{
			HashMap<String,DataOutputStream> chromOut = new HashMap<String,DataOutputStream>();
			BufferedReader in = IO.fetchBufferedReader(bedFile);
			String[] tokens;
			DataOutputStream dos = null;

			//chrom, start, stop, text, score, strand
			//0       1       2     3     4      5
			ArrayList<RegionScoreText> al = new ArrayList<RegionScoreText>(); 
			String currentChromStrand = "";
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() ==0 || line.startsWith("#")) continue;
				tokens = tab.split(line);
				if (tokens.length != 6) {
					System.err.println("Error: the following bed line does not contain 6 tab delimited columns (chr, start, stop, text, score, strand)\n\t-> "+line);
					continue;
				}
				//make chromosome strand text
				String chrStrand = tokens[0]+tokens[5];
				
				//get PrintWriter
				if (chrStrand.equals(currentChromStrand) == false){
					currentChromStrand = chrStrand;
					if(chromOut.containsKey(currentChromStrand)) dos = chromOut.get(chrStrand);
					else {
						//make and set file
						File f = new File(tempDirectory, chrStrand);
						splitData.put(chrStrand, f);
						dos = new DataOutputStream(new BufferedOutputStream( new FileOutputStream(f)));
						chromOut.put(chrStrand, dos);
					}
				}
				
				//write entry, start, stop, text, score
				int start = Integer.parseInt(tokens[1]);
				int stop = Integer.parseInt(tokens[2]);
				float score = Float.parseFloat(tokens[4]);
				dos.writeInt(start);
				dos.writeInt(stop);
				dos.writeFloat(score);
				dos.writeInt(tokens[3].length());
				dos.writeBytes(tokens[3]);
			}
			
			//close writers
			Iterator<DataOutputStream> it = chromOut.values().iterator();
			while (it.hasNext()) it.next().close();
			
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}
		
		return splitData;

	}
	
	/**Split a bed file by chromosome and strand into a HashMap of chromosomeStrand 
	 * (e.g. chr3+, chr3-, or chr3.; can force chr3. if ignoreStrand==true.)
	 * Will automatically add chr3. if missing strand info.*/
	public static HashMap<String,Region[]> parseRegions(File bedFile, boolean ignoreStrand){
		HashMap<String,ArrayList<Region>> chrAls = new HashMap<String,ArrayList<Region>>();
		Pattern tab = Pattern.compile("\\t");
		String line = null;
		try{
			BufferedReader in = IO.fetchBufferedReader(bedFile);
			String[] tokens;

			//chrom, start, stop, text, score, strand; might just be chr start stop
			//0       1       2     3     4      5
			ArrayList<Region> al = new ArrayList<Region>(); 
			String currentChromStrand = "";
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() ==0 || line.startsWith("#")) continue;
				tokens = tab.split(line);
				if (tokens.length < 3) {
					System.err.println("Error: skipping the following bed line, it does not contain requisit number of tab delimited columns (chr, start, stop...)\n\t-> "+line);
					continue;
				}
				//make chromosome strand text
				String chromStrand;
				if (ignoreStrand || tokens.length <6) chromStrand = tokens[0]+".";
				else chromStrand = tokens[0]+tokens[5];
				//fetch ArrayList
				if (currentChromStrand != chromStrand){
					currentChromStrand = chromStrand;
					if (chrAls.containsKey(currentChromStrand)) al = chrAls.get(currentChromStrand);
					else {
						al = new ArrayList<Region>(); 
						chrAls.put(currentChromStrand, al);
					}
				}
				//add entry, start, stop, text, score
				int start = Integer.parseInt(tokens[1]);
				int stop = Integer.parseInt(tokens[2]);
				Region n = new Region(start, stop);
				al.add(n);
			}
		}catch (Exception e){
			e.printStackTrace();
			System.out.println("Malformed bed line? -> "+line);
		}
		//sort and load hash
		HashMap<String,Region[]> chrSpec = new HashMap<String,Region[]>();
		Iterator<String> it = chrAls.keySet().iterator();
		ArrayList<Region> al = null;
		while (it.hasNext()){
			String cs = it.next();
			al = chrAls.get(cs);
			Region[] nsss = new Region[al.size()];
			al.toArray(nsss);
			Arrays.sort(nsss);
			chrSpec.put(cs, nsss);
		}
		return chrSpec;
	}


	/**Split a bed file by chromosome and strand into a HashMap of chromosomeStrand (e.g. chr3+, chr3-, chr3. or chr3 if ignoreStrand==true) : sorted NamedScoredCoordinate[].*/
	public static HashMap<String,RegionScoreText[]> parseBedFile(File bedFile, boolean ignoreStrand){
		HashMap<String,ArrayList<RegionScoreText>> chrAls = new HashMap<String,ArrayList<RegionScoreText>>();
		Pattern tab = Pattern.compile("\\t");
		String line = null;
		try{
			BufferedReader in = IO.fetchBufferedReader(bedFile);
			String[] tokens;

			//chrom, start, stop, text, score, strand
			//0       1       2     3     4      5
			ArrayList<RegionScoreText> al = new ArrayList<RegionScoreText>(); 
			String currentChromStrand = "";
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() ==0 || line.startsWith("#")) continue;
				tokens = tab.split(line);
				if (tokens.length < 3) {
					System.err.println("Error: skipping the following bed line, it does not contain requisit number of tab delimited columns (chr, start, stop...)\n\t-> "+line);
					continue;
				}
				//make chromosome strand text
				String chromStrand;
				if (ignoreStrand ) chromStrand = tokens[0];
				else if (tokens.length==3) chromStrand = tokens[0]+".";
				else chromStrand = tokens[0]+tokens[5];
				//fetch ArrayList
				if (currentChromStrand != chromStrand){
					currentChromStrand = chromStrand;
					if (chrAls.containsKey(currentChromStrand)) al = chrAls.get(currentChromStrand);
					else {
						al = new ArrayList<RegionScoreText>(); 
						chrAls.put(currentChromStrand, al);
					}
				}
				//add entry, start, stop, text, score
				int start = Integer.parseInt(tokens[1]);
				int stop = Integer.parseInt(tokens[2]);
				float score = 0;
				String name = ".";
				if (tokens.length == 4){
					name = tokens[3];
				}
				else if (tokens.length > 4) {
					score = Float.parseFloat(tokens[4]);
					name = tokens[3];
				}
				RegionScoreText n = new RegionScoreText(start, stop, score, name);
				al.add(n);
			}
		}catch (Exception e){
			e.printStackTrace();
			System.out.println("Malformed bed line? -> "+line);
		}
		//sort and load hash
		HashMap<String,RegionScoreText[]> chrSpec = new HashMap<String,RegionScoreText[]>();
		Iterator<String> it = chrAls.keySet().iterator();
		ArrayList<RegionScoreText> al = null;
		while (it.hasNext()){
			String cs = it.next();
			al = chrAls.get(cs);
			RegionScoreText[] nsss = new RegionScoreText[al.size()];
			al.toArray(nsss);
			Arrays.sort(nsss);
			chrSpec.put(cs, nsss);
		}
		return chrSpec;
	}
	
	public String toStringNoStrand(){
		return chromosome+"\t"+start+"\t"+stop+"\t"+name+"\t"+score;
	}
	public String toString(){
		return chromosome+"\t"+start+"\t"+stop+"\t"+name+"\t"+score+"\t"+strand;
	}
	public double getScore() {
		return score;
	}
	public void setScore(double score) {
		this.score = score;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public char getStrand() {
		return strand;
	}

	public void setStrand(char strand) {
		this.strand = strand;
	}
}
