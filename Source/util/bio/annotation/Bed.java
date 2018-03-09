package util.bio.annotation;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import edu.utah.seq.its.Interval1D;
import edu.utah.seq.its.IntervalST;
import edu.utah.seq.useq.data.Region;
import edu.utah.seq.useq.data.RegionScoreText;
import edu.utah.seq.vcf.VCFLookUp;
import edu.utah.seq.vcf.VCFParser;
import edu.utah.seq.vcf.VCFRecord;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class Bed extends Coordinate implements Serializable{

	//fields
	private String name;
	private double score; 
	private char strand;
	private static final long serialVersionUID = 1L;
	public static Pattern END = Pattern.compile(".*END=(\\d+).*");

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
	
	/**Parses a tab delimited txt file containing chr start stop. Skips # lines and spaces.
	 * @param subStart - bases to be subtracted from region starts
	 * @param subEnd - bases to be subtracted from region ends
	 * @param chrStartStopIndexes - three indexes that define the chr start stop columns, for bed use 0,1,2
	 * */
	public static Bed[] parseFilePutLineInNameNoScoreOrStrand(File bedFile, int subStart, int subEnd, int[] chrStartStopIndexes){
		Bed[] bed =null;
		String line = null;
		try{
			BufferedReader in = IO.fetchBufferedReader(bedFile);
			String[] tokens;
			ArrayList<Bed> al = new ArrayList<Bed>();
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() ==0 || line.startsWith("#")) continue;
				tokens = line.split("\\s+");
				if (tokens.length < 3) continue;				
				al.add(new Bed(tokens[chrStartStopIndexes[0]], Integer.parseInt(tokens[chrStartStopIndexes[1]])-subStart, Integer.parseInt(tokens[chrStartStopIndexes[2]])- subEnd, line, 0, '.'));
			}
			bed = new Bed[al.size()];
			al.toArray(bed);
		}catch (Exception e){
			e.printStackTrace();
			System.out.println("Bad line? -> "+line);
		}
		return bed;
	}

	
	/**Chunks the Bed[] by number of bp, then sorts each.*/
	public static ArrayList<Bed[]> splitByBp(Bed[] regions, int numBpsPerSplit) {
		ArrayList<Bed[]> split = new ArrayList<Bed[]>();
		int bpCount = 0;
		ArrayList<Bed> workingRegions = new ArrayList<Bed>();
		for (int i=0; i< regions.length; i++){
			int len = regions[i].getLength();
			if ((bpCount+len) > numBpsPerSplit){
				//close and save
				Bed[] set = new Bed[workingRegions.size()];
				workingRegions.toArray(set);
				split.add(set);
				workingRegions.clear();
				bpCount = 0;
			}
			//add region
			workingRegions.add(regions[i]);
			bpCount += len;
		}
		//save last?
		if (workingRegions.size() !=0){
			Bed[] set = new Bed[workingRegions.size()];
			workingRegions.toArray(set);
			Arrays.sort(set);
			split.add(set);
		}
		return split;
	}
	
	/**Chunks the Bed[] by number of regions, then sorts each.*/
	public static ArrayList<Bed[]> splitByNumber(Bed[] regions, int numRegionsPerChunk) {
		ArrayList<Bed[]> split = new ArrayList<Bed[]>();
		int count = 0;
		ArrayList<Bed> workingRegions = new ArrayList<Bed>();
		for (int i=0; i< regions.length; i++){
			if ((count+1) > numRegionsPerChunk){
				//close and save
				Bed[] set = new Bed[workingRegions.size()];
				workingRegions.toArray(set);
				Arrays.sort(set);
				split.add(set);
				workingRegions.clear();
				count = 0;
			}
			//add region
			workingRegions.add(regions[i]);
			count++;
		}
		//save last?
		if (workingRegions.size() !=0){
			Bed[] set = new Bed[workingRegions.size()];
			workingRegions.toArray(set);
			Arrays.sort(set);
			split.add(set);
		}
		return split;
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
	
	/**Splits the Bed[] by chrom, does not add strand onto the chrom name.
	 * Be sure to sort before calling! */
	public static HashMap<String, Bed[]> splitBedByChrom(Bed[] sortedBed) {
		HashMap<String, Bed[]> chromRegions = new HashMap<String, Bed[]>();
		ArrayList<Bed> recordsAL = new ArrayList<Bed>();

		//set first
		String oldChrom = sortedBed[0].getChromosome();
		recordsAL.add(sortedBed[0]);

		//for each region
		for (int i=1; i< sortedBed.length; i++){
			String testChrom = sortedBed[i].getChromosome();

			//is it the same chrom?
			if (oldChrom.equals(testChrom) == false){
				//close old
				Bed[] v = new Bed[recordsAL.size()];
				recordsAL.toArray(v);
				recordsAL.clear();
				chromRegions.put(oldChrom, v);
				oldChrom = testChrom;
				if (chromRegions.containsKey(oldChrom)) {
					System.err.println("ERROR: problem with spliting Bed by chrom, looks like the array wasn't sorted!\n");
					return null;
				}
			}
			//save info
			recordsAL.add(sortedBed[i]);
		}
		//set last
		Bed[] v = new Bed[recordsAL.size()];
		recordsAL.toArray(v);
		chromRegions.put(oldChrom, v);

		return chromRegions;
	}

	/**Creates an IntervalTree for rapid intersection.  The returning string is just the chr name.
	 * @throws Exception */
	public static HashMap<String,IntervalST<String>> parseRegions(File bedFile) throws Exception{
		HashMap<String,IntervalST<String>> chrAls = new HashMap<String,IntervalST<String>>();
		String line = null;
		BufferedReader in = IO.fetchBufferedReader(bedFile);
		String[] tokens;
		String currChr = "";
		IntervalST<String> al = null;
		
		while ((line = in.readLine()) !=null) {
			line = line.trim();
			if (line.length() ==0 || line.startsWith("#")) continue;
			tokens = Misc.TAB.split(line);
			if (tokens.length < 3) throw new Exception("Error: the following bed line does not contain requisit number of tab delimited columns (chr, start, stop...)\n\t-> "+line);

			//fetch IntervalST
			if (currChr != tokens[0]){
				currChr = tokens[0];
				if (chrAls.containsKey(currChr)) al = chrAls.get(currChr);
				else {
					al = new IntervalST<String>(); 
					chrAls.put(currChr, al);
				}
			}
			//add entry, end is included in search so subtract 1
			int start = Integer.parseInt(tokens[1]);
			int stop = Integer.parseInt(tokens[2])-1;
			al.put(new Interval1D(start, stop), currChr);
		}
		return chrAls;
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
	
	public static boolean scoresSet(File bedFile){
		BufferedReader in = null;
		try {
			in = IO.fetchBufferedReader(bedFile);
			String[] tokens;
			String line;
			while ((line = in.readLine()) !=null) {
				//chrom, start, stop, text, score, strand
				//0       1       2     3     4      5
				line = line.trim();
				if (line.length() ==0 || line.startsWith("#")) continue;
				tokens = Misc.TAB.split(line);
				if (tokens.length > 4) return true;
				return false;
			}
			
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (in != null)
				try {
					in.close();
				} catch (IOException e) {}
		}
		return false;
	}


	/**Split a bed file by chromosome and strand into a HashMap of chromosomeStrand (e.g. chr3+, chr3-, chr3. or chr3 if ignoreStrand==true) : sorted NamedScoredCoordinate[].
	 * If appendChr, chr is appended to chrom names that lack it.*/
	public static HashMap<String,RegionScoreText[]> parseBedFile(File bedFile, boolean ignoreStrand, boolean appendChr){
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
				if (appendChr && chromStrand.startsWith("chr") == false) chromStrand = "chr"+chromStrand;
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
				String name = "";
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
			Misc.printErrAndExit("\nMalformed bed line? -> "+line+"\n Aborting!");
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
	
	/**Splits a vcf file by chromosome into a HashMap of RegionScoreText[] for variants that contain an END=xxx String.
	 * Places the record number in the bed score column and the vcf ID in the bed name column.*/
	public static HashMap<String,RegionScoreText[]> parseVcfFileForENDVars(File vcfFile){
		HashMap<String,ArrayList<RegionScoreText>> chrAls = new HashMap<String,ArrayList<RegionScoreText>>();
		String line = null;
		try{
			BufferedReader in = IO.fetchBufferedReader(vcfFile);
			ArrayList<RegionScoreText> al = new ArrayList<RegionScoreText>(); 
			String currChrom = "";
			Matcher mat;
			int recordNumber = -1;
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() ==0 || line.startsWith("#")) continue;
				recordNumber++;
				//look for end
				mat = END.matcher(line);
				if (mat.matches() == false) continue;
				//pull stop, set start, and fetch chr
				int stop = Integer.parseInt(mat.group(1));
				int start = stop -1;
				String[] fields = Misc.TAB.split(line);
				String chr = fields[0];
				//save it
				if (currChrom.equals(chr) == false) {
					al = chrAls.get(chr);
					if (al == null) al = new ArrayList<RegionScoreText>();
					chrAls.put(chr, al);
				}
				al.add(new RegionScoreText(start, stop, recordNumber, fields[2]));
			}
		}catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nMalformed vcf line? -> "+line+"\n Aborting!");
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

	
	/**Splits a vcf file by chromosome into a HashMap of RegionScoreText[]. Use the padding to expand the 
	 * size of the variant. SNV, INS, or DEL is assigned to the name field, QUAL to the score field.*/
	public static HashMap<String,RegionScoreText[]> parseVcfFile(File vcfFile, int padding, boolean appendChrFixMT){
		HashMap<String,ArrayList<RegionScoreText>> chrAls = new HashMap<String,ArrayList<RegionScoreText>>();
		String line = null;
		try{
			VCFParser vp = new VCFParser();
			BufferedReader in = IO.fetchBufferedReader(vcfFile);
			ArrayList<RegionScoreText> al = new ArrayList<RegionScoreText>(); 
			String currChrom = "";
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() ==0 || line.startsWith("#")) continue;
				VCFRecord vcf = new VCFRecord(line, vp, false, false);
				if (appendChrFixMT) {
					vcf.appendChr();
					vcf.correctChrMTs();
				}
				
				//fetch ArrayList
				if (currChrom != vcf.getChromosome()){
					currChrom = vcf.getChromosome();
					if (chrAls.containsKey(currChrom)) al = chrAls.get(currChrom);
					else {
						al = new ArrayList<RegionScoreText>(); 
						chrAls.put(currChrom, al);
					}
				}
				
				//positions
				int start = vcf.getPosition() - padding;
				if (start<0) start = 0;
				int stop;
				String type;
				if (vcf.isDeletion()){
					stop = vcf.getPosition() + vcf.getReference().length() + 1 + padding;
					type = "DEL";
				}
				else {
					stop = vcf.getPosition() + vcf.getReference().length() + padding;
					if (vcf.isSNP()) type = "SNV";
					else type = "INS";
				}
				RegionScoreText n = new RegionScoreText(start, stop, vcf.getQuality(), type+"_"+ vcf.getReference() +"_"+Misc.stringArrayToString(vcf.getAlternate(), ","));
				al.add(n);
			}
		}catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nMalformed vcf line? -> "+line+"\n Aborting!");
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
	public String toStringNoStrandNoName(){
		return chromosome+"\t"+start+"\t"+stop+"\t"+score;
	}
	public String toString(){
		return chromosome+"\t"+start+"\t"+stop+"\t"+name+"\t"+score+"\t"+strand;
	}
	public String toStringScoreName() {
		return chromosome+"\t"+start+"\t"+stop+"\t"+score+"\t"+name;
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
	/*Quickie method to stat a bed file.*/
	public static void main (String[] args){
		if (args.length ==0 ) Misc.printErrAndExit("\nEnter a vcf file and bp to pad.");
		int padding = Integer.parseInt(args[1]);
		File vcf = new File (args[0]);
		File bed = new File (vcf.getParentFile(), Misc.removeExtension(vcf.getName())+"Pad"+ padding +"bp.bed.gz");
		try {
			Gzipper out = new Gzipper(bed);
			HashMap<String,RegionScoreText[]> chrReg = Bed.parseVcfFile(vcf, padding, false);
			for (String chr: chrReg.keySet()){
				RegionScoreText[] regions = chrReg.get(chr);
				for (RegionScoreText r: regions){
					out.println(r.getBedLine(chr));
				}
			}
			out.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
	}

}
