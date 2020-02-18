package edu.utah.seq.vcf.json;

import java.io.*;
import java.util.regex.*;
import edu.utah.seq.vcf.VCFParser;
import edu.utah.seq.vcf.xml.SimpleVcf;
import util.gen.*;
import java.util.*;

/**
 * Takes a patient vcf file parsed from a Tempus json report and compares it to a vcf generated by reprocessing the raw data.
 * Writes out a final arbitrated vcf containing all the Tempus calls plus non duplicate recall variants with no FILTER flags.
 * 
 * @author david.nix@hci.utah.edu 
 **/
public class TempusVcfComparator {

	//user defined fields
	private File tempusVcf = null;
	private File recallVcf = null;
	private File mergedVcf = null;
	private SimpleVcf[] fVcfs;
	private SimpleVcf[] rVcfs;	
	private int bpPaddingForOverlap = 2;
	private boolean appendChr = false;
	private boolean excludeInherited = false;
	private boolean excludeSomatic = false;

	//counters
	private int numberShortTempus = 0;
	private int numberOtherTempus = 0;
	private int numberRecall = 0;
	private int numberExactMatches = 0;
	private int numberTempusWithOnlyOverlap = 0;
	private int numberModifiedTempusCalls = 0;
	private int numberTempusWithNoMatch = 0;
	private int numberPassingRecallWithNoMatch = 0;
	private int numberInherited = 0;
	private int numberSomatic = 0;
	
	private ArrayList<SimpleVcf> vcfToPrint = new ArrayList<SimpleVcf>();
	private ArrayList<String> headerLines = new ArrayList<String>();
	
	
	

	//constructors
	public TempusVcfComparator(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			//load vcf files
			fVcfs = load(tempusVcf, true);
			rVcfs = load(recallVcf, true);
			numberRecall = rVcfs.length;

			compareVcfs();

			processTempusVcfs();
			
			processRecallVcfs();
			
			printVcfs();
			
			printStats();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem running TempusJson2Vcf app!");
		}
	}

	private void printStats() {
		System.out.println("\nComparator stats:");
		System.out.println( numberRecall +"\t# Recall variants");
		System.out.println( numberShortTempus +"\t# Short Tempus variants");
		System.out.println( numberOtherTempus +"\t# Other Tempus variants");
		System.out.println( numberInherited +"\t# Germline Tempus variants, skippped? "+excludeInherited);
		System.out.println( numberSomatic +"\t# Somatic Tempus variants, skippped? "+excludeSomatic);
		System.out.println( numberExactMatches +"\t# Short with an exact match");
		System.out.println( numberTempusWithOnlyOverlap +"\t# Short with overlap recal variants");
		System.out.println( numberModifiedTempusCalls +"\t# Short recommended for modification");
		System.out.println( numberTempusWithNoMatch +"\t# Short with no match"); 
		System.out.println( numberPassingRecallWithNoMatch +"\t# Passing recall variants with no Short match");
	}

	private void printVcfs() {
		//sort vcf
		SimpleVcf[] vcf = new SimpleVcf[vcfToPrint.size()];
		vcfToPrint.toArray(vcf);
		Arrays.sort(vcf);
		
		try {
			Gzipper out = new Gzipper(mergedVcf);
			//fetch merged header
			String[] header = mergeHeaders(headerLines);
			
			for (String h: header) out.println(h);
			for (SimpleVcf v: vcf) {
				if (v.getFilter().toLowerCase().contains("fail") == false) out.println(v.getVcfLine());
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: problem writing out the merged vcf file. "+mergedVcf);
		}
	}

	private void processRecallVcfs() {
		//for each Recall vcf, call this after processing the Tempus vcfs, skip those with a fail filter field
		for (SimpleVcf r:rVcfs){
			//print it?
			if (r.isPrint() && r.getFilter().toLowerCase().contains("fail") == false) {
				vcfToPrint.add(r);
				if (r.getMatch() == null) numberPassingRecallWithNoMatch++;
			}
		}
	}
	
	/**Merges header lines eliminating duplicates.  Does a bad ID name collision checking, silently keeps first one. 
	 * Returns null if CHROM lines differ. */
	public static String[] mergeHeaders(ArrayList<String> header) {
		
		LinkedHashSet<String> other = new LinkedHashSet<String>();
		LinkedHashSet<String> contig = new LinkedHashSet<String>();
		LinkedHashSet<String> info = new LinkedHashSet<String>();
		LinkedHashSet<String> filter = new LinkedHashSet<String>();
		LinkedHashSet<String> format = new LinkedHashSet<String>();
		TreeSet<String> source = new TreeSet<String>();
		String chromLine = null;

		for (String h: header){
			h=h.trim();
			if (h.startsWith("##contig")){
				if (contig.contains(h) == false) contig.add(h);
			}
			else if (h.startsWith("##INFO")){
				if (info.contains(h) == false) info.add(h);
			}
			else if (h.startsWith("##FILTER")){
				if (filter.contains(h) == false) filter.add(h);
			}
			else if (h.startsWith("##FORMAT")){
				if (format.contains(h) == false) format.add(h);
			}
			else if (h.startsWith("##source=")){
				source.add(h);
			}
			else if (h.startsWith("#CHROM")) {
				if (chromLine == null) chromLine = h;
				//skip this check, else if (chromLine.equals(h) == false) Misc.printErrAndExit("\nERROR: chrom lines differ!\n");;
			}
			else if (other.contains(h) == false) {
				other.add(h);
			}
		}


		//add in filter lines
		filter.add(SimpleVcf.ncFilter);

		//remove ID dups from contig, filter, format, info
		ArrayList<String> contigAL = VCFParser.mergeHeaderIds(contig);
		ArrayList<String> filterAL = VCFParser.mergeHeaderIds(filter);
		ArrayList<String> formatAL = VCFParser.mergeHeaderIds(format);
		ArrayList<String> infoAL = VCFParser.mergeHeaderIds(info);

		ArrayList<String> lines = new ArrayList<String>();
		for (String s : other) lines.add(s);
		for (String s : source) lines.add(s);
		for (String s : contigAL) lines.add(s);
		for (String s : filterAL) lines.add(s);
		for (String s : infoAL) lines.add(s);
		for (String s : formatAL) lines.add(s);
		if (chromLine != null) lines.add(chromLine);

		return Misc.stringArrayListToStringArray(lines);
	}

	private void processTempusVcfs() {
		//for each Tempus record
		for (SimpleVcf f: fVcfs){

			//not a short? just save it
			if (f.isShortVariant() == false) {
				vcfToPrint.add(f);
				numberOtherTempus++;
			}
			else {
				numberShortTempus++;
				//exact match? 
				if (f.getMatch() != null) {
					numberExactMatches++;
					//exact match then add tempus info to recall
					SimpleVcf vcf = f.getMatch();
					vcf.appendID(f);
					vcf.appendINFO(f);
					f.setPrint(false);
				}
				else {
					//So no exact match any overlap?
					if (f.getOverlap().size()!=0) numberTempusWithOnlyOverlap++;

					//No exact or overlap
					else {
						System.err.println("WARNING: No match to this Tempus variant.");
						System.err.println("F:\t"+f.getVcfLine());
						numberTempusWithNoMatch++;
					}
					//always print it
					f.appendFilter("NC");
					vcfToPrint.add(f);
				}

			}
		}
	}

	private void compareVcfs() {
		//slow comparator, could do many things to speed up...
		for (int i=0; i< fVcfs.length; i++){
			SimpleVcf f = fVcfs[i];
			//is it a short variant?
			if (f.isShortVariant() == false) continue;

			for (int j=0; j< rVcfs.length; j++){
				SimpleVcf r = rVcfs[j];
				if (f.compareToExact(r)){
					if (f.getMatch() !=null || r.getMatch() != null) Misc.printErrAndExit("\nERROR: more than one exact match found for \n"+f+"\n"+r);
					f.setMatch(r);
					r.setMatch(f);
				}
				else if (f.compareToOverlap(r)){
					f.getOverlap().add(r);
					r.getOverlap().add(f);
				}
			}
		}

	}

	private SimpleVcf[] load(File vcf, boolean excludeContig) {
		String[] lines = IO.loadFileIntoStringArray(vcf);
		ArrayList<SimpleVcf> al = new ArrayList<SimpleVcf>();
		for (String v: lines){
			if (v.startsWith("#") == false) {
				if (v.contains("inherited")) {
					numberInherited++;
					if (excludeInherited) continue;
				}
				if (v.contains("somatic")) {
					numberSomatic++;
					if (excludeSomatic) continue;
				}
				if (appendChr && v.startsWith("chr") == false) v = "chr"+v;
				al.add(new SimpleVcf(v, bpPaddingForOverlap));
			}
			else {
				if (excludeContig){
					if (v.startsWith("##contig") == false) headerLines.add(v);
				}
				else headerLines.add(v);
			}
		}
		SimpleVcf[] svs = new SimpleVcf[al.size()];
		al.toArray(svs);
		return svs;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new TempusVcfComparator(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		String source = useqVersion+" Args: "+ Misc.stringArrayToString(args, " ");
		System.out.println("\n"+ source +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 't': tempusVcf = new File(args[++i]); break;
					case 'r': recallVcf = new File(args[++i]); break;
					case 'm': mergedVcf = new File(args[++i]); break;
					case 'c': appendChr = true; break;
					case 'g': excludeInherited = true; break;
					case 's': excludeSomatic = true; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//check files
		if (tempusVcf == null || recallVcf == null) Misc.printErrAndExit("\nError: cannot find both of your vcf files to compare?!\n");
		if (mergedVcf == null) Misc.printErrAndExit("\nError: please provide a named file for writing the merged vcf!\n");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Tempus Vcf Comparator: March 2019                      **\n" +
				"**************************************************************************************\n" +
				"TVC compares a Tempus vcf generated with the TempusJson2Vcf to a recalled vcf.\n"+
				"Exact recall vars are so noted and removed. Tempus vcf with no exact but one\n"+
				"overlapping record can be merged with -k. Be sure to vt normalize each before running.\n"+
				"Recall variants failing FILTER are not saved.\n"+

				"\nOptions:\n"+
				"-t Path to a TempusOne vcf file, see the TempusJson2Vcf app.\n"+
				"-r Path to a recalled snv/indel vcf file.\n"+
				"-m Path to named vcf file for saving the results.\n"+
				"-c Append chr if absent in chromosome name.\n"+
				"-g Exclude 'inherited' germline Tempus records from the comparison and merged output.\n"+
				"-s Exclude 'somatic' tumor Tempus records from the comparison and merged output.\n"+

				"\nExample: java -Xmx2G -jar pathToUSeq/Apps/TempusVcfComparator -f TL-18-03CFD6.vcf\n" +
				"     -r /F1/TL-18-03CFD6_recall.vcf.gz -g -c -m /F1/TL-18-03CFD6_merged.vcf.gz -k \n\n" +

				"**************************************************************************************\n");
	}



}