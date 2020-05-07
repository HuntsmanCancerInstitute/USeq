package edu.utah.seq.vcf.sim;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.data.sam.SamLayoutLite;
import edu.utah.seq.query.QueryIndexFileLoader;
import edu.utah.seq.vcf.VCFLookUp;
import edu.utah.seq.vcf.VCFParser;
import edu.utah.seq.vcf.VCFRecord;
import util.bio.seq.Seq;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Given a vcf file of snv and indels, introduces the variants into the overlapping bam record and 
 * outputs fastq for realignment as well as two bam files, a variants hit stats file, and variants 
 * that were skipped (too close, not snv/indel).
 * @author Nix
 * */
public class BamBlaster {

	//fields
	private File bamFile;
	private File saveDirectory;
	private File vcfFile;
	private int maxSizeIndel = 50;
	private int minAlignmentDepth = 25;
	
	//internal fields
	private boolean warn = false;
	private File tempBam;
	private SamReaderFactory readerFactory;
	private SAMFileWriter tempBamWriter;
	private SAMFileWriter finalBamWriter;
	private SAMFileWriter unmodifiedFastqBamWriter;
	private SamReader bamReader;
	private HashMap<String, VCFLookUp> chromVariant;
	private HashSet<String> modifiedNames = new HashSet<String>();
	private PrintWriter varReportOut = null;
	private int minVarDist = 150;
	private static Pattern CIGAR_SUB = Pattern.compile("(\\d+)([MSDHN])");
	
	//per chrom fields
	private String workingChrom;
	private VCFRecord[] workingVcf;
	private Gzipper fastqFirstOut;
	private Gzipper fastqSecondOut;
	private Gzipper fastqUnpairedOut;
	private Gzipper excludedVCF;
	private Gzipper includedVCF;
	private int numberExcludedVcf = 0;
	private int numberProcessedVcf = 0;
	
	public BamBlaster (String[] args){
		long startTime = System.currentTimeMillis();

		processArgs(args);

		//for each chromosome of variants
		System.out.println("Processing variants....");
		for (String c: chromVariant.keySet()){
			workingChrom = c; 
			workingVcf = chromVariant.get(workingChrom).getVcfRecord();
			processVars();
		}
		System.out.println("\t"+numberProcessedVcf +"\t# Processed Vcf records");
		if (numberExcludedVcf !=0) System.out.println("\t"+numberExcludedVcf +"\t# Vcf records excluded (too close (<"+minVarDist+"bp), too big (>"+maxSizeIndel+"), lacking reference first base, or not a SNV/INDEL)");
		else excludedVCF.getGzipFile().deleteOnExit(); 
			
		//clear hashes and close temBamWriter
		chromVariant = null;
		modifiedNames = null;
		workingVcf = null;
		tempBamWriter.close();
		
		//find mates
		System.out.print("\nFinding mates, extracting fastq, and saving alignments (slow)...");
		matchMates();
		
		//cleanup
		try {
			finalBamWriter.close();
			bamReader.close();
			fastqFirstOut.close();
			fastqSecondOut.close();
			fastqUnpairedOut.close();
			excludedVCF.close();
			includedVCF.close();
			varReportOut.close();
			unmodifiedFastqBamWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR closing readers, writers, or gzippers?\n");
		}
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}
	
	private void matchMates(){
		int numPaired = 0;
		int numUnpaired = 0;
		
		//load hash of readName: SAMRecord from modified bam file
		HashMap<String, SAMRecord> readName2Pairs = new HashMap<String, SAMRecord>();
		SAMRecordIterator i = readerFactory.open(tempBam).iterator();
		SAMRecord dummy = new SAMRecord(null);
		
		while (i.hasNext()){
			SAMRecord sam = i.next();
			SAMRecord saved = readName2Pairs.get(sam.getReadName());
			//null so add first
			if (saved == null) readName2Pairs.put(sam.getReadName(), sam);
			//is it the dummy and thus already paired, throw error, this means that three alignments came through with the same name!
			else if (saved == dummy) Misc.printErrAndExit("\nERROR: more than two alignments were found for \n"+sam.getSAMString().trim()+ "\n"+saved.getSAMString().trim());
			//must be a pair so print
			else {
				//print pair and set dummy
				printPairedFastq(sam, saved);
				numPaired++;
				//set the dummy
				readName2Pairs.put(new String (sam.getReadName()), dummy);
			}
		}
		i.close();
		
		//OK now lets walk through entire original bam to parse out remaining mates and alignments that haven't been modified
		SAMRecordIterator it = bamReader.iterator();
		int numAlignments = 0;
		int numAlignmentsToUnmodifiedBam = 0;
		int counter = 0;
		
		while (it.hasNext()) {
			SAMRecord sam = it.next();
			numAlignments++;
			if (counter++ > 1000000){
				System.out.print(".");
				counter = 0;
			}
			//junk alignment? skip it
			if (sam.isSecondaryOrSupplementary() || sam.getReadUnmappedFlag()) continue;
			
			//is it present in the modified hash
			if (readName2Pairs.containsKey(sam.getReadName())){
				//write unmodified
				unmodifiedFastqBamWriter.addAlignment(sam);
				
				SAMRecord saved = readName2Pairs.get(sam.getReadName());
				//is this the dummy?
				if (saved != dummy){
					//is it the mate?
					if (areMates(sam, saved)) {
						//print pair
						printPairedFastq(sam, saved);
						numPaired++;
						//set the dummy
						readName2Pairs.put(new String(sam.getReadName()), dummy);
					}
					//is it the same? if not then there is a third alignment, throw error!
					else if (areSame(sam, saved) == false) Misc.printErrAndExit("\nERROR: more than two alignments upon spooling were found for \n"+sam.getSAMString().trim()+ "\n"+saved.getSAMString().trim());
				}
			}
			//OK it is not in the modified hash so just print it
			else {
				finalBamWriter.addAlignment(sam);
				numAlignmentsToUnmodifiedBam++;
			}
		}
		it.close();

		//now run through the modified hash and extract fastq as unpaired any non dummy record
		for (SAMRecord sam : readName2Pairs.values()){
			if (sam != dummy){
				String name = sam.getReadName();
				printFastq(sam, fastqUnpairedOut, name);
				numUnpaired++;
			}
		}
		
		System.out.println("\n\t"+ numAlignments+"\t# Alignments in ori bam file");
		System.out.println("\t"+ numAlignmentsToUnmodifiedBam+"\t# Unmodified alignments written to filtered.bam file");
		System.out.println("\t"+ numPaired+"\t# Paired fastq reads extracted for realignment");
		System.out.println("\t"+ numUnpaired+"\t# Unpaired fastq reads extracted for realignment");
		
		//delete any of the gzippers?
		if (numPaired == 0) {
			fastqFirstOut.getGzipFile().deleteOnExit();
			fastqSecondOut.getGzipFile().deleteOnExit();
		}
		if (numUnpaired == 0) fastqUnpairedOut.getGzipFile().deleteOnExit();
		
	}

	private void printFastq(SAMRecord sam, Gzipper out, String name) {
		try {
			String seq;
			String qual = null;
			//original qualities present?
			Object o = sam.getAttribute("OQ");
			if (o != null) qual = (String)o;
			else Misc.printErrAndExit("No OQ:Z tag with original quality scores in " + sam.getSAMString() + "\nAborting!");
			
			//negative strand? reverse em
			if (sam.getReadNegativeStrandFlag()) {
				seq = Seq.reverseComplementDNA(sam.getReadString());
				qual = Misc.reverse(qual);
			}
			else seq = sam.getReadString();
			out.print("@");
			out.println(name);
			out.println(seq);
			out.println("+");
			out.println(qual);
			
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem printing out fastq data to "+out.getGzipFile());
		}

	}

	private void printPairedFastq(SAMRecord sam, SAMRecord mate) {
		String name = sam.getReadName();
		
		//order by first second
		SAMRecord first = mate;
		SAMRecord second = sam;
		if (sam.getFirstOfPairFlag()) {
			first = sam;
			second = mate;
		}
		
		//who was modified
		String bb1 = null;
		Object ob = first.getAttribute("BB");
		if (ob != null) bb1= (String)ob;
		String bb2 = null;
		ob = second.getAttribute("BB");
		if (ob != null) bb2= (String)ob;
		
		//both modified?
		if (bb1 != null && bb2 != null){
			//same mod?
			if (bb1.equals(bb2)) name= name+":BB-"+bb1+"_12";
			else{
				//diff mods
				name = name+":BB-"+bb1+"_1"+":BB-"+bb2+"_2";
			}
		}
		//nope just one
		else if (bb1 != null) name = name+":BB-"+bb1+"_1";
		else name = name+":BB-"+bb2+"_2";
		
		if (sam.getFirstOfPairFlag()) {
			printFastq(sam, fastqFirstOut, name);
			printFastq(mate, fastqSecondOut, name);
		}
		else {
			printFastq(sam, fastqSecondOut, name);
			printFastq(mate, fastqFirstOut, name);
		}
	}


	private void processVars() {
		try {
			//for each variant
			int priorEnd = -1000;
			for (VCFRecord vcf: workingVcf){
				//check size, leading base, and for N's!
				if (vcf.getSizeIndel() > maxSizeIndel || vcf.checkLeadingRef() == false || vcf.getAlternate()[0].equals("N")){
					excludedVCF.println(vcf.getOriginalRecord());
					numberExcludedVcf++;
				}
				//check distance
				else if ((vcf.getPosition() - priorEnd) <= minVarDist){
					excludedVCF.println(vcf.getOriginalRecord());
					numberExcludedVcf++;
				}
				else {
					//fetch alignments
					ArrayList<SAMRecord> al = fetchAlignments(vcf);

					//enough coverage?
					if (al.size() < minAlignmentDepth ){
						excludedVCF.println(vcf.getOriginalRecord());
						numberExcludedVcf++;
					}
					else {
						numberProcessedVcf++;
						priorEnd= vcf.getMaxEndPosition();
						//attempt to modify
						modify(al, vcf);	
					}
					
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR printing excluded vcfs to file.\n");
		}
	}
	
	/**Assumes read names are the same.*/
	private boolean areMates(SAMRecord f, SAMRecord s){
		//check if opp pair
		if (f.getFirstOfPairFlag() == s.getFirstOfPairFlag()) return false;
		//check mate positions
		if (f.getMateAlignmentStart() != s.getAlignmentStart()) return false;
		if (s.getMateAlignmentStart() != f.getAlignmentStart()) return false;
		return true;
	}
	
	/**Assumes read names are the same.*/
	private boolean areSame(SAMRecord f, SAMRecord s){
		//check positions
		if (f.getAlignmentStart() != s.getAlignmentStart()) return false;
		//check if paired
		if (f.getReadPairedFlag() != s.getReadPairedFlag()) return false;
		//check if opp pair
		if (f.getFirstOfPairFlag() != s.getFirstOfPairFlag()) return false;
		//check cigar
		if (s.getCigarString().equals(f.getCigarString())  == false) return false;
		return true;
	}
	
	
	
	private void modify(ArrayList<SAMRecord> al, VCFRecord vcf) throws IOException {
		//substitute base qualities with original
		for (SAMRecord r: al){
			Object o = r.getAttribute("OQ");
			if (o == null) throw new IOException("Faild to find an OQ original quality tag in this sam alignment "+r.getSAMString());
			String qual = (String)o;
			r.setBaseQualityString(qual);
		}
		
		int numberOverlaps = al.size();
		int numberModified = 0;
		boolean included = false;
		if (vcf.isSNP()) {
			for (SAMRecord r: al) {
				//already modified?
				String qualName = getQualifiedName(r);
				if (modifiedNames.contains(qualName)) continue;
				if (addSnv(r, vcf)) {
					modifiedNames.add(qualName);
					numberModified++;
					tempBamWriter.addAlignment(r);
					included = true;
				}
			}
		}
		else if (vcf.isDeletion()){
			for (SAMRecord r: al) {
				//already modified?
				String qualName = getQualifiedName(r);
				if (modifiedNames.contains(qualName)) continue;
				if (addDeletion(r, vcf)) {
					modifiedNames.add(qualName);
					numberModified++;
					tempBamWriter.addAlignment(r);
					included = true;
					r.setOriginalBaseQualities(r.getBaseQualities());
				}
			}
		}
		else if (vcf.isInsertion()){
			for (SAMRecord r: al) {
				//already modified?
				String qualName = getQualifiedName(r);
				if (modifiedNames.contains(qualName)) continue;
				if (addInsertion(r, vcf)) {
					modifiedNames.add(qualName);
					numberModified++;
					tempBamWriter.addAlignment(r);
					included = true;
					r.setOriginalBaseQualities(r.getBaseQualities());
				}
			}
		}
		else {
			if (warn) System.err.println("WARNING: only SNVs and INDELs. Ignoring and writing to file:\n"+vcf);
			excludedVCF.println(vcf);
			numberExcludedVcf++;
		}

		if (included) includedVCF.println(vcf);
		//stats
		varReportOut.println(numberModified+"\t"+numberOverlaps+"\t"+vcf);
	}


	private boolean addInsertion(SAMRecord r, VCFRecord vcf) {
		int pos = vcf.getPosition();
		String refBases = vcf.getReference();
		String mutBases = vcf.getAlternate()[0];
		//make layout
		SamLayoutLite sl = new SamLayoutLite(r);
		//attempt mod
		if (sl.insertBases(pos, refBases, mutBases, warn)){
			String[] sq = sl.getSequenceAndQualtities();
			r.setReadString(sq[0]);
			r.setBaseQualityString(sq[1]);
			r.setAttribute("BB", "INS_"+pos+"_"+refBases+"_"+mutBases);
			return true;
		}
		if (warn) {
			System.err.println(vcf);
			System.err.println(r.getSAMString().trim());
		}
		return false;
	}


	private boolean addDeletion(SAMRecord r, VCFRecord vcf) {
		int pos = vcf.getPosition();
		String refBases = vcf.getReference();
		String mutBases = vcf.getAlternate()[0];
		//make layout
		SamLayoutLite sl = new SamLayoutLite(r);
		//attempt mod
		if (sl.markDeletion(pos, refBases, mutBases, warn)){
			String[] sq = sl.getSequenceAndQualtities();
			r.setReadString(sq[0]);
			r.setBaseQualityString(sq[1]);
			r.setAttribute("BB", "DEL_"+pos+"_"+refBases+"_"+mutBases);
			return true;
		}
		
		if (warn) {
			System.err.println(vcf);
			System.err.println(r.getSAMString().trim());
		}
		return false;
	}

	private boolean addSnv(SAMRecord r, VCFRecord vcf) {
		int pos = vcf.getPosition();
		char refBase = vcf.getReference().charAt(0);
		char mutBase = vcf.getAlternate()[0].charAt(0);
		//make layout
		SamLayoutLite sl = new SamLayoutLite(r);
		//attempt mod
		if (sl.changeBase(pos, refBase, mutBase, warn)){
			r.setReadString(new String(sl.getSequence()));
			r.setAttribute("BB", "SNV_"+pos+"_"+refBase+"_"+mutBase);
			return true;
		}
		else if (warn) {
			System.err.println(vcf);
			System.err.println(r.getSAMString().trim());
		}
		return false;
	}




	/**Checks that the alignment actually touches down on at least one base of the region to avoid spanners and masked alignments.
	 * @throws IOException */
	public ArrayList<SAMRecord> fetchAlignments (VCFRecord vcf) throws IOException{
		ArrayList<SAMRecord> al = new ArrayList<SAMRecord>();

		int[] interbaseStartStop = QueryIndexFileLoader.fetchEffectedBps(Misc.TAB.split(vcf.getOriginalRecord()), false);
		if (interbaseStartStop == null) throw new IOException("WARNING: Failed to parse the effected bps for : "+vcf.getOriginalRecord());
		int start = interbaseStartStop[0];
		int end = interbaseStartStop[1];

		
		SAMRecordIterator i = bamReader.queryOverlapping(workingChrom, (start+1), end); 
		while (i.hasNext()) {
			SAMRecord sam = i.next();
			//unmapped? not primary? 
			if (sam.getReadUnmappedFlag() || sam.isSecondaryOrSupplementary()) continue;
			//fetch blocks of actual alignment
			ArrayList<int[]> blocks = fetchAlignmentBlocks(sam.getCigarString(), sam.getUnclippedStart()-1);
			//check to see if any intersect the region, needed since the queryOverlap returns spanners
			for (int[] b : blocks){
				if (end < b[0] || start > b[1]) continue;
				al.add(sam);
				break;
			}
		}
		i.close();
		return al;
	}
	
	private String getQualifiedName(SAMRecord sam){
		StringBuilder sb = new StringBuilder(sam.getReadName());
		//primary or secondary
		if (sam.getNotPrimaryAlignmentFlag()) sb.append("_S");
		else sb.append("_P");
		//unpaired, paired1, paired2
		if (sam.getReadPairedFlag()== false) sb.append("3");
		else {
			if (sam.getFirstOfPairFlag()) sb.append("1");
			else sb.append("2");
		}
		//strand, prob isn't necessary
		if (sam.getReadNegativeStrandFlag()) sb.append("-");
		else sb.append("+");
		return sb.toString();
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BamBlaster(args);
	}

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bamFile = new File(args[++i]); break;
					case 'r': saveDirectory = new File(args[++i]); break;
					case 'v': vcfFile = new File(args[++i]); break;
					case 's': maxSizeIndel = Integer.parseInt(args[++i]); break;
					case 'm': minVarDist = Integer.parseInt(args[++i]); break;
					case 'd': minAlignmentDepth = Integer.parseInt(args[++i]); break;
					case 'i': warn = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//check save dir
		if (saveDirectory == null) Misc.printErrAndExit("\nCan't find your save directory? Aborting!\n");
		saveDirectory.mkdirs();
		
		//parse vcf file
		if (vcfFile == null || vcfFile.exists() == false) Misc.printErrAndExit("\nCan't find or read your vcf file? Aborting!\n"+vcfFile);
		VCFParser vcfParser = new VCFParser(vcfFile, true, false, false);
		chromVariant = vcfParser.getChromosomeVCFRecords();
		if (chromVariant == null || chromVariant.size() == 0)  Misc.printErrAndExit("\nCan't parse variants from your vcf file?\n"+vcfFile);

		//check bam file
		if (bamFile == null || bamFile.exists() == false) Misc.printErrAndExit("\nCan't find your bam file? Aborting!\n"+bamFile);
		readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		bamReader = readerFactory.open(bamFile);
		if (bamReader.hasIndex() == false){
			try {
				bamReader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			Misc.printErrAndExit("\nCan't find your bam index? Aborting!\n"+bamFile);
		}
		
		//make files
		String name = Misc.removeExtension(bamFile.getName());
		tempBam = new File (saveDirectory, name+"_temp.bam");
		tempBam.deleteOnExit();
		File finalBam = new File (saveDirectory, name+"_filtered.bam");
		File unmodifiedFastqBam = new File (saveDirectory, name+"_unmodified.bam");
		
		//make File just to delete it, this is created by Picard's sort sam
		File tempSortedBai = new File (saveDirectory, name+"_temp_sorted.bai");
		tempSortedBai.deleteOnExit();
		
		//make writers
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		writerFactory.setTempDirectory(saveDirectory);
		tempBamWriter = writerFactory.makeBAMWriter(bamReader.getFileHeader(), false, tempBam);
		writerFactory.setCreateIndex(true);
		finalBamWriter = writerFactory.makeBAMWriter(bamReader.getFileHeader(), true, finalBam);
		unmodifiedFastqBamWriter = writerFactory.makeBAMWriter(bamReader.getFileHeader(), true, unmodifiedFastqBam);
		
		//make gzippers and log writer
		try {
			fastqFirstOut = new Gzipper(new File(saveDirectory, name+"_BB_1.fastq.gz"));
			fastqSecondOut = new Gzipper(new File(saveDirectory, name+"_BB_2.fastq.gz"));
			fastqUnpairedOut = new Gzipper(new File(saveDirectory, name+"_BB_unpaired.fastq.gz"));
			String vcfName = Misc.removeExtension(vcfFile.getName());
			excludedVCF = new Gzipper(new File(saveDirectory, vcfName+"_excluded.vcf.gz"));
			excludedVCF.println(vcfParser.getStringComments());
			includedVCF = new Gzipper(new File(saveDirectory, vcfName+"_included.vcf.gz"));
			includedVCF.println(vcfParser.getStringComments());
			//make writer for variant info
			File varReport = new File(saveDirectory, name+"_varRep.txt");
			varReportOut = new PrintWriter (new FileWriter(varReport));
			varReportOut.println("#ModifiedAlignments\t#OverlappingAlignments\tVCFRecord");
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem with creating gzippers? Aborting!\n");
		}
	}	
	
	/**Assumes interbase coordinates for start and returned blocks.*/
	public static ArrayList<int[]> fetchAlignmentBlocks(String cigar, int start){
		//for each cigar block
		Matcher mat = CIGAR_SUB.matcher(cigar);
		ArrayList<int[]> blocks = new ArrayList<int[]>();
		while (mat.find()){
			String call = mat.group(2);
			int numberBases = Integer.parseInt(mat.group(1));
			//a match
			if (call.equals("M")) {
				blocks.add(new int[]{start, start+numberBases});
			}
			//just advance for all but insertions which should be skipped via the failure to match
			start += numberBases;
		}
		return blocks;
	}

	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               Bam Blaster : May 2020                            **\n" +
				"**************************************************************************************\n" +
				"Injects SNVs and INDELs from a vcf file into bam alignments. These and their mates are\n"+
				"extracted as fastq for realignment. For SNVs, only alignment bases that match the\n"+
				"reference and have a CIGAR of M are modified. Not all alignments can be modified.\n"+
				"Secondary/supplemental/not proper are skipped. One var per alignment. Variants within\n"+
				"read length distance of prior are ignored and saved to file for iterative processing.\n"+
				"Be sure to normalize and decompose your vcf file (e.g.https://github.com/atks/vt).\n" +
				"INDELs first base must be reference. Use the ExactBamMixer or BamMixer to add\n"+
				"realignments (e.g. 10%) with the unmodified.bams (e.g. 90%). Use the VCFVariantMaker\n"+
				"to generate random vcf variants or pull a VCF from Clinvar/ Cosmic.\n\n"+

				"Required:\n"+
				"-b Path to a coordinate sorted bam file with index.\n"+
				"-v Path to a trimmed, normalized, decomposed vcf variant file, zip/gz OK.\n"+
				"-r Full path to a directory to save the results.\n" +
				"-s Max size INDEL, defaults to 50\n"+
				"-d Min alignment depth, defaults to 25\n"+
				"-m Min distance between variants, defaults to 150\n"+
				"-i Verbose debugging output\n"+

				"\nExample: java -Xmx10G -jar pathTo/USeq/Apps/BamBlaster -b ~/BMData/na12878.bam\n"+
				"    -r ~/BMData/BB0 -v ~/BMData/clinvar.pathogenic.SnvIndel.vcf.gz \n\n" +

				"**************************************************************************************\n");
	}
}
