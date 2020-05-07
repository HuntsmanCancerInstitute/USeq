package edu.utah.seq.vcf.sim;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.query.QueryIndexFileLoader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Pulls mpileup lines that overlap each vcf record and calculates AF and z-score stats.*/
public class BamMixerLoader implements Runnable{

	//fields
	private boolean failed = false;
	private SamReader modifiedReader = null;
	private SamReader unModifiedReader = null;
	private ArrayList<String> vcfLines = new ArrayList<String>();
	private ExactBamMixer ebm = null;
	private Gzipper[] samWriters = null;
	private Gzipper[] vcfWriters = null;
	private Gzipper targetVarResults;
	private int loaderId;
	private int minNumAltReads = 0;

	private static final Pattern BB = Pattern.compile(":BB-");
	private static final Pattern numUnder = Pattern.compile("^\\d[\\d_]+");
	private static final Pattern trailingBB = Pattern.compile(":BB-");
	private Random random = new Random(1); 
	private double[] targetFractions = null;

	public BamMixerLoader (File modifiedBam, File unModifiedBam, ExactBamMixer ebm, int loaderId) throws IOException{
		this.ebm = ebm;
		this.loaderId = loaderId;
		minNumAltReads = ebm.getMinNumAltReads();

		//make alignment readers
		modifiedReader = ebm.getReaderFactory().open(modifiedBam);
		unModifiedReader = ebm.getReaderFactory().open(unModifiedBam);

		//make writers
		targetFractions = ebm.getFractions();
		samWriters = new Gzipper[targetFractions.length];
		vcfWriters = new Gzipper[targetFractions.length];
		for (int i=0; i< samWriters.length; i++) {
			File samOut = new File (ebm.getSaveDirectory(), targetFractions[i]+"_"+loaderId+".sam.gz");		
			samOut.deleteOnExit();
			samWriters[i] = new Gzipper(samOut);
			File vcfOut = new File (ebm.getSaveDirectory(), targetFractions[i]+"_"+loaderId+".vcf.gz");
			vcfOut.deleteOnExit();
			vcfWriters[i] = new Gzipper(vcfOut);
		}
		targetVarResults = new Gzipper(new File(ebm.getSaveDirectory(), loaderId+"_VarRes.txt.gz"));
		targetVarResults.getGzipFile().deleteOnExit();

		//save headers?
		if (loaderId == 0) {
			//sam
			String samHeader = unModifiedReader.getFileHeader().getSAMString();
			for (int i=0; i< targetFractions.length; i++) {
				samWriters[i].print(samHeader);
				targetVarResults.print("TargetAF\tActualAF\t(Mod/UnModAlignPairs)\t");
			}
			targetVarResults.println("VcfRecord");
			//vcf
			for (Gzipper g: vcfWriters) g.println(ebm.getVcfHeader());
		}

	}



	public void run() {	
		try {
			
			//get next chunk of work
			while (ebm.loadVcfRecords(vcfLines)){ 

				//for each record
				for (String record: vcfLines) processRecord(record);
				
				//cleanup
				vcfLines.clear();
			}
		} catch (Exception e) {
			failed = true;
			e.printStackTrace();
			Misc.printErrAndExit("\nError in loader "+loaderId);

		} finally {
			try {
				//close readers
				modifiedReader.close();
				unModifiedReader.close();
				//close writers
				for (Gzipper g: samWriters) g.close();
				for (Gzipper g: vcfWriters) g.close();				
				targetVarResults.close();
			} catch (IOException e){
				failed = true;
				e.printStackTrace();
				Misc.printErrAndExit("\nError closing readers or writers for loader "+loaderId);
			}
		}
	}

	private void processRecord(String record) throws IOException {
		//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR
		String[] vcfFields = Misc.TAB.split(record);

		int[] interbaseStartStop = QueryIndexFileLoader.fetchEffectedBps(vcfFields, false);
		if (interbaseStartStop == null) throw new IOException("WARNING: Failed to parse the effected bps for : "+record);

		//fetch modified alignments, not all will be truly modified
		HashMap<String, ArrayList<SAMRecord>> modifiedPairs = fetchAlignmentPairs(modifiedReader, interbaseStartStop, vcfFields[0], true);
		String[] modNames = new String[modifiedPairs.size()];
		int index =0;
		for (String key : modifiedPairs.keySet()) modNames[index++] = key;
		Misc.randomize(modNames, random);	

		//fetch unModifiedAlignments
		HashMap<String, ArrayList<SAMRecord>> unModifiedPairs = fetchAlignmentPairs(unModifiedReader, interbaseStartStop, vcfFields[0], false);
		double numUnmodified = unModifiedPairs.size();
		ArrayList<String> readPairNamesWritten = new ArrayList<String>();
		//for each target fraction
		for (int i=0; i< targetFractions.length; i++) {

			//calculate number of modified pairs to add
			int numMod = (int)Math.round(targetFractions[i] * numUnmodified);
			if (numMod < minNumAltReads) {
				targetVarResults.print(targetFractions[i]+"\t0\tSkippedTooFewMod(" +numMod+ "/"+ (int)numUnmodified+")\t");
			}
			else if (numMod > modifiedPairs.size()) {
				targetVarResults.print(targetFractions[i]+"\t0\tSkippedNotEnoughMod(" +modifiedPairs.size()+ "/"+ (int)numUnmodified+")\t");
			}
			else {
				double actAF = (double)numMod / (numUnmodified);
				targetVarResults.print(targetFractions[i]+"\t"+ Num.formatNumber(actAF, 4)+ "\t(" +numMod+ "/"+ (int)numUnmodified+")\t");
				vcfWriters[i].println(record);
				readPairNamesWritten.clear();

				//write out modified pairs, randomized
				for (int m=0; m< numMod; m++) {
					ArrayList<SAMRecord> al = modifiedPairs.get(modNames[m]);
					readPairNamesWritten.add(modNames[m]);
					for (SAMRecord sam: al) samWriters[i].print(sam.getSAMString());
				}
				
				//write out unModified minus those already written out
				Set<String> keys = unModifiedPairs.keySet();			
				keys.removeAll(readPairNamesWritten);
				for (String rm: keys) {			
					ArrayList<SAMRecord> al = unModifiedPairs.get(rm);
					for (SAMRecord sam: al) samWriters[i].print(sam.getSAMString()); 
				}	
			}
		}
		//close line
		targetVarResults.println(record);
	}



	public HashMap<String, ArrayList<SAMRecord>> fetchAlignmentPairs(SamReader reader, int[] interbaseStartStop, String chr, boolean modified) {			
		//use hash to collect paired reads, only want one
		HashMap<String, ArrayList<SAMRecord>> nameRecord = new HashMap<String, ArrayList<SAMRecord>>();

		int start = interbaseStartStop[0];
		int end = interbaseStartStop[1];
		
		//won't fetch mate that doesn't intersect coordinates!
		SAMRecordIterator i = reader.queryOverlapping(chr, (start+1), end); 

		while (i.hasNext()) {
			SAMRecord sam = i.next();
			//unmapped? not primary? 
			if (sam.getReadUnmappedFlag() || sam.isSecondaryOrSupplementary() ) continue;

			//fetch blocks of actual alignment
			ArrayList<int[]> blocks = BamBlaster.fetchAlignmentBlocks(sam.getCigarString(), sam.getUnclippedStart()-1);

			//check to see if any intersect the region, needed since the queryOverlap returns spanners
			for (int[] b : blocks){
				if (end < b[0] || start > b[1]) continue;

				//OK overlaps so add it
				//clean up name and add tag?  Only with modified
				if (modified) {
					stripNameNumber(sam);
					String[] tokens = BB.split(sam.getReadName());
					if (tokens.length != 1) sam.setAttribute("BB", tokens[1]);
				}

				ArrayList<SAMRecord> al = nameRecord.get(sam.getReadName());
				if (al == null) {
					al = new ArrayList<SAMRecord>();
					nameRecord.put(sam.getReadName(), al);
				}
				al.add(sam);
				break;
			}
		}
		i.close();
		return nameRecord;
	}

	/**Removes leading and trailing info from BamBlaster
	 * e.g. 0_HWI-D00294:322:CATY4ANXX:6:1101:1166:34209:BB-INS_178536299_C_CT_2  ->  HWI-D00294:322:CATY4ANXX:6:1101:1166:34209 */
	public static void stripNameNumber(SAMRecord sam){
		String name = sam.getReadName();
		int start = 0; 
		int stop = name.length();

		Matcher mat = numUnder.matcher(name);
		if (mat.find()) start = mat.end();

		Matcher trailing = trailingBB.matcher(name);
		if (trailing.find()) stop = trailing.start();

		String newName = name.substring(start, stop);
		sam.setReadName(newName);
	}
	public boolean isFailed() {
		return failed;
	}
	public Gzipper[] getSamWriters() {
		return samWriters;
	}
	public Gzipper getTargetVarResults() {
		return targetVarResults;
	}
	public Gzipper[] getVcfWriters() {
		return vcfWriters;
	}

}
