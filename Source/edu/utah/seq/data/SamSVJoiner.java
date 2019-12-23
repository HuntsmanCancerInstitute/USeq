package edu.utah.seq.data;

import java.io.*;
import java.util.regex.*;

import util.bio.parsers.MultiFastaParser;
import util.bio.seq.Seq;
import util.gen.*;

import java.util.*;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMSequenceRecord;
import edu.utah.seq.data.sam.PicardSortSam;
import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.useq.data.IntersectingRegions;
import edu.utah.seq.useq.data.Region;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class SamSVJoiner{
	//user defined fields
	private File bamSpanner;
	private File[] bamOthers;
	private File bamResult;

	//internal fields
	private HashMap<String, SAMRecord> spanners = new HashMap<String, SAMRecord>();
	private Gzipper samOut;

	//alignment counts for sam files
	private long numberSpanners = 0;
	private long numberSpannersAdded = 0;
	private long numberOthers = 0;

	//constructors
	public SamSVJoiner(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);

		//load spanners into memory, ugly
		System.out.println("Loading paired read spanner file into memory...");
		loadSpanners();

		//parsing other
		System.out.println("Parsing other alignment files and removing spanners that match read name...");
		parseOthers();
		
		//add remaining spanners
		System.out.println("Adding remaining spanners...");
		addSpanners();
		
		System.out.println("Sorting...");
		new PicardSortSam (samOut.getGzipFile(), bamResult);
		
		
		System.out.println("Stats:");
		System.out.println(numberOthers+"\t# Other Alignments");
		System.out.println(numberSpanners+"\t# Alignment spanners");
		System.out.println(numberSpannersAdded+"\t# Alignment spanners added to others");
		
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}
	
	private void addSpanners() {
		try {
			//write out all the remaining records
			for (SAMRecord sr: spanners.values()) samOut.println(sr.getSAMString().trim());
			numberSpannersAdded = spanners.size();
			//close the gzipper
			samOut.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
	}

	public void parseOthers(){
		try{
			SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			for (File bamFile: bamOthers){
				//make reader
				SamReader reader = factory.open(bamFile);	
				SAMRecordIterator iterator = reader.iterator();

				int counter =0;
				//for each record, should be only two alignments per read name 
				while (iterator.hasNext()){
					SAMRecord samRecord = iterator.next();
					String pair;
					if (samRecord.getFirstOfPairFlag()) pair = "1";
					else pair = "2";
					String name = samRecord.getReadName()+pair;
					//remove it from the hash if present
					spanners.remove(name);
					//print current record;
					samOut.println(samRecord.getSAMString().trim());
					numberOthers++;
					//print status blip
					if (++counter == 100000){
						System.out.print(".");
						counter = 0;
					}
				}
				reader.close();
			}

		} catch (Exception e){
			System.err.println("\nError parsing other alignment files.\n");
			e.printStackTrace();
			System.exit(1);
		}
	}

	

	public void loadSpanners(){
		try{
			//make reader
			SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamSpanner);
			SAMRecordIterator iterator = reader.iterator();
			
			//make gzipper and write header to it
			String sam = Misc.removeExtension(bamResult.getName()) +".sam.gz";
			File samRes = new File(bamResult.getParentFile(), sam);
			samRes.deleteOnExit();
			samOut = new Gzipper (samRes);
			samOut.println(reader.getFileHeader().getTextHeader().trim());
		
			int counter =0;
			//for each record, should be only two alignments per read name 
			while (iterator.hasNext()){
				SAMRecord samRecord = iterator.next();
				String pair;
				if (samRecord.getFirstOfPairFlag()) pair = "1";
				else pair = "2";
				String name = samRecord.getReadName()+pair;
				spanners.put(name, samRecord);
				numberSpanners++;
				//print status blip
				if (++counter == 100000){
					System.out.print(".");
					counter = 0;
				}
			}
			reader.close();
			System.out.println();
			
		} catch (Exception e){
			System.err.println("\nError parsing spanner alignment file.\n");
			e.printStackTrace();
			System.exit(1);
		}
	}
	

	
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SamSVJoiner(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'o': forExtraction = new File(args[++i]); break;
					case 's': bamSpanner = new File(args[++i]); break;
					case 'r': bamResult = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (bamResult == null || bamResult.getName().endsWith(".bam") == false) Misc.printErrAndExit("\nPlease provide a results file with the '.bam' extension.\n");
		
		File[][] tot = new File[4][];
		tot[0] = IO.extractFiles(forExtraction,".sam");
		tot[1] = IO.extractFiles(forExtraction,".sam.gz");
		tot[2] = IO.extractFiles(forExtraction,".sam.zip");
		tot[3] = IO.extractFiles(forExtraction,".bam");
		bamOthers = IO.collapseFileArray(tot);
		if (bamOthers == null || bamOthers.length ==0 || bamOthers[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.sam(.zip/.gz) or xxx.bam file(s)!\n");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                SamSVJoiner: Feb 2013                           **\n" +
				"**************************************************************************************\n" +
				"Merges the -o alignments and uses their read names to filter the -s alignments for\n"+
				"those that don't match and appends these to the merge.\n"+

				"\nOptions:\n"+
				"-s Alignment file containing SAM records to parse (e.g. spanner -n30 alignments).\n"+
				"-o Alignment file or directory containing SAM/BAM files to be merged. xxx.sam(.gz/.zip)\n"+
				"     or xxx.bam are OK. Any sort order is OK. (e.g. soft, single, spanner alignments)\n" +
				"-r Merged bam results file.\n"+

				"\nExample: java -Xmx25G -jar pathToUSeq/Apps/SamSVJoiner -r merge.bam -s spanner.bam\n"+
				"     -o AllAligns/ \n\n"+


				"**************************************************************************************\n");

	}

}
