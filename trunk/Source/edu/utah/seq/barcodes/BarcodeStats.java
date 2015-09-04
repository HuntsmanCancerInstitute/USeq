package edu.utah.seq.barcodes;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.parsers.BarParser;
import edu.utah.seq.useq.ArchiveInfo;
import edu.utah.seq.useq.SliceInfo;
import edu.utah.seq.useq.USeqUtilities;
import edu.utah.seq.useq.data.*;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import util.bio.annotation.Bed;
import util.bio.annotation.ExportIntergenicRegions;
import util.gen.*;

/**Counts stats on barcodes
 * @author Nix
 * */
public class BarcodeStats {

	//user fields
	private File[] samFiles;
	private Histogram histogram;
	private long numAlign1;
	private long numAlign2Plus;
	private int minFamilySizeForSaving = 10;
	private String familyTag = "FZ";
	private SAMFileWriter passingBamWriter;

	/**For stand alone app.*/
	public BarcodeStats(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();
		
		//process args
		processArgs(args);

		
		parseAlignmentFiles();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");
	}


	
	

	public void parseAlignmentFiles(){
		try {
			//make reader and writer factories
			SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
			writerFactory.setCreateIndex(true);
			writerFactory.setTempDirectory(samFiles[0].getParentFile());
			
			System.out.println("Processing:");
			
			for (File samFile: samFiles){
				System.out.print(samFile.getName()+"\t");
				//reset counters
				histogram = new Histogram(0, 100, 100);
				numAlign1 = 0;
				numAlign2Plus = 0;
				
				//make reader and writer
				SamReader samReader = factory.open(samFile);
				File pass = new File (samFile.getParentFile(), "pass"+familyTag+minFamilySizeForSaving+"Min_"+samFile.getName());
				passingBamWriter = writerFactory.makeBAMWriter(samReader.getFileHeader(), true, pass);
				
				//iterate through records saving those to bam if pass minFamilySize
				SAMRecordIterator it = samReader.iterator();
				parse(it);
				
				//cleanup
				it.close();
				samReader.close();
				passingBamWriter.close();
				
				//stats
				System.out.println(numAlign1+"\t"+numAlign2Plus);
				histogram.printScaledHistogram();
				System.out.println();
			}

		} catch (Exception e){
			e.printStackTrace();
			System.exit(1);
		}

	}

	

	private void parse(SAMRecordIterator it) {
		while (it.hasNext()){
			SAMRecord sam = it.next();
			Object fz = sam.getAttribute(familyTag);
			if (fz == null) Misc.printErrAndExit("Problem finding FZ tag for "+sam.getSAMString().trim());
			int count = (Integer)fz;
			histogram.count(count);
			if (count == 1) numAlign1++;
			else numAlign2Plus++;
			if (count >= minFamilySizeForSaving) passingBamWriter.addAlignment(sam);
		}
		
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BarcodeStats(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': forExtraction = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		File[][] tot = new File[4][];
		tot[0] = IO.extractFiles(forExtraction,".sam");
		tot[1] = IO.extractFiles(forExtraction,".sam.gz");
		tot[2] = IO.extractFiles(forExtraction,".sam.zip");
		tot[3] = IO.extractFiles(forExtraction,".bam");
		samFiles = IO.collapseFileArray(tot);
		if (samFiles == null || samFiles.length ==0 || samFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.sam(.zip/.gz) or xxx.bam file(s)!\n");
			

	}	

	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Barcode Stats : Aug 2015                          **\n" +
				"**************************************************************************************\n" +
				"Beta\n\n" +

				"Required Options:\n"+
				"-f Full path to a sorted bam file or directory containing such. Multiple files are\n"+
				"    processed independantly.\n" +

				
				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/BarcodeStats -f /Data/SamFiles/ \n\n"+

				"**************************************************************************************\n");

	}
		

}
