
package edu.utah.seq.data;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.io.*;
import java.util.regex.*;
import util.gen.*;


/**Converts hg19 to b37 by cutting off the chr from the chromosome name. Also swaps out the header. Doesn't do anything with non standard chroms.
 * @author david.nix@hci.utah.edu 
 **/
public class BamHg19ToB37Converter{
	//fields
	private File[] bamFiles;
	private File bamFileWithGoodHeader;
	
	
	
	//constructors
	public BamHg19ToB37Converter(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			doWork();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void doWork() throws IOException{
		//make IO
		SamReaderFactory readerFactory = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX).validationStringency(ValidationStringency.SILENT);
		//open file to pull header from
		SamReader header = readerFactory.open(bamFileWithGoodHeader);
		
		for (File sam: bamFiles){
			SamReader in = readerFactory.open(sam);
			File samFile = new File (sam.getParentFile(), "b37_"+Misc.removeExtension(sam.getName())+".sam.gz");
			Gzipper out = new Gzipper(samFile);
			String h = header.getFileHeader().getSAMString().trim();
			out.println(h);
			parseSam(in, out);
			in.close();
			out.close();
		}
		
		System.out.println();

	}

	public boolean parseSam(SamReader sr, Gzipper out){
		SAMRecord sam = null;
		try {
			SAMRecordIterator it = sr.iterator();
			int dotCounter = 0;
			while (it.hasNext()) {
				if (++dotCounter > 1000000){
					System.out.print(".");
					dotCounter = 0;
				}
				sam = it.next();
				
				//have to pull back to txt to avoid index headaches
				String samString = sam.getSAMString().trim();
				String[] tokens = Misc.TAB.split(samString);
				if (tokens[2].contains("M")) tokens[2] = "MT";
				else tokens[2] = tokens[2].replace("chr", "");
				
				out.print(tokens[0]);
				for (int i=1; i< tokens.length; i++) {
					out.print("\t");
					out.print(tokens[i]);
				}
				out.println();
			}

		} catch (Exception e) {
			if (sam != null) System.err.println("Problem processing -> "+sam.getSAMString());
			e.printStackTrace();
			return false;
		}
		return true;
	}



	

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BamHg19ToB37Converter(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		System.out.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'e': bamFileWithGoodHeader = new File(args[++i]); break;
					case 'b': bamFiles = IO.extractFiles(new File(args[++i]), "bam"); break; 
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//pull file
		if (bamFiles == null || bamFiles.length ==0) Misc.printErrAndExit("\nError: cannot find any hg19 xxx.bam files to convert?\n");
		if (bamFileWithGoodHeader == null ) Misc.printErrAndExit("\nError: cannot find your b37 bam file to use in swapping headers?\n");
		


	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Bam Hg19 to B37 Converter: Aug 2016                   **\n" +
				"**************************************************************************************\n" +
				"Cuts off the chr from each reference chromosome, converts chrM to MT, and swaps out\n"+
				"the header to convert hg19 alignments to b37 alignments.\n"+

				"\nOptions:\n"+
				"-b Bam files to covert to b37, a directory with such or a single file.\n" +
				"-e A bam file with a good b37 header to add to the converted hg19 alignments.\n"+

				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/BamHg19B37Converter -b . -e ~/b37.bam\n\n" +

				"**************************************************************************************\n");

	}	
}
