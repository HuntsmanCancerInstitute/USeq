package edu.utah.seq.base;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import net.sf.samtools.*;

import util.gen.*;

/**Collects summary statistics on nmers related to per base intensities for all 4 channels at each position.
 * Must run the modified Picard IlluminaBaseCallsToSam app on the pipeline output to place the intensities in each record.
 * @author Nix
 * */
public class BamNMerIntensityParser {

	//user fields
	private File[] samFiles;
	private File resultsFile;
	private int nMerSize = 5;
	private byte minimumBaseQuality = 20;

	//internal
	private short[] aIntensities;
	private short[] cIntensities;
	private short[] gIntensities;
	private short[] tIntensities;
	private HashMap<String, NMer> nMers = new HashMap<String, NMer>();
	private Pattern badBases = Pattern.compile(".*[^GATC].*");
	private long numberReadsProcessed = 0;


	//constructor
	public BamNMerIntensityParser(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();
		//process args
		processArgs(args);

		//split sam files by nmer
		System.out.println("Parsing SAM/BAM files by N-mer...");
		parseFiles();

		//write out nmer data
		printNMers();
		
		//save NMers hash?
		
		System.out.println(numberReadsProcessed+"\tReads processed.\n");

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	public void printNMers(){
		try {
			//make results writer
			if (resultsFile == null){
				String name;
				name = Misc.removeExtension(samFiles[0].getCanonicalPath()) + "_"+nMerSize+"Mers.txt.gz";
				resultsFile = new File (name);
			}
			Gzipper gzip = new Gzipper(resultsFile);
			
			//add header
			StringBuilder sb = new StringBuilder();
			sb.append("NMer");
			//for every position in the nMer
			for (int i=0; i<nMerSize; i++){
				//A
				sb.append("\t");
				sb.append(i);
				sb.append("_A_Mean");
				sb.append("\t");
				sb.append(i);
				sb.append("_A_SD");
				//C
				sb.append("\t");
				sb.append(i);
				sb.append("_C_Mean");
				sb.append("\t");
				sb.append(i);
				sb.append("_C_SD");
				//G
				sb.append("\t");
				sb.append(i);
				sb.append("_G_Mean");
				sb.append("\t");
				sb.append(i);
				sb.append("_G_SD");
				//T
				sb.append("\t");
				sb.append(i);
				sb.append("_T_Mean");
				sb.append("\t");
				sb.append(i);
				sb.append("_T_SD");
			}
			gzip.println(sb.toString());
			
			//print nmers
			for (String seq: nMers.keySet()){
				sb = new StringBuilder();
				NMer nmer = nMers.get(seq);
				sb.append(seq);
				sb.append("\t");
				sb.append(nmer.toString());
				gzip.println(sb.toString());
			}
			gzip.close();
			
		} catch (IOException e) {
			System.out.println("\nError writing out results!\n");
			e.printStackTrace();
			System.exit(1);
		}
	}



	public void parseFiles(){
		for (File samFile: samFiles){
			System.out.print("\t"+samFile.getName());

			SAMFileReader samReader = null;
			int counter =0;

			try {
				samReader = new SAMFileReader(samFile);
				SAMRecordIterator it = samReader.iterator();

				while (it.hasNext()) {
					SAMRecord sam = it.next();

					//print status blip
					if (++counter == 200000){
						System.out.print(".");
						counter = 0;
					}
					
					numberReadsProcessed++;

					//fetch data
					String sequence = sam.getReadString();
					aIntensities = sam.getSignedShortArrayAttribute("IA");
					cIntensities = sam.getSignedShortArrayAttribute("IC");
					gIntensities = sam.getSignedShortArrayAttribute("IG");
					tIntensities = sam.getSignedShortArrayAttribute("IT");
					byte[] qualities = sam.getBaseQualities();

					//count each nmer
					int end = sequence.length()-nMerSize+1;
					for (int i=0; i< end; i++){
						int endNMer = i+nMerSize;
						//fetch seq
						String nMerSeq = sequence.substring(i, endNMer);

						//look for bad qualities
						if (minimumBaseQuality !=0){
							for (int x=i; x< endNMer; x++){
								if (qualities[x] < minimumBaseQuality) {
									//System.out.println("\nBadBase "+nMerSeq);
									continue;
								}
							}
						}
						
						//look for bad bases
						Matcher mat = badBases.matcher(nMerSeq);
						if (mat.matches()) continue;

						//attempt to fetch old
						NMer nMer = nMers.get(nMerSeq);
						if (nMer == null){
							nMer = new NMer((byte)nMerSeq.length());
							nMers.put(new String(nMerSeq), nMer);
						}

						//count intensities
						nMer.count(i, endNMer, aIntensities, cIntensities, gIntensities, tIntensities);
					}


				}
				System.out.println();
			} catch (Exception e){
				System.err.println("\nError parsing sam file!\n");
				e.printStackTrace();
				System.exit(1);
			}
			System.out.println();
		}

	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BamNMerIntensityParser(args);
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
					case 'f': forExtraction = new File(args[++i]); break;
					case 'q': minimumBaseQuality = Byte.parseByte(args[++i]); break;
					case 'r': resultsFile = new File(args[++i]); break;
					case 'n': nMerSize = Integer.parseInt(args[++i]); break;
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
		if (samFiles == null || samFiles.length ==0 || samFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.sam(.zip/.gz) file(s)!\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Bam NMer Intensity Parser : April 2012                 **\n" +
				"**************************************************************************************\n" +
				"Parses a BAM file from a modified Picard IlluminaBaseCallsToSam run on raw Illumina\n" +
				"sequencing data to extract information regarding N mers.\n\n" +

				"Options:\n"+
				"-f Full path to a bam file or directory containing such. Multiple files are merged.\n"+
				"-r Full path file name to save results, defaults to a derivative of -f\n"+
				"-n Length of the N mer, defaults to 5.\n" +
				"-q Minimum base quality score, defaults to 20. Only N-mers where all bases pass the\n"+
				"       threshold are scored.\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/BamIntensityParser -f /Data/BamFiles/\n" +
				"       -n 7 -q 30\n\n"+

		"**************************************************************************************\n");

	}


}
