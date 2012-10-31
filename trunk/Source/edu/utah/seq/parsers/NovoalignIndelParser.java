
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import edu.expr.*;
import util.gen.*;
import java.util.*;
import util.bio.seq.Seq;

/**Parses a Novoalign alignment txt file for indels creating a consensus allele table.
 * Final positions are in interbase coordinates (0 start, stop excluded).
 * @author david.nix@hci.utah.edu 
 **/
public class NovoalignIndelParser {
	//fields
	private File[] dataFiles;
	private File resultsDirectory;
	private Pattern tab = Pattern.compile("\\t");
	private Pattern whiteSpace = Pattern.compile("\\s+");
	private Pattern plus = Pattern.compile("\\+");
	private Pattern minus = Pattern.compile("-");
	private float minimumPosteriorProbabilityThreshold = 13;
	private float minimumBaseQualityThreshold = 13;
	private int minimumUniqueReadsThreshold = 2;
	private HashMap<String,Allele> alleles = new HashMap<String,Allele>();
	private boolean skipUnderScoredChromosomes = true;
	private boolean filterForChrLines = false;

	//constructors
	public NovoalignIndelParser(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);

		System.out.println("\nParsing and filtering...");
		//for each file, parse and load alleles HashMap
		for (int i=0; i< dataFiles.length; i++){
			System.out.println("\t"+dataFiles[i].getName());
			parseFile(dataFiles[i]); 
		}

		//print results
		if (alleles.size() !=0) printData();
		else System.out.println("\nNo indels found passing thresholds.\n");
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	public void printData(){
		try{
			//alleler
			File allelerFile = new File (resultsDirectory, "indels.alleler");
			PrintWriter alleler = new PrintWriter( new FileWriter(allelerFile));
			alleler.println("#MinPostProb = "+minimumPosteriorProbabilityThreshold+"\n#MinBaseQualities = "+minimumBaseQualityThreshold+"\n#MinNumberUniqueAlignments = "+minimumUniqueReadsThreshold);
			alleler.println("#Chr\tStart\tStop\tINDEL\tPosteriorProbabilities\tBaseScores\t#OverlappingAlignments\t#UniqueOverlappingAlignments");
			
			//insertions
			File insertionsFile = new File (resultsDirectory, "insertions.bed");
			PrintWriter insertions = new PrintWriter( new FileWriter(insertionsFile));
			insertions.println("#MinPostProb = "+minimumPosteriorProbabilityThreshold+"\n#MinBaseQualities = "+minimumBaseQualityThreshold+"\n#MinNumberUniqueAlignments = "+minimumUniqueReadsThreshold);
			insertions.println("#Chr\tStart-1\tStop+1\tbasesInserted\t#UniqueOverlappingAlignments\tBlank");
		
			File deletionFile = new File (resultsDirectory, "deletions.bed");
			PrintWriter deletions = new PrintWriter( new FileWriter(deletionFile));
			deletions.println("#MinPostProb = "+minimumPosteriorProbabilityThreshold+"\n#MinBaseQualities = "+minimumBaseQualityThreshold+"\n#MinNumberUniqueAlignments = "+minimumUniqueReadsThreshold);
			deletions.println("#Chr\tStart\tStop\tBlank\t#UniqueOverlappingAlignments\tBlank");
			Iterator<String> it = alleles.keySet().iterator();
			while (it.hasNext()){
				String name = it.next();
				//check number of unique reads
				Allele allele = alleles.get(name);
				String details = allele.alleleTableLine(minimumUniqueReadsThreshold);
				if (details !=null){
					String chrom = tab.split(name)[0];
					alleler.println(chrom+"\t"+details);
					String in = allele.fetchInsertionBedLine();
					if (in != null) insertions.println(chrom+"\t"+in);
					else {
						String del = allele.fetchDeletionBedLine();
						deletions.println(chrom+"\t"+del);
					}
				}
			}
			alleler.close();
			insertions.close();
			deletions.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public boolean parseFile(File workingFile){
		try{
			//get reader
			BufferedReader in = IO.fetchBufferedReader(workingFile);
			String line;
			String[] tokens = null;

			/*read in lines
			text type read qual matchClass alignScore pprob chrom position ori
			 0    1     2   3     4            5        6    7       8     9  
			                                                     
			 */
			int counter = 0;
			while ((line = in.readLine()) !=null){
				if (++counter == 25000){
					System.out.print(".");
					counter = 0;
				}
				if (line.startsWith("#")) continue;
				tokens = tab.split(line);
				if (tokens.length < 14) continue;
				//look for indels
				if (tokens[13].contains("-") == false && tokens[13].contains("+") == false) continue;
				//check probability threshold
				float probScore = Float.parseFloat(tokens[6]);
				if (probScore < minimumPosteriorProbabilityThreshold) continue;
				//parse chromosome
				String chromosome = tokens[7].substring(1);
				//any underscores designating this as a splice junction match?
				if (skipUnderScoredChromosomes && chromosome.contains("_")) continue;
				//parse strand
				String strand = "+";
				if (tokens[9].equals("R")) strand = "-";
				//convert quality scores
				int[] baseQualities = Seq.convertScores(tokens[3]);
				//coordinate in 1base
				int coordinate = Integer.parseInt(tokens[8]);
				//parse modification(s), these are single base insertions or deletions
				String[] modifications = whiteSpace.split(tokens[13]);
				Allele allele = null;
				for (int i=0; i< modifications.length; i++){
					//insertion?
					if (modifications[i].contains("+")){
						String[] posBase = plus.split(modifications[i]);
						int base = Integer.parseInt(posBase[0]);
						//check score
						int score;
						if (strand.equals("+")) score = baseQualities[base-1];
						else score = baseQualities[baseQualities.length-base];
						if (score < minimumBaseQualityThreshold) continue;
						//make start/stop
						int startStop = coordinate+base-2;
						//make text
						String name = chromosome+"\t"+ startStop+posBase[1];
						//fetch or make allele
						if (alleles.containsKey(name)) allele = alleles.get(name);
						else {
							allele = new Allele(startStop, startStop, posBase[1], null);
							allele.setBaseQualities(new StringBuilder());
							allele.setPosteriorProbabilites(new StringBuilder());
							allele.setReadStartPositions(new StringBuilder());
							alleles.put(name, allele);
						}
						//add
						allele.getBaseQualities().append(score+",");
						allele.getPosteriorProbabilites().append(probScore+",");
						allele.getReadStartPositions().append(strand+coordinate+",");
					}
					//deletion
					else if (modifications[i].contains("-")){
						String[] posBase = minus.split(modifications[i]);
						int base = Integer.parseInt(posBase[0]);
						//check flanking scores with respect to plus strand
						int scoreL;
						int scoreR;
						if (strand.equals("+")) {
							scoreL = baseQualities[base-2];
							scoreR = baseQualities[base-1];
						}
						else {
							scoreL = baseQualities[baseQualities.length-base+2];
							scoreR = baseQualities[baseQualities.length-base+1];
						}
						if (scoreL < minimumBaseQualityThreshold || scoreR < minimumBaseQualityThreshold) continue;
						//make start/stop
						int start = coordinate+base-2;
						int stop = start+1;
						//make text
						String name = chromosome+"\t"+ start+"_"+stop+posBase[1];
						//fetch or make allele
						if (alleles.containsKey(name)) allele = alleles.get(name);
						else {
							allele = new Allele(start, stop, "", null);
							allele.setBaseQualities(new StringBuilder());
							allele.setPosteriorProbabilites(new StringBuilder());
							allele.setReadStartPositions(new StringBuilder());
							alleles.put(name, allele);
						}
						//add
						allele.getBaseQualities().append(scoreL+"/"+scoreR+",");
						allele.getPosteriorProbabilites().append(probScore+",");
						allele.getReadStartPositions().append(strand+coordinate+",");
					}
				}
			}
			in.close();
			return true;
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new NovoalignIndelParser(args);
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
					case 'f': forExtraction = new File(args[++i]); break;
					case 'r': resultsDirectory = new File(args[++i]); break;
					case 'p': minimumPosteriorProbabilityThreshold = Float.parseFloat(args[++i]); break;
					case 'b': minimumBaseQualityThreshold = Float.parseFloat(args[++i]); break;
					case 'u': minimumUniqueReadsThreshold = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//dir
		if (resultsDirectory == null) Misc.printErrAndExit("\nError: cannot find your results directory?");
		resultsDirectory.mkdir();
		
		//pull files
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction,".txt");
		tot[1] = IO.extractFiles(forExtraction,".txt.zip");
		tot[2] = IO.extractFiles(forExtraction,".txt.gz");

		dataFiles = IO.collapseFileArray(tot);
		if (dataFiles == null || dataFiles.length==0) dataFiles = IO.extractFiles(forExtraction);
		if (dataFiles == null || dataFiles.length ==0 || dataFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.txt(.zip/.gz) file(s)!\n");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         Novoalign Indel Parser: June 2010                        **\n" +
				"**************************************************************************************\n" +
				"Parses Novoalign alignment xxx.txt(.zip/.gz) files for consensus indels, something\n" +
				"currently not supported by the maq apps. Generates a consensus indel allele file,\n" +
				"interbase coordinates, for running through the Alleler application. Also creates two\n" +
				"bed files for the insertions and deletions.\n" +

				"\nOptions:\n"+
				"-f The full path directory/file text of your Novoalign xxx.txt(.zip or .gz) file(s).\n" +
				"-r Full path directory for saving the results.\n"+
				"-p Minimum alignment posterior probability (-10Log10(prob)) of being incorrect,\n" +
				"      defaults to 13 (0.05). Larger numbers are more stringent.\n"+
				"-b Minimum effected indel base quality score(s), ditto, defaults to 13.\n" +
				"-u Minimum number of unique reads covering indel, defaults to 2.\n"+


				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/NovoalignIndelParser -f /Novo/Run7/\n" +
				"     -r /Novo/Run7/indelAlleleTable.txt -p 20 -b 20 -u 3 \n" +

		"**************************************************************************************\n");

	}	
}
