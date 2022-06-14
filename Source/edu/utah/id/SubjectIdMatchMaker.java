package edu.utah.id;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;


public class SubjectIdMatchMaker {
	
	//user defined fields
	private File subjectRegistryFile = null;
	private File querySubjectFile = null;
	private File matchResultsFile = null;
	private boolean addCoreIdsToRegistry = false;
	private boolean addQuerySubjectsToRegistry = false;
	private File updatedRegistry = null;
	
	//internal
	private Subject[] allSubjects = null;
	private Subject[] testSubjects = null;
	private int numberThreads = 0;
	private int minSubjectsPerChunk = 10;
	private MatcherEngine[] matchers = null;
	public int numberTopMatchesToReturn = 3;
	public double missingOneKeyPenalty = 0.1;
	public double missingAdditionalKeyPenalty = 1;
	private double maxEditScoreForMatch = 0.12;
	private HashMap<String,Subject> coreIdSubject = null;
	private CoreId coreId = new CoreId();
	
	public SubjectIdMatchMaker (String[] args) {
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);
			
			//load test subjects, will throw error if malformed
			IO.p("Loading test subjects to match against the registry... ");
			testSubjects = loadSubjectData(querySubjectFile, false);
			IO.pl(testSubjects.length);
			
			//load registry subjects
			IO.p("Loading registry... ");
			allSubjects = loadSubjectData(subjectRegistryFile, true);
			loadIdSubjectHash();
			IO.pl(allSubjects.length);

			//make a matcher for each chunk
			int numPerCore = fetchMinPerCore();

			Subject[][] split = chunk(allSubjects, numPerCore);
			matchers = new MatcherEngine[split.length];
			IO.pl("Launching "+split.length+" lookup threads...");
			for (int i=0; i< matchers.length; i++)  matchers[i] = new MatcherEngine(split[i], this);
			
			//run the comparison
			ExecutorService executor = Executors.newFixedThreadPool(matchers.length);
			for (MatcherEngine l: matchers) executor.execute(l);
			executor.shutdown();
			while (!executor.isTerminated()) {}

			//check loaders 
			for (MatcherEngine m: matchers) {
				if (m.isFailed()) throw new IOException("ERROR: Matcher engine issue! \n");
			}
			
			//for each test report best matches
			printResults();
			
			
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
			IO.pl("\nDone! "+Math.round(diffTime)+" Sec\n");
		} catch (Exception e) {
			IO.el("\nERROR running the SubjectIdMatchMaker, aborting. ");
			e.printStackTrace();
			matchResultsFile.delete();
			if (updatedRegistry!= null) updatedRegistry.delete();
			System.exit(1);
		}
	}
	
	
	
	
	private void loadIdSubjectHash() throws IOException {
		 coreIdSubject = new HashMap<String,Subject>();
		 for (int i=0; i< allSubjects.length; i++) {
			 String[] coreIds = allSubjects[i].getCoreIds();
			 if (coreIds == null) throw new IOException("ERROR: the registry subject["+i+"] has no coreId, see the -a option to create these." );
			 for (String id: coreIds) {
				 if (coreIdSubject.containsKey(id)) throw new IOException("ERROR: the registry subject["+i+"] has a duplicate coreId with a prior subject, see "+id );
				 coreIdSubject.put(id, allSubjects[i]);
			 }
		 }
	}

	private void printResults() throws IOException {
		
		//# OriginalSubject NumMatches Type Score lowestScoringSubject1 type score lowestScoringSubject2 type....
		//# NumMatches = number registry subjects with lowest score that are <= maxEditScore
		//# Type = Match or NonMatch
		//# Score = sum of the key's LevenshteinEditDistance/KeyLength + missingOneKeyPenalty for the first missing key + missingAdditionalKeyPenalty for each additional key  
		
		PrintWriter out = new PrintWriter( new FileWriter(matchResultsFile));
		out.println("# OriginalSubject\tNumMatches\tType1\tScore1\tLowestScoringRegistrySubject1\tType2\tScore2\tLowestScoringRegistrySubject2....");
		out.println("# NumMatches = number registry subjects with lowest score that are <= maxEditScore");
		out.println("# Type = Match or NonMatch");
		out.println("# Score = sum of the key's LevenshteinEditDistance/KeyLength + missingOneKeyPenalty for the first missing key + missingAdditionalKeyPenalty for each additional key");
		out.println("# Subject = LastName FirstName DobMonth DobDay DobYear Gender Mrn CoreIds HciPersonIds");
		//here
		
		for (Subject tp: testSubjects) {
			IO.pl();
			IO.pl("Orig\t"+ tp.toStringPretty());
			
			//set scores and sort
			Subject[] topMatches = tp.getTopMatches();
			double[] topMatchScores = tp.getTopMatchScores();
			for (int i=0; i< topMatches.length; i++) topMatches[i].setScore(topMatchScores[i]);
			Arrays.sort(topMatches);
			for (int i=0; i< topMatches.length; i++) IO.pl(Num.formatNumber(topMatches[i].getScore(), 3)+"\t"+topMatches[i].toStringPretty());
		}
		
		out.close();
	}

	private int fetchMinPerCore() {
		double numAllSubjects = allSubjects.length;
		for (int i=numberThreads; i >=1; i--) {
			int numPerChunk = (int)Math.round(numAllSubjects/(double)i);
			if (numPerChunk >= minSubjectsPerChunk) return numPerChunk;
		}
		return minSubjectsPerChunk;
	}

	private Subject[] loadSubjectData(File dataFile, boolean addCoreId) throws IOException {
		String line = null;
		BufferedReader in = IO.fetchBufferedReader(dataFile);
		ArrayList<Subject> pAL = new ArrayList<Subject>();
		while ((line = in.readLine())!= null) {
			if (line.length()==0 || line.startsWith("#"))continue;
			pAL.add(new Subject(line, addCoreId, coreId));
		}
		in.close();
		Subject[] p = new Subject[pAL.size()];
		pAL.toArray(p);
		return p;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SubjectIdMatchMaker(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		try {
			IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
			Pattern pat = Pattern.compile("-[a-zA-Z]");
			
			for (int i = 0; i<args.length; i++){
				Matcher mat = pat.matcher(args[i]);
				if (mat.matches()){
					char test = args[i].charAt(1);
					try{
						switch (test){
						case 'r': subjectRegistryFile = new File(args[++i]).getCanonicalFile(); break;
						case 'q': querySubjectFile = new File(args[++i]); break;
						case 'o': matchResultsFile = new File(args[++i]); break;
						case 't': numberThreads = Integer.parseInt(args[++i]); break;
						case 'm': numberTopMatchesToReturn = Integer.parseInt(args[++i]); break;
						case 'p': missingOneKeyPenalty = Double.parseDouble(args[++i]); break;
						case 'c': addCoreIdsToRegistry = true; break;
						case 'a': addQuerySubjectsToRegistry = true; break;
						case 'u': updatedRegistry = new File(args[++i]); break;
						case 's': maxEditScoreForMatch = Double.parseDouble(args[++i]); break;
						default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
						}
					}
					catch (Exception e){
						e.printStackTrace();
						Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
					}
				}
			}
			//check registry file
			if (subjectRegistryFile == null || subjectRegistryFile.canRead() == false) Misc.printErrAndExit("ERROR: failed to find the subject registry file "+subjectRegistryFile);
			
			//check query subject file
			if (querySubjectFile == null || querySubjectFile.canRead() == false) Misc.printErrAndExit("ERROR: failed to find the query subject file "+querySubjectFile);
			
			//check output file
			if (matchResultsFile == null) Misc.printErrAndExit("ERROR: failed to find the output match results file "+matchResultsFile);
			
			//check modifiers of registry
			if (addCoreIdsToRegistry || addQuerySubjectsToRegistry) {
				if (updatedRegistry==null) Misc.printErrAndExit("ERROR: looks like you've like to modify the registry but haven't provided a file path indicating where to save it? Complete -u .");
			}
			if (updatedRegistry!=null) {
				if (addCoreIdsToRegistry==false && addQuerySubjectsToRegistry==false) Misc.printErrAndExit("ERROR: looks like you'd like to save a modified registry but haven't indicated what modifications to make? Add -a and or -c .");
				if (updatedRegistry.getCanonicalPath().equals(subjectRegistryFile.getCanonicalPath())) Misc.printErrAndExit("ERROR: the exiting registry and updated registry file cannot be the same.");
			}
			//threads
			int numProc = Runtime.getRuntime().availableProcessors() - 1;
			if (numberThreads == 0 || numberThreads > numProc) numberThreads = numProc;		
			
			//print params
			printOptions();
			
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing arguments!");
		}
	}
	
	private void printOptions() throws IOException {
		String opt = "Options:\n"+
				"-r Registry file "+ subjectRegistryFile +"\n"+
				"-q Query file "+ querySubjectFile +"\n"+
				"-o Output match results "+ matchResultsFile.getCanonicalFile()+"\n\n"+
				
				"-a Add query subjects to registry? "+ addQuerySubjectsToRegistry+"\n"+
				"-c Add coreIds to registry subjects? "+ addCoreIdsToRegistry+"\n"+
				"-u Updated registry file "+ updatedRegistry+"\n"+
				"-s Max edit score for match "+ maxEditScoreForMatch+"\n"+
				"-p First missing key score penalty "+ missingOneKeyPenalty+ "\n"+
				"-k Subsequent missing key score penalty "+ missingAdditionalKeyPenalty+ "\n"+
				"-t Number threads "+ numberThreads+ "\n"+
				"-m Number of top matches to return "+ numberTopMatchesToReturn;
		
		IO.pl(opt);
	}




	/**Splits an object[] into chunks containing the minNumEach. Any remainder is evenly distributed over the prior.
	 * Note this is by reference, the array is not copied. */
	public static Subject[][] chunk (Subject[] s, int minNumEach){
		//watch out for cases where the min can't be met
		int numChunks = s.length/minNumEach;
		if (numChunks == 0) return new Subject[][]{s};
		
		double numLeftOver = (double)s.length % (double)minNumEach;
		
		int[] numInEach = new int[numChunks];
		for (int i=0; i< numChunks; i++) numInEach[i] = minNumEach;
		
		while (numLeftOver > 0){
			for (int i=0; i< numChunks; i++) {
				numInEach[i]++;
				numLeftOver--;
				if (numLeftOver == 0) break;
			}
		}
		//build chunk array
		Subject[][] chunks = new Subject[numChunks][];
		int index = 0;
		//for each chunk
		for (int i=0; i< numChunks; i++){
			//create container and fill it
			Subject[] sub = new Subject[numInEach[i]];
			for (int j=0; j< sub.length; j++) sub[j] = s[index++];
			chunks[i] = sub;
		}
		return chunks;
	}


	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                            Subject Id Match Maker : April 2022                   **\n" +
				"**************************************************************************************\n" +
				"Attempts to match each subject's PHI keys (FirstLastName, DoB, Gender, MRN) against a\n"+
				"registry of the same.  Uses a sum of each key's LevenshteinEditDistance/Length as the\n"+
				"distance metric with penalties for missing keys. Outputs each query, the top hits,\n"+
				"and an indication if the hit passes the max score cutoff.  If indicated, will asign\n"+
				"coreIds to queries not matched and write out an updated registry.\n"+
				
				"\nRequired:\n"+
				"-r File containing a registry of subjects, tab delimited txt file(.gz/.zip OK), one\n"+
				"      subject per line: lastName firstName dobMonth(1-12) dobDay(1-31)\n"+
				"      dobYear(1900-2050) gender(M|F) mrn coreIds HCIPatientIds. The last two columns\n"+
				"      are optional. Comma delimit multiple ids. Use space for missing info. Example:\n"+
				"      Biden Joseph 11 20 1942 M 19485763 vd3ec3XR,xp2kC2Gb 7474732\n"+ 
				"-q File containing queries to match to the registry, ditto.\n"+
				"-o File to write out the match results\n"+
				
				"\nOptional:\n"+
				"-a Add query subjects that failed to match, to the registry and assign them a coreId,\n"+
				"      requires -u\n"+
				"-c Add coreIds to registry subjects that lack them, requires -u\n"+
				"-u File to write the updated registry, requires -a and or -c\n"+
				"-s Max edit score for match, defaults to 0.12\n"+
				"-p Score penalty for a single missing key, defaults to 0.1\n"+
				"-k Score penatly for additional missing keys, defaults to 1\n"+
				"-t Number of threads to use, defaults to all\n"+
				"-m Number of top matches to return per query, defaults to 3\n"+
				
				"\nExample: java -jar pathToUSeq/Apps/SubjectIdMatchMaker -r transGen.simm.txt.zip \n"+
				"\n**************************************************************************************\n");
	}

	public Subject[] getTestSubjects() {
		return testSubjects;
	}
	public int getNumberTopMatchesToReturn() {
		return numberTopMatchesToReturn;
	}
	public double getMissingOneKeyPenalty() {
		return missingOneKeyPenalty;
	}
	public CoreId getCoreId() {
		return coreId;
	}
	public double getMissingAdditionalKeyPenalty() {
		return missingAdditionalKeyPenalty;
	}

}
