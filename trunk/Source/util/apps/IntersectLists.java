package util.apps;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.*;
import util.gen.*;

/**
 * Compares two lists of words for intersection. Uses random permutation to assess fold enrichment and significance.
 */
public class IntersectLists {
	//fields
	private int numberPermutations = 1000;
	private int totalNumberOfSourceList;
	private File[] aFiles;
	private File[] bFiles;
	private File listAFile;
	private File listBFile;
	private boolean printLists;
	private HashSet<String> listA;
	private HashSet<String> listB;
	private int numListA;
	private int numListB;
	private int numIntersect;
	//permutation results
	private String coorPVal;
	private String antiCorrPVal;
	private String averageIntersection;
	private String foldEnrichment;

	public  IntersectLists(String[] args){
		//process the arguments
		processArguments(args);
		
		//make buffer to hold results
		StringBuffer sb = new StringBuffer();
		sb.append("#Intersect\tNameA\t#A\tFractionAIntersect\tNameB\t#B\tFractionBIntersect\tFoldEnrichment\tCorrPVal\tAntiCorrPVal\n");
		
		//for each set of files
		for (int x=0; x< aFiles.length; x++){
			listAFile = aFiles[x];
			for (int y=0; y<bFiles.length; y++){
				listBFile = bFiles[y];
				//load lists
				listA = IO.loadFileIntoHashSet(listAFile);
				listB = IO.loadFileIntoHashSet(listBFile);
				numListA = listA.size();
				numListB = listB.size();

				//run thru testNames and look for in one
				String[] testNames = Misc.hashSetToStringArray(listB);
				int num = testNames.length;
				ArrayList<String> found = new ArrayList<String>();
				for (int i=0; i<num; i++){
					if (listA.contains(testNames[i])) {
						found.add(testNames[i]);
						listA.remove(testNames[i]);
						listB.remove(testNames[i]);
					}
				}
				numIntersect = found.size();

				//random permutations

				permutate();

				//save results
				double fractIntOne = (double)numIntersect/(double) numListA;
				double fractIntTwo = (double)numIntersect/(double) numListB;
				sb.append(numIntersect); sb.append("\t");
				sb.append(Misc.removeExtension(listAFile.getName())); sb.append("\t");
				sb.append(numListA); sb.append("\t");
				sb.append(Num.formatNumber(fractIntOne, 4)); sb.append("\t");
				sb.append(Misc.removeExtension(listBFile.getName())); sb.append("\t");
				sb.append(numListB); sb.append("\t");
				sb.append(Num.formatNumber(fractIntTwo, 4)); sb.append("\t");
				sb.append(foldEnrichment+"("+numIntersect+"/"+averageIntersection+")"); sb.append("\t");
				sb.append(coorPVal); sb.append("\t");
				sb.append(antiCorrPVal); sb.append("\t\n");

				//print names?
				if (printLists){
					System.out.println("\n***************************************************\nCommon names in "+ listAFile.getName() + " and "+listBFile.getName()+":\n"+ Misc.stringArrayListToString(found, ", "));
					System.out.println("\n\nUncommon names in "+ listAFile.getName() + ":\n" + Misc.stringArrayToString((Misc.hashSetToStringArray(listA)), ", "));
					System.out.println("\n\nUncommon names in "+ listBFile.getName() + ":\n" + Misc.stringArrayToString((Misc.hashSetToStringArray(listB)), ", "));
				}
			}
		}
		System.out.println("\nResults:\n"+sb);
	}

	public void permutate(){
		//make master set and a setA
		String[] masterSet = new String[totalNumberOfSourceList];
		String[] listA = new String[numListA];
		for (int i=0; i< totalNumberOfSourceList; i++) {
			masterSet[i] = i+"";
			if (i<numListA) listA[i] = i+"";
		}
		//vars
		int totalInt = 0;
		int numMatchRealInt = 0;
		//for each permutation
		for (int i=0; i< numberPermutations; i++){
			//randomize master
			Misc.randomize(masterSet, System.currentTimeMillis());
			//draw a setB
			HashSet<String> listB = new HashSet<String>(numListB);
			for (int x =0; x< numListB; x++) listB.add(masterSet[x]);
			//intersect
			int numInt = 0;
			for (int x=0; x< listA.length; x++){
				if (listB.contains(listA[x])) numInt++;
			}			
			if (numInt >= numIntersect) numMatchRealInt++;
			totalInt += numInt;
		}

		//calculate summary stats
		if (numMatchRealInt == 0){
			coorPVal = "<1/"+numberPermutations;
			antiCorrPVal = "1";
		}
		else {
			double pval = (double)numMatchRealInt/(double)numberPermutations;
			coorPVal = Num.formatNumber(pval, 4);
			antiCorrPVal = Num.formatNumber((1-pval), 4);
		}

		averageIntersection = "";
		foldEnrichment = "";
		if (totalInt !=0) {
			double aveInt = (double)totalInt/(double)numberPermutations;
			averageIntersection = Num.formatNumber(aveInt, 4);
			double foldEnr = (double)numIntersect/aveInt;
			foldEnrichment = Num.formatNumber(foldEnr, 4);
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new IntersectLists(args);
	}	

	/**This method will process each argument and assign new varibles*/
	public void processArguments(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'a': aFiles = IO.extractFiles(new File(args[++i])); break;
					case 'b': bFiles = IO.extractFiles(new File(args[++i])); break;
					case 't': totalNumberOfSourceList = Integer.parseInt(args[++i]); break;
					case 'n': numberPermutations = Integer.parseInt(args[++i]); break;
					case 'p': printLists = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for files
		if (aFiles == null || aFiles[0].exists()== false || bFiles == null || bFiles[0].exists()== false){
			Misc.printExit("\nCannot find your list files?");
		}
		//look for totalNumberOfSource
		if (totalNumberOfSourceList == 0){
			Misc.printExit("\nEnter the total number of items in the master list from which your two lists were drawn.\n");
		}

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Intersect Lists: Dec 2008                             **\n" +
				"**************************************************************************************\n" +
				"IL intersects two lists (of genes) and using randomization, calculates the\n" +
				"significance of the intersection and the fold enrichment over random. Note, duplicate\n" +
				"items are filtered from each list prior to analysis.\n\n" +

				"-a Full path file text for list A (or directory containing), one item per line.\n" +
				"-b Full path file text for list B (or directory containing), one item per line.\n" +
				"-t The total number of unique items from which A and B were drawn.\n"+
				"-n Number of permutations, defaults to 1000.\n"+
				"-p Print the intersection sets (common, unique to A, unique to B) to screen.\n"+
				"\n"+

				"Example: java -Xmx1500M -jar pathTo/Apps/IntersectLists -a /Data/geneListA.txt -b \n" +
				"     /Data/geneListB.txt -t 28356 -n 10000\n\n" +

		"**************************************************************************************\n");

	}

}
