package trans.roc;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import util.bio.annotation.*;
import util.bio.parsers.gff.*;
import util.gen.*;

/**
 * Finds RocWindows upstream and withing the neighborhood of the ATG start site of a patterned gene list.
 */
public class ParsePatternedWindows {
	
	//fields
	private File resultsDirectory;
	private File gffFile;
	private File rocWindowDirectory;
	private File patternedGeneListFile;
	private File nonPatternedGeneListFile;
	private GeneGroup[] patternedGeneGrps;
	private GeneGroup[] nonPatternedGeneGrps;
	private int sizeUpstream = 2500;
	private double posCutOff = 0;
	private double negCutOff = 0;
	
	public ParsePatternedWindows(String[] args){
		processArgs(args);
		System.out.println("\nLaunching...");
		
		//process gff file retrieving an array of GeneGroup
		System.out.println("\tProcessing GFF file...");
		DmelRel4Extractor ex = new DmelRel4Extractor();
		ex.setChromosomeNameAppender("chr");
		ex.extract(gffFile, true);
		ArrayList geneGrpsAL = ex.getGeneGroupArrayList();
		int num = geneGrpsAL.size();
		//convert to a hashMap
		HashMap map = new HashMap(num);
		GeneGroup g;
		for (int i=0; i<num; i++){
			g = (GeneGroup)geneGrpsAL.get(i);
			map.put(g.getName(), g);
		}
		//make positive list
		patternedGeneGrps = extractGeneGroups(map, IO.loadFileIntoStringArray(patternedGeneListFile));
		System.out.println("\tNumber of Patterned Genes: "+patternedGeneGrps.length);
		
		//make negative list
		nonPatternedGeneGrps = extractGeneGroups(map, IO.loadFileIntoStringArray(nonPatternedGeneListFile));
		System.out.println("\tNumber of Non Patterned Genes: "+nonPatternedGeneGrps.length);

		//process RocWindow[] array
		RocWindow[] allRocWindows = loadRocWindows(IO.extractFiles(rocWindowDirectory));
		System.out.println("\tNumber of RocWindows: "+allRocWindows.length);
		
		//filter for a cutoff
		RocWindow[] pwmPassWindows = filterRocWindows(allRocWindows, posCutOff, true);
		System.out.println("\tNumber that Pass Pwm Positive Cutoff: "+pwmPassWindows.length);
		RocWindow[] posWindows = fetchProximalWindows(patternedGeneGrps, pwmPassWindows, sizeUpstream);
		System.out.println("\tNumber that also pass the proximity filter: "+posWindows.length);
		
		//negative windows
		RocWindow[] pwmFailWindows = filterRocWindows(allRocWindows, negCutOff, false);
		System.out.println("\n\tNumber that pass pwm negative cutoff "+ pwmFailWindows.length);
		
		//find those near a nonPatterned gene
		RocWindow[] negWinNearNonPat = fetchProximalWindows(nonPatternedGeneGrps, pwmFailWindows, sizeUpstream);
		System.out.println("\tNumber near a non patterned gene "+ negWinNearNonPat.length);
		
		//find those near a Patterned gene
		RocWindow[] negWinNearPat = fetchProximalWindows(patternedGeneGrps, pwmFailWindows, sizeUpstream *2);
		System.out.println("\tNumber near a patterned gene "+ negWinNearPat.length);

		//subtract those near pat
		//make hashset of negPosWindows and find difference
		HashSet negPat = Misc.loadHashSet(negWinNearPat);
		ArrayList realNeg = new ArrayList();
		for (int i=0; i< negWinNearNonPat.length; i++){
			if (negPat.contains(negWinNearNonPat[i]) == false) realNeg.add(negWinNearNonPat[i]);
		}
		RocWindow[] negWindows = new RocWindow[realNeg.size()];
		realNeg.toArray(negWindows);
		System.out.println("\tFinal number of pwm fail windows near a non patterned gene but not near a patterned gene "+negWindows.length);
		printChroms(negWindows);
		
		//save results to file
		//sgr
		File posSGR = new File(resultsDirectory,"posCenterScores.sgr");
		File negSGR = new File(resultsDirectory,"negCenterScores.sgr");
		writeSGRFile(posWindows, posSGR);
		writeSGRFile(negWindows, negSGR);
		//txt
		File posTxt = new File(resultsDirectory,"posWindows.roc");
		File negTxt = new File(resultsDirectory,"negWindows.roc");
		File allTxt = new File(resultsDirectory,"allWindows.roc");
		writeTxtRocFile(posWindows, posTxt);
		writeTxtRocFile(negWindows, negTxt);
		writeTxtRocFile(allRocWindows, allTxt);
		
		

		
		
		
	}
	
	/**Writes an sgr file.*/
	public static void writeSGRFile(RocWindow[] r, File file){
		try{
			PrintWriter out = new PrintWriter(new FileWriter(file));
				for (int i=0; i< r.length; i++){
					out.println(r[i].getChromosome()+"\t"+r[i].getMiddle()+"\t"+r[i].getScore());
				}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	/**Writes a chrom, start, stop, score file.*/
	public static void writeTxtRocFile(RocWindow[] r, File file){
		try{
			PrintWriter out = new PrintWriter(new FileWriter(file));
				for (int i=0; i< r.length; i++){
					out.println(r[i].getChromosome()+"\t"+r[i].getStart()+"\t"+r[i].getEnd()+"\t"+r[i].getMiddle()+"\t"+r[i].getScore());
				}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	/**reads a txt roc window (chrom, start, stop, middle(not used), score).*/
	public static RocWindow[] loadTxtRocFile(File file){
		try{
			BufferedReader in = new BufferedReader(new FileReader(file));
			String line;
			String[] tokens;
			ArrayList al = new ArrayList();
			while ((line = in.readLine()) !=null){
				tokens = line.split("\\s+");
				if (tokens.length == 5){
					al.add(new RocWindow(tokens[0], 
							Integer.parseInt(tokens[1]), 
							Integer.parseInt(tokens[2]), 
							Double.parseDouble(tokens[4])     
							));
				}
			}
			RocWindow[] win = new RocWindow[al.size()];
			al.toArray(win);
			return win;
		} catch (Exception e){
			e.printStackTrace();
		}
		return null;
	}
	
	
	public static void printChroms(RocWindow[] r){
		HashSet h = new HashSet();
		int num = r.length;
		for (int i=0; i<num; i++){
			h.add(r[i].getChromosome());
		}
		System.out.println("Chroms: "+h);
	}
	
	/**Scans an array of RocWindow for those that are greater than (or less than if boolean is false) or
	 * equal to the cutOff. Returns those that pass.*/
	public static RocWindow[] filterRocWindows(RocWindow[] windows, double cutOff, boolean greaterThan){
		int num = windows.length;
		ArrayList al = new ArrayList();
		if (greaterThan){
			for (int i=0; i<num; i++){
				if (windows[i].getScore()>= cutOff) al.add(windows[i]);
			}
		}
		else {
			for (int i=0; i<num; i++){
				if (windows[i].getScore()<= cutOff) al.add(windows[i]);
			}
		}
		RocWindow[] r = new RocWindow[al.size()];
		al.toArray(r);
		return r;
	}
	
	/**Loads all the RocWindows into one from serialized ArrayLists of RocWindows.*/
	public static RocWindow[] loadRocWindows(File[] arrayLists){
		RocWindow[] rocs = null;
		ArrayList all = (ArrayList)IO.fetchObject(arrayLists[0]);
		for (int i=1; i< arrayLists.length; i++){
			ArrayList toAdd = (ArrayList)IO.fetchObject(arrayLists[i]);
			all.addAll(  toAdd);
		}
		rocs = new RocWindow[all.size()];
		all.toArray(rocs);
		return rocs;
	}
	
	/**Pulls out a sub set of GeneGroups with the names in the geneList.*/
	public static GeneGroup[] extractGeneGroups(HashMap geneGroups, String[] geneList){
		//pull out those in hotGeneListFile
		int numHot = geneList.length;
		ArrayList hotGeneGrps = new ArrayList(numHot);
		for (int i=0; i<numHot; i++){
			if (geneGroups.containsKey(geneList[i])) hotGeneGrps.add(geneGroups.get(geneList[i]));
			else System.out.println("\tWarning, couldn't find '"+geneList[i]+"' in gff derived gene groups.");
		}
		GeneGroup[] geneGrps = new GeneGroup[hotGeneGrps.size()];
		hotGeneGrps.toArray(geneGrps);
		return geneGrps;
	}
	
	public static RocWindow[] fetchProximalWindows(GeneGroup[] geneGrps, RocWindow[] windows, int sizeUpstream){
		HashSet hits = new HashSet();
		for (int i=0; i< geneGrps.length; i++){
			//for each gene calc distance to win
			int ATG = geneGrps[i].getGeneRep().getStartATGPosition();
			for (int j=0; j<windows.length; j++){
				//same chromosome?
				if (geneGrps[i].getChromosome().equals(windows[j].getChromosome())){
					int gap;
					if (geneGrps[i].getOrientation() == 1) gap = ATG-windows[j].getMiddle();
					else gap = windows[j].getMiddle() - ATG;
					if (gap >= -350 && gap <= sizeUpstream) hits.add(windows[j]);
				}
			}
		}
		//return array
		RocWindow[] pos = new RocWindow[hits.size()];
		hits.toArray(pos);
		Arrays.sort(pos);
		return pos;

	}
	

	

	
	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**              Parse Patterned Windows:     June 2005                               **\n" +
				"**************************************************************************************\n" +
				"Finds RocWindows upstream and withing the neighborhood of the ATG start site of a \n" +
				"patterned gene list.\n\n" +
				
				"-g Full path file text for the DmelRel4.1 Stripped GFF3 file.\n" +
				"-p Full path file text for the RocWindow[] directory.\n" +
				"-c Full path file text for the patterned CG names file\n"+
				"-n Full path file text for the non patterned CG names file\n"+
				"-r Full path file text for the directory in which to write results files.\n"+
				"-b Size of neighborhood in bp, default is 2000 \n"+
				"-x Positive cutoff \n"+
				"-y Negative cutoff \n\n"+
				
		"**************************************************************************************\n");
	}
	
	/**This method will process each argument and assign new varibles*///stripGFF
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'g': gffFile = new File(args[i+1]); i++; break;
					case 'r': resultsDirectory = new File(args[i+1]); i++; break;
					case 'p': rocWindowDirectory = new File(args[i+1]); i++; break;
					case 'c': patternedGeneListFile = new File(args[i+1]); i++; break;
					case 'n': nonPatternedGeneListFile = new File(args[i+1]); i++; break;
					case 'b': sizeUpstream =Integer.parseInt(args[i+1]); i++; break;
					case 'x': posCutOff = Double.parseDouble(args[i+1]); i++; break;
					case 'y': negCutOff = Double.parseDouble(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter request: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
		if (gffFile == null || rocWindowDirectory == null){ //|| dgcFile == null){
			System.out.println("\nPlease enter the required files!\n");
			System.exit(0);
		}
		
	}
	
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ParsePatternedWindows (args);
	}
}
