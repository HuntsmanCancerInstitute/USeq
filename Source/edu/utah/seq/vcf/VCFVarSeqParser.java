package edu.utah.seq.vcf;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import util.gen.*;

/**Splits variants into those that intersect and those that do not.
 * @author Nix
 * */
public class VCFVarSeqParser {

	//user fields
	private HashMap<String, File> originalVcfNameFile = null;
	private File originalDir = null;
	private File forVarseqImportDir = null;
	private File fromVarseqExportDir = null;
	private File finalFilteredDir = null;
	private String[] filtersToExclude = null;
	private boolean dropSamples = false;

	//constructor
	public VCFVarSeqParser(String[] args){
		try {
			//start clock
			long startTime = System.currentTimeMillis();

			//process args
			processArgs(args);

			//process dirs
			processDirs();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
			System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing VCF files.");
		}
	}

	private void processDirs() throws FileNotFoundException, IOException {
		
		//original dir is needed for everything
		if (originalDir == null || originalDir.exists()==false ) Misc.printErrAndExit("Cannot find your -o originial dir containing vcf files to parse "+originalDir);

		//figure out what they want to do
		if (forVarseqImportDir != null) {
			//pull original vcf files
			File[] ori = fetchVcfFiles(originalDir);
			originalVcfNameFile = new HashMap<String, File>();
			for (File f: ori) originalVcfNameFile.put(Misc.removeExtension(f.getName()), f);
			forVarseqImportDir.mkdirs();
			IO.pl("Converting the originals into vcfs for import into VarSeq, saving files to: "+forVarseqImportDir);
			convertOriginalsForVarSeq();
		}
		
		else if (fromVarseqExportDir != null && fromVarseqExportDir.exists() && finalFilteredDir!= null) {
			finalFilteredDir.mkdirs();
			IO.pl("Parsing VarSeq export vcfs to pull maching original vcf records, saving files to: "+finalFilteredDir);
			parseOriRecordsMatchingVarSeq();
		}
		
		else {
			Misc.printErrAndExit("Looks like there is a problem with the dirs. Provide complete option sets \n\t1) -o and -i "
					+ "(to generate vcf files for VarSeq import)  OR \n\t2) -o -e -f (to pull ori vcf records that match the "
					+ "filtered VarSeq files)." );
		}


	}

	private void parseOriRecordsMatchingVarSeq() throws IOException {
		File[] vsFiltered = fetchVcfFiles(fromVarseqExportDir);
		
		//load hash with VRIDs
		HashMap<String, HashSet<Integer>> vridsToKeep = loadVRIDs(vsFiltered);
		
		//for each file name
		for (String fileName: vridsToKeep.keySet()) {
			//fetch the original
			File originalVcf = new File(originalDir, fileName+".vcf.gz");
			if (originalVcf.exists() == false) originalVcf = new File(originalDir, fileName+".vcf");
			if (originalVcf.exists() == false) throw new IOException ("\nFailed to find the original vcf for "+fileName +" in "+originalDir);
			
			//walk it printing only those line numbers in the hash
			File finalVcf = new File (finalFilteredDir, fileName+".VSFilt.vcf.gz");
			Gzipper out = new Gzipper(finalVcf);
			BufferedReader in = IO.fetchBufferedReader(originalVcf);
			HashSet<Integer> lineIndexesToKeep = vridsToKeep.get(fileName);
			
			String line;
			int lineNumber = 0;
			int numFound = 0;
			int numTotal = 0;
			
			//for each line in the file
			while ((line = in.readLine()) != null){
				line = line.trim();
				//header? just print out
				if (line.startsWith("#") == false) {
					numTotal++;
					if (lineIndexesToKeep.contains(new Integer(lineNumber))) {
						out.println(line);
						numFound++;
					}

				}
				else out.println(line);
				lineNumber++;
			}

			//close the IO
			in.close();
			out.close();
			
			if (numFound != lineIndexesToKeep.size()) throw new IOException("The # VRIDs found "+numFound+" doesn't match the number exported "+lineIndexesToKeep.size()+" for "+fileName);

			IO.pl("\t"+numTotal+"\t->\t"+numFound+"\t"+fileName);


		}

	}

	private HashMap<String, HashSet<Integer>> loadVRIDs(File[] vsFiltered) throws IOException {
		HashMap<String, HashSet<Integer>> vridsToKeep = new HashMap<String, HashSet<Integer>>();
		
		//Pattern pat = Pattern.compile(".*VRID=([\\w,\\.]+);.+");
		Pattern pat = Pattern.compile(".*VRID=([\\w,\\.]+).*");
				
		for (File f: vsFiltered) {
			
			BufferedReader in = IO.fetchBufferedReader(f);
			String line;
			while ((line = in.readLine()) != null){
				line = line.trim();
				if (line.startsWith("#") == false) {
					//#CHROM POS ID REF ALT QUAL FILTER INFO ......
					//   0    1   2  3   4   5     6     7
					String[] tokens = Misc.TAB.split(line);
					
					//fetch the vrids 1102562_Anno_Hg38.anno_29,1103990_Anno_Hg38.anno_48
					Matcher mat= pat.matcher(tokens[7]);
					if (mat.matches()) {
						String[] vrids = Misc.COMMA.split(mat.group(1));
						for (String s: vrids) {
							int lastUnderscoreIndex = s.lastIndexOf("_");
							String name = s.substring(0, lastUnderscoreIndex);
							String number = s.substring(lastUnderscoreIndex+1);
							//IO.pl(name +" "+number);
							HashSet<Integer> ints = vridsToKeep.get(name);
							if (ints == null) {
								ints = new HashSet<Integer>();
								vridsToKeep.put(name, ints);
							}
							ints.add(Integer.parseInt(number));
						}
					}
					else throw new IOException ("Failed to find any VRIDs in '"+line+"' from "+f);
					
				}
			}
			//close the IO
			in.close();

		}
		return vridsToKeep;
	}

	
	private void convertOriginalsForVarSeq() throws FileNotFoundException, IOException {
		//for each ori file
		for (String oriName: originalVcfNameFile.keySet()) addID(oriName);


	}

	/*Appends VRID=xxx onto each INFO cell.*/
	private void addID(String oriName) throws FileNotFoundException, IOException {
		File forVarSeq = new File (forVarseqImportDir, oriName+".vvsp.vcf.gz");
		Gzipper out = new Gzipper(forVarSeq);
		BufferedReader in = IO.fetchBufferedReader(originalVcfNameFile.get(oriName));

		int recordIndex = 0;
		int numExcluded = 0;
		int numSaved = 0;
		String line;
		boolean addInfo = true;
		//for each line in the file
		while ((line = in.readLine()) != null){
			line = line.trim();
			//header? just print out
			if (line.startsWith("#")) {
				if (addInfo && line.startsWith("##INFO=")) {
					addInfo = false;
					out.println("##INFO=<ID=VRID,Number=.,Type=String,Description=\"Variant Record ID for VarSeq export.\">");
				}
				if (dropSamples && line.startsWith("#CHROM")) out.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
				else out.println(line);
			}
			//data line
			else {	
				//#CHROM POS ID REF ALT QUAL FILTER INFO ......
				//   0    1   2  3   4   5     6     7
				String[] tokens = Misc.TAB.split(line);
				//filter it?
				boolean save = true;
				if (filtersToExclude !=null) {
					for (String f: filtersToExclude) {
						if (tokens[6].contains(f)) {
							numExcluded++;
							save = false;
							break;
						}
					}
				}
				
				if (save) {
					tokens[7] = tokens[7]+";VRID="+oriName+"_"+recordIndex;
					if (dropSamples) out.println(tokens[0]+"\t"+ tokens[1]+"\t"+ tokens[2]+"\t"+ tokens[3]+"\t"+ tokens[4]+"\t"+ tokens[5]+"\t"+ tokens[6]+"\t"+ tokens[7]);
					else out.println(Misc.stringArrayToString(tokens, "\t"));
					numSaved++;
				}
			}
			recordIndex++;
		}

		//close the IO
		in.close();
		out.close();

		if (filtersToExclude ==null) IO.pl("\t"+recordIndex+"\tOri records\t"+originalVcfNameFile.get(oriName).getName());
		else IO.pl("\t"+numSaved+" saved\t"+numExcluded+" filtered\t"+originalVcfNameFile.get(oriName).getName());

	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFVarSeqParser(args);
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
					case 'o': originalDir = new File(args[++i]); break;
					case 'i': forVarseqImportDir = new File(args[++i]); break;
					case 'e': fromVarseqExportDir = new File(args[++i]); break;
					case 'f': finalFilteredDir = new File(args[++i]); break;
					case 'x': filtersToExclude = Misc.COMMA.split(args[++i]); break;
					case 'd': dropSamples = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
	}	

	public static File[] fetchVcfFiles(File dir) {
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(dir,".vcf");
		tot[1] = IO.extractFiles(dir,".vcf.gz");
		tot[2] = IO.extractFiles(dir,".vcf.zip");
		return IO.collapseFileArray(tot);
	}



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              VCF VarSeq Parser: May 2020                         **\n" +
				"**************************************************************************************\n" +
				"VVSP either creates vcf files for import into VarSeq by inserting a unique ID into the\n"+
				"INFO field of each record OR, VVSP takes vcf file(s) exported from VarSeq and pulls the\n"+
				"original unmodified vcf record via the unique ID. Use this tool when you want to\n"+
				"just filter variants using VarSeq and retain the unmodified original vcf records.\n"+
				"If you aren't going to filter on sample info, e.g. with t/n somatic variants, use -d\n"+
				"Be carefult to not modify the original vcfs prior to rematching with the VarSeq export\n"+
				"vcf(s)!\n"+
				
				"\nWhen exporting vcf files from VarSeq, include the VRID INFO tag AND where possible\n"+
				"only export Affected Samples, e.g. Field Selection-> Variant Info-> check VRID. On the\n"+
				"Options page select the Affected Samples bubble.\n"+

				"\nParams (either complete -o -i -x OR -o -e -f):\n"+
				"-o Directory containing original vcfs to modify for VarSeq import\n"+
				"     (xxx.vcf(.gz/.zip OK))\n"+
				"-i Directory to save the modified vcf files for VarSeq import\n"+
				"-e Directory containing the exported VarSeq vcf file(s)\n"+
				"-f Directory to save the final filtered original vcf records\n"+
				"-x List of FILTER keys to exclude from the original vcfs, comma delimited, no spaces,\n"+
				"     case sensitive\n"+
				"-d Drop sample info for VarSeq vcf file import.\n"+

				"\nExamples:\n"+
				"   1) java -Xmx1G -jar pathTo/USeq/Apps/VCFVarSeqParser -o StartingOriVcfs/ \n" +
				"       -i VcfsForVarSeqImport -x CallFreq,VCFBkz -d \n"+
				"   2) java -Xmx1G -jar pathTo/USeq/Apps/VCFVarSeqParser -o StartingOriVcfs/ \n" +
				"       -e VarSeqExportedVcfs/ -f FinalFilteredOriVcfs/\n"+

				"\n**************************************************************************************\n");

	}
}
