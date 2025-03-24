package edu.utah.seq.vcf;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import edu.utah.seq.vcf.fdr.VCFFdrEstimator;
import util.gen.*;

/**Takes an OncoKB annotated MAF file and appends the OncoKB column values to the original VCF as INFO fields
 * @author Nix
 * */
public class OncoKB2VCF {

	//user fields
	private File vcfFile = null;
	private File mafFile = null;
	private File okbVcfFile = null;
	private String[] mafFieldsToAdd = null;

	//vcf
	private TreeSet<String> vcfHeader = new TreeSet<String>();
	private String vcfChromLine = null; // #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR
	private String vcfFileFormat = null; // ##fileformat=VCFv4.1
	private BufferedReader vcfIn = null;
	private HashSet<String> vcfIds = new HashSet<String>();
	private File headerOkbVcfFile;
	private File bodyOkbVcfFile;

	//maf
	private LinkedHashMap<String, String> knownOkbNameDesc = new LinkedHashMap<String, String>();
	private LinkedHashMap<String, Integer> okbNameIndexes = new LinkedHashMap<String, Integer>();
	private LinkedHashSet<String> foundOkbNames = new LinkedHashSet<String>();
	private HashMap<String, String[]> mafIdLine = new HashMap<String, String[]>();
	private BufferedReader mafIn = null;
	private int vcfIdIndex = -1;
	private int annotatedIndex = -1;


	//constructor
	public OncoKB2VCF(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);

		//load the vcf and maf headers
		loadVcfHeader();
		loadMafHeader();
		
		joinRecords();
		writeOutVcfHeader();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	private void writeOutVcfHeader() {
		IO.pl("Writing out vcf header and final vcf...");
		Gzipper out = null;

		try {
			out = new Gzipper(headerOkbVcfFile);
			out.println(vcfFileFormat);
			for (String s: foundOkbNames) vcfHeader.add(this.knownOkbNameDesc.get(s));
			for (String s: vcfHeader) out.println(s);
			out.println(vcfChromLine);
			out.close();
			
			//concatinate gz files
			ArrayList<File> toConcat = new ArrayList<File>();
			toConcat.add(headerOkbVcfFile);
			toConcat.add(bodyOkbVcfFile);
			IO.concatinateFiles(toConcat , okbVcfFile);
			IO.pl("\t"+okbVcfFile);
			
		} catch (Exception e) {
			e.printStackTrace();
			headerOkbVcfFile.delete();
			Misc.printErrAndExit("\nProblem writing new vcf header ");
		} finally {
			out.closeNoException();
		}

		
	}

	private void joinRecords () {
		IO.pl("Adding OncoKB annotations to vcf records...");
		String vcfLine = null;
		String mafLine = null;
		String[] vcfFields = null;
		String[] mafFields = null;
		Gzipper out = null;
		int numMissing = 0;

		try {
			out = new Gzipper(bodyOkbVcfFile);
			//for each vcfLine
			while ((vcfLine = vcfIn.readLine())!= null) {
				
				//#CHROM-0 POS-1 ID-2 REF-3 ALT-4 QUAL-5 FILTER-6 INFO-7	FORMAT	NORMAL	TUMOR
				vcfFields = Misc.TAB.split(vcfLine);
				
				//check IDs
				String vcfId = vcfFields[2];
				//not needed?
				if (vcfIds.contains(vcfId)) throw new Exception("\nERROR, seeing duplicate VCF ID, these must be unique. "+vcfId);
				vcfIds.add(vcfId);
				
				mafFields = mafIdLine.get(vcfId);
				
				if (mafFields == null) {
					IO.el("\tWARNING, failed to find the Vcf: "+vcfId+" in the maf? Skipping.");
					numMissing++;
					if (numMissing > 9) throw new Exception("\nToo many missing vcfIds in maf, aborting!");
					out.println(vcfLine);
				}
				else {
					//was it annotated
					String annotated = mafFields[annotatedIndex];
					if (annotated.equals("FALSE")) out.println(vcfLine);
					else out.println(mergeInfo(vcfFields, mafFields));
				}
			}
			if (vcfChromLine == null || vcfFileFormat == null) throw new Exception("Failed to parse one or both "+vcfFileFormat+" "+vcfChromLine);
			
		} catch (Exception e) {
			e.printStackTrace();
			bodyOkbVcfFile.delete();
			Misc.printErrAndExit("\nJoining records.");
		} finally {
			IO.closeNoException(vcfIn);
			out.closeNoException();
		}
	}

	private String mergeInfo(String[] vcfFields, String[] mafFields) {
		ArrayList<String> al = new ArrayList<String>();
		// look for requested annotations
		for (String s: mafFieldsToAdd) {
			int index = okbNameIndexes.get(s);
			if (mafFields.length > index && mafFields[index].length()!=0) {
				String cleaned = clean(mafFields[index]);
				if (cleaned.length()!=0) {
					al.add(s+"="+cleaned);
					foundOkbNames.add(s);
				}
			}
		}
		if (al.size()!=0) {
			vcfFields[7] = vcfFields[7]+";"+Misc.stringArrayListToString(al, ";");
		}
		return Misc.stringArrayToString(vcfFields, "\t");
	}

	//no whitespace, semicolons, or equals-signs
	private static Pattern toDrop = Pattern.compile("[\\s;=]+");
	
	private static String clean(String toClean) {
		if (toClean.contains("Unknown")) return "";
		toClean = toClean.trim();
		return toDrop.matcher(toClean).replaceAll("_");
	}

	private void loadMafHeader () {
		IO.pl("Loading the maf header...");
		String line = null;
		
		try {
			mafIn = IO.fetchBufferedReader(mafFile);
			//load the header
			while ((line = mafIn.readLine())!= null) {
				if (line.startsWith("Hugo_Symbol")) {
					String[] f = Misc.TAB.split(line);
					for (int i=0; i< f.length; i++) {
						if (f[i].equals("vcf_id")) vcfIdIndex = i;
						else if (f[i].equals("ANNOTATED")) annotatedIndex = i;
						else if (knownOkbNameDesc.get(f[i]) !=null) okbNameIndexes.put(f[i], i);
					}
					break;
				}
			}
			if (vcfIdIndex == -1 || annotatedIndex == -1) throw new Exception("Failed to parse the 'vcf_id' or 'ANNOTATED' from "+line);
			
			//load the data
			while ((line = mafIn.readLine())!= null) {
				String[] mafFields = Misc.TAB.split(line);
				String mafId = mafFields[vcfIdIndex];
				mafIdLine.put(mafId, mafFields);
			}
			
			mafIn.close();
		} catch (Exception e) {
			IO.closeNoException(mafIn);
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing maf header for "+mafFile);
		} 
	}

	private void loadVcfHeader () {
		IO.pl("Loading the vcf header...");
		String line = null;

		try {
			vcfIn = IO.fetchBufferedReader(vcfFile);
			while ((line = vcfIn.readLine())!= null) {
				if (line.startsWith("#")) {

					if (vcfFileFormat == null && line.contains("##fileformat")) vcfFileFormat = line;
					else if (line.startsWith("#CHROM")) {
						vcfChromLine = line;
						break;
					}
					else vcfHeader.add(line);
				}
			}
			if (vcfChromLine == null || vcfFileFormat == null) throw new Exception("Failed to parse one or both "+vcfFileFormat+" "+vcfChromLine);
		} catch (Exception e) {
			IO.closeNoException(vcfIn);
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing vcf header for "+vcfFile);
		} 
	}


	/*From https://github.com/oncokb/oncokb-annotator, last reviewed 1 Jan 2025, check annually!*/
	private static String[] info = new String[]{
			"ANNOTATED", "##INFO=<ID=ANNOTATED,Number=1,Type=String,Description='OncoKB: TRUE, FALSE - Whether the variant is annotated by OncoKB successfully.'>",
			"ONCOKB_HUGO_SYMBOL", "##INFO=<ID=ONCOKB_HUGO_SYMBOL,Number=1,Type=String,Description='OncoKB: When annotating genomic change, we obtained gene hugo symbol from GenomeNexus. This can be cross-referenced with your own gene name.'>",
			"ONCOKB_PROTEIN_CHANGE", "##INFO=<ID=ONCOKB_PROTEIN_CHANGE,Number=1,Type=String,Description='OncoKB: When annotating genomic change, we obtained alteration protein change from GenomeNexus. This can be cross-referenced with your own protein change.'>",
			"ONCOKB_CONSEQUENCE", "##INFO=<ID=ONCOKB_CONSEQUENCE,Number=1,Type=String,Description='OncoKB: When annotating genomic change, we obtained alteration consequence from GenomeNexus. This can be cross-referenced with your own consequence/Variant Class.'>",
			"GENE_IN_ONCOKB", "##INFO=<ID=GENE_IN_ONCOKB,Number=1,Type=String,Description='OncoKB: TRUE, FALSE - Whether the gene has been curated by the OncoKB Team.'>",
			"VARIANT_IN_ONCOKB", "##INFO=<ID=VARIANT_IN_ONCOKB,Number=1,Type=String,Description='OncoKB: TRUE, FALSE - Whether the variant has been curated by the OncoKB Team. Note: when a variant does not exist, it may still have annotations.'>",
			"MUTATION_EFFECT", "##INFO=<ID=MUTATION_EFFECT,Number=1,Type=String,Description='OncoKB: Gain-of-function, Likely_Gain-of-function, Loss-of-function, Likely_Loss-of-function, Switch-of-function, Likely_Switch-of-function, Neutral, Likely_Neutral, Inconclusive, Unknown - The biological effect of a mutation/alteration on the protein function that gives rise to changes in the biological properties of cells expressing the mutant/altered protein compared to cells expressing the wildtype protein.'>",
			"MUTATION_EFFECT_CITATIONS", "##INFO=<ID=MUTATION_EFFECT_CITATIONS,Number=1,Type=String,Description='OncoKB: PMID, Abstract, Website link - All citations related to the biological effect.'>",
			"ONCOGENIC", "##INFO=<ID=ONCOGENIC,Number=1,Type=String,Description='OncoKB: Oncogenic, Likely_Oncogenic, Likely_Neutral, Inconclusive, Unknown, Resistance - In OncoKB, oncogenic is defined as referring to the ability to induce or cause cancer as described in the second edition of The Biology of Cancer by Robert Weinberg (2014).'>",
			"LEVEL_1", "##INFO=<ID=LEVEL_1,Number=1,Type=String,Description='OncoKB: FDA-recognized biomarker predictive of response to an FDA-approved drug in this indication'>",
			"LEVEL_2", "##INFO=<ID=LEVEL_2,Number=1,Type=String,Description='OncoKB: Standard care biomarker recommended by the NCCN or other professional guidelines predictive of response to an FDA-approved drug in this indication'>",
			"LEVEL_3A", "##INFO=<ID=LEVEL_3A,Number=1,Type=String,Description='OncoKB: Compelling clinical evidence supports the biomarker as being predictive of response to a drug in this indication'>",
			"LEVEL_3B", "##INFO=<ID=LEVEL_3B,Number=1,Type=String,Description='OncoKB: Standard care or investigational biomarker predictive of response to an FDA-approved or investigational drug in another indication'>",
			"LEVEL_4", "##INFO=<ID=LEVEL_4,Number=1,Type=String,Description='OncoKB: Compelling biological evidence supports the biomarker as being predictive of response to a drug'>",
			"LEVEL_R1", "##INFO=<ID=LEVEL_R1,Number=1,Type=String,Description='OncoKB: Standard care biomarker predictive of resistance to an FDA-approved drug in this indication'>",
			"LEVEL_R2", "##INFO=<ID=LEVEL_R2,Number=1,Type=String,Description='OncoKB: Compelling clinical evidence supports the biomarker as being predictive of resistance to a drug'>",
			"HIGHEST_LEVEL", "##INFO=<ID=HIGHEST_LEVEL,Number=1,Type=String,Description='OncoKB: LEVEL_1, LEVEL_2, LEVEL_3A, LEVEL_3B, LEVEL_4, LEVEL_R1, LEVEL_R2 -The highest level of evidence for therapeutic implications. Order: LEVEL_R1 > LEVEL_1 > LEVEL_2 > LEVEL_3A > LEVEL_3B > LEVEL_4 > LEVEL_R2'>",
			"HIGHEST_SENSITIVE_LEVEL", "##INFO=<ID=HIGHEST_SENSITIVE_LEVEL,Number=1,Type=String,Description='OncoKB: LEVEL_1, LEVEL_2, LEVEL_3A, LEVEL_3B, LEVEL_4 - The highest sensitive level of evidence for therapeutic implications. Order: LEVEL_1 > LEVEL_2 > LEVEL_3A > LEVEL_3B > LEVEL_4'>",
			"HIGHEST_RESISTANCE_LEVEL", "##INFO=<ID=HIGHEST_RESISTANCE_LEVEL,Number=1,Type=String,Description='OncoKB: LEVEL_R1, LEVEL_R2 - The highest resistance level of evidence for therapeutic implications. Order: LEVEL_R1 > LEVEL_R2'>",
			"TX_CITATIONS", "##INFO=<ID=TX_CITATIONS,Number=1,Type=String,Description='OncoKB: PMID, Abstract, Website link - All citations related to therapeutic implications.'>",
			"LEVEL_Dx1", "##INFO=<ID=LEVEL_Dx1,Number=1,Type=String,Description='OncoKB: Tumor type the level of evidence is assigned to - The leveled diagnostic implications.'>",
			"LEVEL_Dx2", "##INFO=<ID=LEVEL_Dx2,Number=1,Type=String,Description='OncoKB: Tumor type the level of evidence is assigned to - The leveled diagnostic implications.'>",
			"LEVEL_Dx3", "##INFO=<ID=LEVEL_Dx3,Number=1,Type=String,Description='OncoKB: Tumor type the level of evidence is assigned to - The leveled diagnostic implications.'>",
			"HIGHEST_DX_LEVEL", "##INFO=<ID=HIGHEST_DX_LEVEL,Number=1,Type=String,Description='OncoKB: LEVEL_Dx1, LEVEL_Dx2, LEVEL_Dx3 - The highest level of evidence for diagnostic implications.'>",
			"DX_CITATIONS", "##INFO=<ID=DX_CITATIONS,Number=1,Type=String,Description='OncoKB: PMID, Abstract, Website link - All citations related to diagnostic implications.'>",
			"LEVEL_Px1", "##INFO=<ID=LEVEL_Px1,Number=1,Type=String,Description='OncoKB: Tumor type the level of evidence is assigned to - The leveled prognostic implications.'>",
			"LEVEL_Px2", "##INFO=<ID=LEVEL_Px2,Number=1,Type=String,Description='OncoKB: Tumor type the level of evidence is assigned to - The leveled prognostic implications.'>",
			"LEVEL_Px3", "##INFO=<ID=LEVEL_Px3,Number=1,Type=String,Description='OncoKB: Tumor type the level of evidence is assigned to - The leveled prognostic implications.'>",
			"HIGHEST_PX_LEVEL", "##INFO=<ID=HIGHEST_PX_LEVEL,Number=1,Type=String,Description='OncoKB: LEVEL_Px1, LEVEL_Px2, LEVEL_Px3 - The highest level of evidence for prognostic implications.'>",
			"PX_CITATIONS", "##INFO=<ID=PX_CITATIONS,Number=1,Type=String,Description='OncoKB: PMID, Abstract, Website link - All citations related to prognostic implications.'>",
			"GENE_SUMMARY", "##INFO=<ID=GENE_SUMMARY,Number=1,Type=String,Description='OncoKB: Brief overview of the gene and its role in cancer'>",
			"VARIANT_SUMMARY", "##INFO=<ID=VARIANT_SUMMARY,Number=1,Type=String,Description='OncoKB:  Variant summary describes the variant oncogenicity, last review if it is VUS'>",
			"TUMOR_TYPE_SUMMARY", "##INFO=<ID=TUMOR_TYPE_SUMMARY,Number=1,Type=String,Description='OncoKB: Tumor type summary describes the therapeutic implication that applies to the indication'>",
			"DIAGNOSTIC_SUMMARY", "##INFO=<ID=DIAGNOSTIC_SUMMARY,Number=1,Type=String,Description='OncoKB: Diagnostic summary that applies to the indication, for hematologic malignancies only'>",
			"PROGNOSTIC_SUMMARY", "##INFO=<ID=PROGNOSTIC_SUMMARY,Number=1,Type=String,Description='OncoKB: Prognostic summary that applies to the indication, for hematologic malignancies only'>",
			"MUTATION_EFFECT_DESCRIPTION", "##INFO=<ID=MUTATION_EFFECT_DESCRIPTION,Number=1,Type=String,Description='OncoKB: The mutation effect description provides a brief overview of the biological and oncogenic effect of the VPS and includes appropriate references to peer-reviewed literature.'>"
	};


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new OncoKB2VCF(args);
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
					case 'm': mafFile = new File(args[++i]); break;
					case 'v': vcfFile = new File(args[++i]); break;
					case 'o': okbVcfFile = new File(args[++i]); break;
					case 'a': mafFieldsToAdd = Misc.COMMA.split(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//add the info fileds
		IO.pl("\nRecognized OncoKB annotations...");
		for (int i=0; i< info.length; i+=2) {
			knownOkbNameDesc.put(info[i], info[i+1]);
			IO.pl("\t"+info[i]+" -> "+info[i+1]);
		}
		
		if (mafFieldsToAdd == null) mafFieldsToAdd = new String[] {"MUTATION_EFFECT","ONCOGENIC","HIGHEST_LEVEL"};
		boolean failed = false;
		for (String s: mafFieldsToAdd) {
			if (knownOkbNameDesc.containsKey(s) == false) {
				failed = true;
				IO.pl("\tERROR: failed to find "+s);
			}
		}
		if (failed) Misc.printErrAndExit("Correct MAF column names for annotation.");
		
		if (mafFile == null || mafFile.exists()==false) Misc.printErrAndExit("\nERROR: failed to find your maf file! "+mafFile);
		if (vcfFile == null || vcfFile.exists()==false) Misc.printErrAndExit("\nERROR: failed to find your vcf file! "+vcfFile);
		String name = Misc.removeExtension(vcfFile.getName());
		if (okbVcfFile == null) okbVcfFile = new File (vcfFile.getParentFile(), name+".okb.vcf.gz");
		if (okbVcfFile.getName().endsWith(".gz") == false) Misc.printErrAndExit("\nERROR: your output vcf file must end in .gz "+okbVcfFile);
		headerOkbVcfFile = new File (vcfFile.getParentFile(), "tmpHeader_"+name+".okb.vcf.gz");
		bodyOkbVcfFile = new File (vcfFile.getParentFile(), "tmpBody_"+name+".okb.vcf.gz");
		headerOkbVcfFile.deleteOnExit();
		bodyOkbVcfFile.deleteOnExit();
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              OncoKB 2 VCF : March 2025                           **\n" +
				"**************************************************************************************\n" +
				"Converts OncoKB annotations from a maf file into vcf INFO fields and appends them to \n"+
				"the original vcf. The vcf ID column must contain unique values. Both files must be\n"+
				"sorted identically.\n"+

				"\nRequired Params:\n"+
				"-m MAF file annotated by the OncoKB MafAnnotator.py app (.gz/.zip OK)\n"+
				"-v Original VCF file that was converted to maf and then annotated.\n"+
				"-o (Optional) Output VCF file with OncoKB annotations, it should end in xxx.vcf.gz\n"+
				"-a (Optional) Comma delimited list, no spaces, of OncoKB annotations to transfer.\n"+
				"      Defaults to MUTATION_EFFECT,ONCOGENIC,HIGHEST_LEVEL\n"+

				"\nExample: java -Xmx1G -jar pathTo/USeq/Apps/OncoKB2VCF -v avatarP7.vcf.gz\n" +
				"       -m avatarP7.okb.maf.gz \n"+

				"\n**************************************************************************************\n");

	}
}
