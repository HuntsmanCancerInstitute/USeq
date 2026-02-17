package edu.utah.seq.run.tempus.v3;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.json.JSONArray;
import org.json.JSONObject;

import edu.utah.hci.bioinfo.smm.Subject;
import edu.utah.hci.bioinfo.smm.SubjectMatchMaker;
import edu.utah.seq.pmr.PMRSearch;
import edu.utah.seq.vcf.json.tempusv3.TempusV3Json2Vcf;
import edu.utah.seq.vcf.json.tempusv3.TempusV3JsonCollection;
import edu.utah.seq.vcf.json.tempusv3.TempusV3JsonSummary;
import edu.utah.seq.vcf.json.tempusv3.TempusV3Order;
import edu.utah.seq.vcf.json.tempusv3.TempusV3Specimen;
import util.gen.IO;
import util.gen.Misc;

/**Works with v3.3+ Tempus json test results, otherwise they are ignored, use TempusDataWrangler for legacy analysis.* */
public class TempusPathoPrinter {

	//user defined fields
	private File[] jsonFiles = null;
	private File resultsDir = null;
	
	private HashMap<String, String> icd10CodeDesc = null;
	private HashMap<String, String> icdOCodeMorphology = null;
	private HashMap<String, String> icdOCodeTopology = null;
	private ArrayList<String> keys = new ArrayList<String>();

	public TempusPathoPrinter (String[] args) {
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);
			
			parseJsonFiles();
			
			saveKeyFile();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			IO.pl("\nDone! "+Math.round(diffTime)+" Min\n");
			
		} catch (Exception e) {
			IO.el("\nERROR running the TempusPathoPrinter, aborting");
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void saveKeyFile() {
		File f = new File(resultsDir, "testKeys.txt");
		IO.writeArrayList(keys, f);
	}

	private void parseJsonFiles() throws Exception {
		
		IO.pl("\nParsing json files...");
		for (File json: jsonFiles) {
			IO.pl("\t"+json.getName());
			TempusV3JsonSummary sum = TempusV3Json2Vcf.parseJsonNoVariants(json);
			JSONObject jo = new JSONObject();
			//add test id
			String rand = Misc.getRandomString(10);
			keys.add(sum.getTempusOrder().getAccessionId()+"\t"+rand);
			jo.put("test_order_id", rand);
			for (TempusV3Specimen s : sum.getTempusSpecimens()) {
				//is this a tumor sample
				String sc = s.getSampleCategory().toLowerCase();
				if (sc.contains("tumor")==false && sc.contains("heme")==false ) continue;
				//primary sample site
				jo.put("sample_site", s.getPrimarySampleSite());
				//path diagnosis
				String diag = s.getOriginPathLabDiagnosis();
				if (diag.trim().length()==0) {
					IO.el("\t"+json.getName()+" missing diag, skipping");
					continue;
				}
				//any w/o?
				diag = diag.replace(" c/w ", " consistent with ");
				// any w/
				diag = diag.replace(" w/ ", " with ");
				// and w/o
				diag = diag.replace(" w/o ", " with out ");
				jo.put("original_path_lab_diagnosis", diag);
				//these codes are not part of the pathology report, assigned by Tempus?
				//icd 10 diagnosis codes
				if (s.getTempusIcd10Code() != null) {
					for (String code: Misc.COMMA.split(s.getTempusIcd10Code())) {
						JSONObject o = new JSONObject();
						o.put("type", "ICD-10");
						o.put("code", code);
						String decoded = icd10CodeDesc.get(code);
						if (decoded !=null) o.put("description", decoded);
						else IO.el("\t\tFailed to find ICD-10 Code '"+code+"' Add it! https://icdlist.com/icd-10/look-up");
						jo.append("icd_codes", o);
					}
				}	
				//ICD-O morphology - describes the cell type (or histology) of the tumor, together with the behavior (malignant or benign)
				if (s.getTempusIcdOCodeMorphology() != null) {
					for (String code: Misc.COMMA.split(s.getTempusIcdOCodeMorphology())) {
						JSONObject o = new JSONObject();
						o.put("type", "ICD-O Morphology");
						o.put("code", code);
						String decoded = icdOCodeMorphology.get(code);
						if (decoded !=null)  o.put("description", decoded);
						else IO.el("\t\tFailed to find ICD-O Morphology Code '"+code+"' Add it! https://cancercenter.ai/icd-o-pathology-codes/");
						jo.append("icd_codes", o);
					}
				}	

				//ICD-O topology - describes the anatomical site of origin (or organ system) of the tumor
				if (s.getTempusIcdOCodeTopography() != null) {
					for (String code: Misc.COMMA.split(s.getTempusIcdOCodeTopography())) {
						JSONObject o = new JSONObject();
						o.put("type", "ICD-O Topology");
						o.put("code", code);
						String decoded = icdOCodeTopology.get(code);
						if (decoded !=null)  o.put("description", decoded);
						else IO.el("\t\tFailed to find ICD-O Topology Code '"+code+"' Add it! https://cancercenter.ai/icd-o-pathology-codes/");
						jo.append("icd_codes", o);
					}
				}	
			}
			IO.writeString(jo.toString(3), new File(resultsDir, rand+".json"));
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new TempusPathoPrinter(args);
	}		

	/**This method will process each argument and assign new variables
	 * @throws Exception */
	public void processArgs(String[] args) throws Exception{

			IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
			Pattern pat = Pattern.compile("-[a-zA-Z]");
			File icd10File = null;
			File icdMorphologyFile = null;
			File icdTopologyFile = null;
			File jsonDir = null;
			for (int i = 0; i<args.length; i++){
				Matcher mat = pat.matcher(args[i]);
				if (mat.matches()){
					char test = args[i].charAt(1);
					try{
						switch (test){
						case 'j': jsonDir = new File(args[++i]); break;
						case 's': resultsDir = new File(args[++i]); break;
						case 'i': icd10File = new File(args[++i]); break;
						case 'm': icdMorphologyFile = new File(args[++i]); break;
						case 't': icdTopologyFile = new File(args[++i]); break;
						default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
						}
					}
					catch (Exception e){
						e.printStackTrace();
						Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
					}
				}
			}
			
			if (jsonDir == null || jsonDir.exists()== false) throw new IOException("ERROR: failed to find your dir of json files to parse.");
			jsonFiles = IO.extractFiles(jsonDir, "json");
			
			//icd files?
			IO.pl("Loading ICD 10 codes...");
			icd10CodeDesc = IO.loadFileIntoHash(icd10File, 0, 1);
			if (icd10CodeDesc == null) throw new IOException("ERROR: failed to pars ICD Diagnosis file "+icd10File);
			IO.pl("Loading ICD-O morphology codes...");
			icdOCodeMorphology = IO.loadFileIntoHash(icdMorphologyFile, 0, 1);
			if (icdOCodeMorphology == null) throw new IOException("ERROR: failed to pars ICD morphology file "+icdMorphologyFile);
			IO.pl("Loading ICD-O topology codes for SpecimenSites...");
			icdOCodeTopology = IO.loadFileIntoHash(icdTopologyFile, 0, 1);
			if (icdOCodeTopology == null) throw new IOException("ERROR: failed to pars ICD topology file "+icdTopologyFile);
			
			if (resultsDir == null) throw new IOException("ERROR: failed to find your dir to save the parsed output.");
			resultsDir.mkdirs();
						
	}

	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                           Tempus Patho Printer : Jan 2025                        **\n" +
				"**************************************************************************************\n" +
				"\n"+
				
				"\nOptions:\n"+
				"-j Directory containing Tempus v3.3+ json reports.\n" +
				"-s Directory to save the results.\n"+
				"-i File containing ICD-10 diagnosis codes and their descriptors. Tab delimited.\n"+
				"-m File containing ICD-O morphology codes. Ditto.\n"+
				"-t File containing ICD-O topology codes. Ditto.\n"+
				
				"\nExample: java -jar pathToUSeq/Apps/TempusPathoPrinter \n"+

				"\n**************************************************************************************\n");
	}
	
}
