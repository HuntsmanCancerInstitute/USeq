package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import util.gen.*;

/** Wraps a smallish Pipeline.jar template for QC, annotating, and run directory generation.
 * @author Nix
 * */
public class ArupPipelineWrapper {

	//fields
	private String jobId = "";
	private String sampleId = "";
	private String submitter = "";
	private String analysisType = "";
	private String minimumReadDepth = "";
	
	//defaults
	private boolean uploadVarsToNGSWeb = false;
	private String snpEffGenome = "hg19_ucsc_20150427";
	
	//files
	private File pJar = null;
	private File pProperties = null;
	private File bedForCoverageQC = null;
	private File bedForVarCalling = null;
	private File bedForBadRegions = null;
	private File fastaReference = null;
	private File unfilteredBam = null;
	private File finalBam = null;
	private File finalVcf = null;
	private File webRootForLinks = null;
	private File outputDirectory = null;
	
	//internal
	private String xmlTemplate = null;
	private boolean deleteTempVcf = false;
	
	
	//constructor
	public ArupPipelineWrapper(String[] args){	
		long startTime = System.currentTimeMillis();
		
		processArgs(args);
		
		buildXmlTemplate();
		
		executePipelineJob();
		
		deleteTempFiles();
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}

	private void executePipelineJob() {
		String[] cmd = null;
		try {
			//write out tempTemplate
			File template = new File (outputDirectory, "pipelineTemplate.xml");
			if (IO.writeString(xmlTemplate, template) == false) Misc.printErrAndExit("Problem writing -> "+template);
			
			//build and execute cmd
			cmd = new String[]{"java", "-jar", "-Xmx2G", pJar.getCanonicalPath(), "-props", pProperties.getCanonicalPath(), template.getCanonicalPath()};			
			
			System.out.println("\nExecuting:\n"+Misc.stringArrayToString(cmd, " "));

			String cmdLineOutput = IO.executeCommandReturnError(cmd);
			//check and write output to stdout
			if (cmdLineOutput.contains("ERROR")) Misc.printErrAndExit(cmdLineOutput+"\n\nERROR found in Pipeline.jar output, see above. Aborting!");
			System.out.println("\nPipelineOutput:\n"+cmdLineOutput);
			
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("ERROR: executing "+Misc.stringArrayToString(cmd, " "));
		}
	}
	
	//remove temp files, the pipeline.jar operator isn't working.
	public void deleteTempFiles(){
		System.out.println("\nRemoving these temp files:");
		File workingDir = new File (System.getProperty("user.dir"));
		File[] toExamine = workingDir.listFiles();
		for (File f: toExamine){
			boolean d = false;
			String n = f.getName();
			if (n.startsWith("pipeinstancelog")) d=true;
			else if (n.contains(".DOC.sample")) d=true;
			else if (n.contains("allDepths.")) d=true;
			else if (n.startsWith("nocalls.")) d=true;
			else if (n.startsWith("snpeff.")) d=true;
			else if (n.contains("plice")) d=true;
			if (d){
				System.out.println("\t"+n);
				f.deleteOnExit();
			}
		}
		//delete the temp uncompressed vcf (required by Pipeline.jar)
		if (deleteTempVcf) System.out.println("\t"+finalVcf.getName());
	}


	private void buildXmlTemplate() {
		StringBuilder sb = new StringBuilder();
		sb.append("<Pipeline>\n\n");
		sb.append(fetchFieldXml());
		sb.append(fetchQCXml());
		sb.append(fetchVariantPoolsXml());
		sb.append(fetchAnnotatorXml());
		if (uploadVarsToNGSWeb) sb.append(fetchVariantUploadXml());
		sb.append(fetchVariantOutputXml());
		sb.append(fetchCleanupReviewDirLinkGenXml());
		if (webRootForLinks != null) fetchWebRootLinksXml();
		sb.append("\n</Pipeline>\n");
		xmlTemplate = sb.toString();
	}
	
	public String fetchFieldXml(){
		String s = 
			"<bedForCoverageQC class='buffer.BEDFile' filename='"+bedForCoverageQC+"' />\n" +
			"<bedForVars class='buffer.BEDFile' filename='"+ bedForVarCalling +"' />\n" +
			"<bedForBadRegions class='buffer.BEDFile' filename='"+ bedForBadRegions +"' />\n" +
			"<reference class='buffer.ReferenceFile' filename='"+ fastaReference +"' />\n" +
			"<finalBAM class='buffer.BAMFile' filename='"+ finalBam +"' />\n" +
			"<unfilteredBAM class='buffer.BAMFile' filename='"+ unfilteredBam +"' />\n" +
			"<finalVariants class='buffer.VCFFile' filename='"+ finalVcf +"' />\n"+
			"\n";
		return s;
	}
	
	public String fetchAnnotatorXml(){
		String s =
				"<SnpEffAnnotate class='operator.snpeff.SnpEffGeneAnnotate' snpeff.genome='"+snpEffGenome +"' > \n"+
				"	<VariantPool /> \n"+
				"</SnpEffAnnotate> \n"+
				
				"<NovelFilter class='operator.variant.NovelFilter'> \n"+
				"	<input> \n"+
				"		<VariantPool /> \n"+
				"	</input> \n"+
				"	<output> \n"+
				"		<NovelVariants class='buffer.variant.VariantPool' /> \n"+
				"	</output> \n"+
				"</NovelFilter> \n"+
				
				"<dbSNP class='operator.variant.DBSNPAnnotator'> \n"+
				"	<VariantPool /> \n"+
				"</dbSNP> \n"+
				
				"<VarsToGenes class='operator.gene.GeneListFromVariants'> \n"+
				"	<VariantPool /> \n"+
				"	<Genes class='buffer.GeneList' /> \n"+
				"</VarsToGenes> \n"+
				
				"<TKGAnnotate class='operator.variant.TGPTabixAnnotator'> \n"+
				"	<VariantPool /> \n"+
				"</TKGAnnotate> \n"+
				
				//This requires installation of List/MoreUtils.pm and Math/Round.pm and setting them in the $PERL5LIB path.
				//It also requires hard coding the path to the splice models in the scoreSpliceSites perl app, see $dir.  
				//NOT good! 
				"<SpliceAnnotate class='operator.variant.SplicingPredictionAnnotator'> \n"+
				"	<VariantPool /> \n"+
				"</SpliceAnnotate> \n"+
				
				"<ExAC63KAnnotate class='operator.variant.ExAC63KExomesAnnotator' justLoadOverall='true' > \n"+
				"	<VariantPool /> \n"+
				"</ExAC63KAnnotate> \n"+
				
				"<UK10KAnnotate class='operator.variant.UK10KAnnotator'> \n"+
				"	<VariantPool /> \n"+
				"</UK10KAnnotate> \n"+
				
				"<ESPAnnotate class='operator.variant.ESP6500Annotator'> \n"+
				"	<VariantPool /> \n"+
				"</ESPAnnotate> \n"+
				
				"<ARUPFreqAnno class='operator.variant.ARUPDBAnnotate'> \n"+
				"	<VariantPool /> \n"+
				"</ARUPFreqAnno> \n"+
				
				"<HGMDVar class='operator.variant.HGMDVarAnnotator'> \n"+
				"	<VariantPool /> \n"+
				"</HGMDVar> \n"+
				
				"<COSMICAnnotate class='operator.variant.COSMICAnnotator'>  \n"+
				"	<VariantPool /> \n"+
				"</COSMICAnnotate> \n"+
				
				//Sift, PolyPhen, MutTast, Gerp.... Very OLD like 2012.  Is this needed given snpEff?
				"<dbNSFPAnnotate class='operator.variant.DBNSFPAnnotator'> \n"+
				"	<VariantPool /> \n"+
				"</dbNSFPAnnotate> \n"+
				
				"<BadRegions class='plugins.annotators.regions.BadRegionAnnotator'> \n"+
				"	<VariantPool /> \n"+
				"	<bedForBadRegions /> \n"+
				"</BadRegions> \n"+
				
				//Gene annotators
				"<HGMDGene class='operator.gene.HGMDAnnotator'> \n"+
				"	<Genes /> \n"+
				"</HGMDGene> \n"+

				"<OMIM class='operator.gene.OMIMAnnotator'> \n"+
				"	<Genes /> \n"+
				"</OMIM> \n"+
				
				//Disease, function, expression... New from Keith
				"<dbNSFPGene class='operator.gene.DBNSFPGeneAnnotator'> \n"+
				"	<Genes /> \n"+
				"</dbNSFPGene> \n"+
				"\n";
				return s;
	}
	
	public String fetchQCXml(){
		String s =
			"<ParallelQCOps class='operator.ParallelOperator'> \n"+
			
			"<NoCalls class='operator.gatk.CallableLoci' min.depth='"+ minimumReadDepth +"'> \n"+
			" <input> \n"+
			"    <reference /> \n"+
			"    <finalBAM /> \n"+
			"    <bedForCoverageQC /> \n"+
			" </input> \n"+
			" <output> \n"+
			"  <NoCallCSV class='buffer.CSVFile' filename='nocalls.bed' /> \n"+
			" </output> \n"+
			"</NoCalls> \n"+
			
			"<ComputeFinalBamMetrics class='operator.qc.BamMetrics' > \n"+
			" <input> \n"+
			"  <finalBAM /> \n"+
			" </input> \n"+
			" <output> \n"+
			"  <finalBamMetrics class='buffer.BAMMetrics' filename='final.bam.metrics.txt' /> \n"+
			" </output> \n"+
			"</ComputeFinalBamMetrics> \n"+
			
			"<ComputeRawBamMetrics class='operator.qc.BamMetrics'> \n"+
			" <input> \n"+
			"  <unfilteredBAM /> \n"+
			" </input> \n"+
			" <output> \n"+
			"  <rawBamMetrics class='buffer.BAMMetrics' filename='final.bam.metrics.txt' /> \n"+
			" </output> \n"+
			"</ComputeRawBamMetrics> \n"+
			
			"<rawdoc class='operator.gatk.DepthOfCoverage'> \n"+
			" <input> \n"+
			"  <reference /> \n"+
			"  <unfilteredBAM /> \n"+
			"  <bedForCoverageQC /> \n"+
			" </input> \n"+
			" <output> \n"+
			"   <rawdocmetrics class='buffer.DOCMetrics' /> \n"+
			" </output> \n"+
			"</rawdoc> \n"+
			
			"<finaldoc class='operator.gatk.DepthOfCoverage'> \n"+
			" <input> \n"+
			"  <reference /> \n"+
			"  <finalBAM /> \n"+
			"  <bedForCoverageQC /> \n"+
			" </input> \n"+
			" <output> \n"+
			"   <finaldocmetrics class='buffer.DOCMetrics' /> \n"+
			" </output> \n"+
			"</finaldoc> \n"+
			
			"</ParallelQCOps> \n"+
			
			"<CallableDepth class='operator.qc.DepthsForNoCalls'> \n"+
			"    <input> \n"+
			"       <reference /> \n"+
			"       <NoCallCSV /> \n"+
			"       <finalBAM /> \n"+
			"    </input> \n"+
			"    <output> \n"+
			"       <NoCallDepths class='buffer.CSVFile' filename='noCallDepths.csv' /> \n"+
			"    </output> \n"+
			"</CallableDepth>\n"+
			"\n";
		return s;
	}
	
	public String fetchVariantPoolsXml(){
		String s = 
				// Create variant pool from vcf
				"<VariantPool class='buffer.variant.VariantPool'>\n" +
				"   <finalVariants />\n" +
				"</VariantPool>\n" +
				
				// Create variants pool for qc analysis 
				"<VariantPoolQCFiltered class='operator.variant.FilterPoolByBED'>\n" +
				"  <input>\n" +
				"       <bedForCoverageQC />\n" +
				"       <VariantPool />\n" +
				"    </input>\n" +
				"    <output>\n" +
				"        <filteredVars class='buffer.variant.VariantPool' />\n" +
				"    </output>\n" +
				"</VariantPoolQCFiltered>\n"+
				"\n";
		return s;
	}
	
	public String fetchVariantOutputXml(){
		String s =
				
		//Write variants to csv file and include bad.region 
		"<ViewerFile class='plugins.writers.varviewer.VarViewerWriter' anno.keys='bad.region'> \n"+
		"	<VariantPool /> \n"+
		"	<Genes /> \n"+
		"	<ViewerCSV class='buffer.CSVFile' filename='"+ outputDirectory.getName()+"_annotated.csv' /> \n"+
		"</ViewerFile> \n"+
				
		//Create Final JsonQC Section (needs to happen after variant calls)
		"<QCtoJSON class='operator.qc.QCtoJSON'> \n"+
		"	<rawBamMetrics /> \n"+
		"	<rawdocmetrics /> \n"+
		"	<finalBamMetrics /> \n"+
		"	<finaldocmetrics /> \n"+
		"	<bedForCoverageQC /> \n"+
		"	<filteredVars /> \n"+
		"	<NoCallCSV /> \n"+
		"	<QCOutputFile class='buffer.TextBuffer' filename='qc.json' /> \n"+
		"</QCtoJSON> \n"+
				"\n";
		return s;
	}
	
	public String fetchVariantUploadXml(){
		String s =
		//Upload variants to NGS.Web
		"<VarUploader class='operator.variant.VariantUploader' sampleID='"+ sampleId+"'> \n"+
		"	<VariantPool /> \n"+
		"</VarUploader> \n\n";
				
		return s;
	}
	
	public String fetchWebRootLinksXml(){
		String s =
				"<CreateBAMLink class='operator.LinkCreator' sample='"+ outputDirectory.getName()+"' web.root='"+ webRootForLinks+"' result.dir='links/'> \n"+
				"        <finalBAM /> \n"+
				"        <bedForVars /> \n"+
				"</CreateBAMLink> \n\n";
				
		return s;
	}
	
	public String fetchCleanupReviewDirLinkGenXml(){
		String s = null;
		try {
			s = //Organize and Clean Up Section
			"<ReviewDirGenerator class='operator.ReviewDirGenerator' destination.dir='"+ outputDirectory.getCanonicalPath()+
					"' sample='"+ outputDirectory.getName()+"' submitter='"+ submitter+
					"' analysis.type='"+ analysisType+"'> \n"+
			"    <finalVariants /> \n"+
			"    <ViewerCSV /> \n"+
			"    <InstanceLog class='buffer.InstanceLogFile' /> \n"+
			"    <finalBAM /> \n"+
			"    <QCOutputFile /> \n"+
			"    <bedForVars /> \n"+
			"</ReviewDirGenerator> \n\n";
			
			//delete files, hmm for some reason this is not working, skipping and running own cleanup
			/*
			"<Cleanup class='operator.RemoveFile'> \n"+
			"	<realignTargets /> \n"+
			"	<gatkDepthFiles class='buffer.GlobFileBuffer' filename='.*.DOC.sample_.*' /> \n"+
			"	<gatkNoCallDepthFiles class='buffer.GlobFileBuffer' filename='nocallDepths.sample_.*' /> \n"+
			"	<snpeffFiles class='buffer.GlobFileBuffer' filename='snpeff.*' /> \n"+
			"	<noCallFiles class='buffer.GlobFileBuffer' filename='nocalls.*' /> \n"+
			"	<noCall2Files class='buffer.GlobFileBuffer' filename='noCallDepths.*' /> \n"+
			"	<spliceFiles class='buffer.GlobFileBuffer' filename='.*plice.*' /> \n"+
			"</Cleanup> \n\n";
			*/

		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError fetching clean up and review dir links.\n");
		}
		return s;
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ArupPipelineWrapper(args);
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
					case 'o': jobId = args[++i]; break;
					case 's': sampleId = args[++i]; break;
					case 'm': submitter = args[++i]; break;
					case 'y': analysisType = args[++i]; break;
					case 'w': webRootForLinks = new File(args[++i]); break;
					case 'e': snpEffGenome = args[++i]; break;
					case 't': minimumReadDepth = args[++i]; break;
					case 'l': uploadVarsToNGSWeb = true; break;
					case 'j': pJar = new File(args[++i]); break;
					case 'p': pProperties = new File(args[++i]); break;
					case 'q': bedForCoverageQC = new File(args[++i]); break;
					case 'b': bedForVarCalling = new File(args[++i]); break;
					case 'a': bedForBadRegions = new File(args[++i]); break;
					case 'r': fastaReference = new File(args[++i]); break;
					case 'u': unfilteredBam = new File(args[++i]); break;
					case 'f': finalBam = new File(args[++i]); break;
					case 'v': finalVcf = new File(args[++i]); break;
					case 'd': outputDirectory = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}	
		}
		
		
		
		//check output dir and if needed set sampleId
		if (outputDirectory !=null) {
			outputDirectory.mkdirs();
			if (sampleId.length() == 0) sampleId = outputDirectory.getName();
		}
		//check root directory if needed
		if (webRootForLinks!= null){
			if (webRootForLinks.exists() == false) webRootForLinks.mkdirs();
			//links dir?
			File links = new File (webRootForLinks, "links");
			if (links.exists() == false) links.mkdirs();
		}
		
		//look for required fields and files
		checkPrintFields();
		checkPrintFiles();
		checkPropFiles();
		checkForGzippedVcf();
		
		
	}

	private void checkForGzippedVcf() {
		String vcfName = finalVcf.getName();
		if (vcfName.endsWith(".gz")){
			File uncomp = new File (finalVcf.getParentFile(), vcfName.substring(0, vcfName.length()-3));
			if (uncomp.exists() == false) {
				deleteTempVcf = true;
				if (IO.uncompress(finalVcf, uncomp) == null) Misc.printErrAndExit("ERROR: failed to uncompress -> "+finalVcf);
			}
			finalVcf = uncomp;
		}
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           ArupPipelineWrapper: March 2016                        **\n" +
				"**************************************************************************************\n" +
				"This app wraps ARUP's pipeline.jar app for generating QC metrics, annotating variants,\n"+
				"and lastly creates the review directory.\n"+
				
				"\nParams:\n"+
				"  -o Job ID\n"+
				"  -m Submitter\n"+
				"  -y Analysis type\n"+
				"  -w Provide a root path for web links if you'd like to make them.\n"+
				"  -t Minimum alignment depth\n"+
				"  -s Sample ID, defaults to name of output directory.\n"+
				"  -d Path to the output directory\n"+
				"  -j Path to the pipeline.jar application\n"+
				"  -p Path to the pipeline properties file\n"+
				"  -q Path to the bed file for coverage QC\n"+
				"  -b Path to the bed file for variant calling\n"+
				"  -a Path to the bed file for flagging poor quality regions\n"+
				"  -r Path to the fasta reference file w/ index and dict\n"+
				"  -u Path to the unfiltered bam file\n"+
				"  -f Path to the filtered bam file\n"+
				"  -v Path to the final vcf file\n"+
				"  -e SnpEff genome, defaults to hg19_ucsc_20150427\n"+
				"  -l Upload variants to NGSWeb, defaults to not uploading\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/ArupPipelineWrapper -o MyJobNix3 -m DNix \n"+
				"    -j ~/BioApps/Pipeline-1.0-SNAPSHOT-jar-with-dependencies.jar -y TestAnaly -w \n"+
				"    ~/WebLinks -t 300 -d Results -p small_pipeline_properties.xml -q \n"+
				"    0758221_compPad25bp_v1.bed -b 0758221_v1.bed -a probBps0.05Pad1bp.bed -r \n"+
				"    ~/HCIAtlatl/data/Human/B37/human_g1k_v37_decoy.fasta -u CNV36B_unfiltered.bam -f \n"+
				"    CNV36B_final.bam -v CNV36B_snvIndel.vcf\n\n" +

				"**************************************************************************************\n");

	}
	
	public void checkPrintFields(){
		LinkedHashMap<String,String> nameField = new LinkedHashMap<String,String>();
		nameField.put("Job ID", jobId);
		nameField.put("Sample ID", sampleId);
		nameField.put("Submitter", submitter);
		nameField.put("Analysis Type", analysisType);
		nameField.put("Minimum Align Depth", minimumReadDepth);
		nameField.put("SnpEff Genome", snpEffGenome);
		//set output of booleans
		nameField.put("Upload Vars to NGSWeb", uploadVarsToNGSWeb+ "");
		
		boolean missingField = false;
		System.out.println("Fields:");
		for (String name: nameField.keySet()){
			String f = nameField.get(name);
			if (f == null || f.length() ==0) missingField = true;
			System.out.println(name+ "\t"+ f);
		}
		if (missingField) Misc.printErrAndExit("\nMissing Fields! See above.");
	}
	
	public void checkPrintFiles(){
		LinkedHashMap<String,File> nameFile = new LinkedHashMap<String,File>();
		nameFile.put("Pipeline.jar", pJar);
		nameFile.put("Pipeline properties", pProperties);
		nameFile.put("Coverage QC bed", bedForCoverageQC);
		nameFile.put("Variant calling bed", bedForVarCalling);
		nameFile.put("Poor quality region bed", bedForBadRegions);
		nameFile.put("Fasta genome reference", fastaReference);
		nameFile.put("Unfiltered bam", unfilteredBam);
		nameFile.put("Final bam", finalBam);
		nameFile.put("Final vcf", finalVcf);
		if (webRootForLinks != null) nameFile.put("Web links dir", webRootForLinks);
		nameFile.put("Final output dir", outputDirectory);
		
		boolean missingFile = false;
		System.out.println("\nResources (name exists path):");
		for (String name: nameFile.keySet()){
			File f = nameFile.get(name);
			if (f.exists() == false) missingFile = true;
			System.out.println(name+"\t"+f.exists()+"\t"+f);
		}
		if (missingFile) Misc.printErrAndExit("\nMissing resources! See above.");
	}
	
	private void checkPropFiles() {
		System.out.println("\nChecking your properties file -> "+pProperties);
		String[] prop = IO.loadFileIntoStringArray(pProperties);
		Pattern val = Pattern.compile("<entry key.+>(/.*)</entry>");
		boolean missingFile = false;
		for (String s: prop){
			Matcher mat = val.matcher(s);
			if (mat.matches()){
				File test = new File (mat.group(1));
				if (test.exists()) System.out.println("Found\t"+test);
				else {
					System.out.println("Misssing\t"+test);
					missingFile = true;
				}
			}
		}
		if (missingFile) Misc.printErrAndExit("\nFailed to find all of the files in your properties file, see above.\n");
	}

}
