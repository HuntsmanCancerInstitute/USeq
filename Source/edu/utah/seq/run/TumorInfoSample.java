package edu.utah.seq.run;

import java.io.File;
import java.util.HashMap;

import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class TumorInfoSample {
	
	//fields
	private File sampleDirectory = null;
	private HashMap<String, File> dirNameFile = null;
	private String msi = null;
	private String msiLine = null;
	private String loh = null;
	private String lohLine = null;
	private String copyRatio = null;
	private String tmb = null;
	private String tmbLine = null;
	
	private String numSomVars = null;
	private String numGermVars = null;

	private static final String svcDir = "SomaticVariantCalls";
	private static final String msiExt = "_Msi";
	private static final String mantisExt = "_Mantis.txt";
	private static final String illExt = "_Illumina";
	private static final String tmbExt = "_Tmb.txt";
	private static final String bkzExt = "_VCFBkz.log";
	
	private static final String caDir = "CopyAnalysis";
	private static final String lohExt = "_LoH";
	private static final String lohLogExt = "_LoH.log";
	private static final String crExt = "_GATKCopyRatio";
	private static final String crFrcExt = "fractionCopyRatio.txt";
	
	private static final String germDir = "GermlineVariantCalling";
	private static final String illAnnDir = "_Illumina_Anno";
	private static final String annVcfParser = "_AnnotatedVcfParser.log";
	
	public static final String toStringHeader = "Name\tMSI\tTMB\tFracCopyRatio\tFracLoH\t#SomCalls\t#GermCalls";
	
	
	public TumorInfoSample (File sampDir) {
		this.sampleDirectory = sampDir;
		dirNameFile = IO.extractDirectories(sampleDirectory);
		
		parseMsi();
		parseTmbSomVar();
		parseLoh();
		parseCopyRatio();
		parseNumGermVars();
	}
	

	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		//name
		sb.append(sampleDirectory.getName()); sb.append("\t");
		//msi
		sb.append(msi); sb.append("\t");
		//tmb
		String tmbTrim = null;
		if (tmb != null) {
			double tmbNum = Double.parseDouble(tmb);
			tmbTrim = Num.formatNumber(tmbNum, 1);
		}
		sb.append(tmbTrim); sb.append("\t");
		//copyRatio
		String crTrim = null;
		if (copyRatio != null) {
			double crNum = Double.parseDouble(copyRatio);
			crTrim = Num.formatNumber(crNum, 3);
		}
		sb.append(crTrim); sb.append("\t");
		//loh
		String lohTrim = null;
		if (loh != null) {
			double lohNum = Double.parseDouble(loh);
			lohTrim = Num.formatNumber(lohNum, 3);
		}
		sb.append(lohTrim); sb.append("\t");
		//# som
		sb.append(numSomVars); sb.append("\t");
		//# germ
		sb.append(numGermVars);
		return sb.toString();
	}
	
	private void parseCopyRatio() {
		//any Copy Analysis?
		if (dirNameFile.containsKey(caDir)==false) return;
		
		//pull GATKCopyRatio dir
		File[] dirs = IO.extractOnlyDirectories(dirNameFile.get(caDir));
		File crDir = null;
		for (File f: dirs) {
			if (f.getName().endsWith(crExt)) {
				crDir = f;
				break;
			}
		}
		if (crDir == null) return;
		
		//pull Results dir
		File resDir = new File(crDir, "Results");
		if (resDir.exists() == false) return;
		
		//pull frac
		File[] resFiles = IO.extractOnlyFiles(resDir);
		File fracFile = null;
		for (File f: resFiles) {
			if (f.getName().endsWith(crFrcExt)) {
				fracFile = f;
				break;
			}
		}
		if (fracFile == null) return;
		
		//parse results file
		String[] fracResults = IO.loadFile(fracFile);
		String fracLine = null;
		for (String r: fracResults) {
			if (r.contains("CNV frac")) {
				fracLine = r;
				break;
			}
		}
		if (fracLine == null) return;
		
		//0.410371273294813 CNV frac
		String[] tokens = Misc.WHITESPACE.split(fracLine);
		copyRatio = tokens[0];
	}

	private void parseLoh() {
		//any SomaticVariantCalls
		if (dirNameFile.containsKey(caDir)==false) return;
		
		//pull LoH dir
		File[] dirs = IO.extractOnlyDirectories(dirNameFile.get(caDir));
		File lohDir = null;
		for (File f: dirs) {
			if (f.getName().endsWith(lohExt)) {
				lohDir = f;
				break;
			}
		}
		if (lohDir == null) return;
		
		//pull Log dir
		File logDir = new File(lohDir, "Logs");
		if (logDir.exists() == false) return;
		
		//pull loh log
		File[] logFiles = IO.extractOnlyFiles(logDir);
		File lohFile = null;
		for (File f: logFiles) {
			if (f.getName().endsWith(lohLogExt)) {
				lohFile = f;
				break;
			}
		}
		if (lohFile == null) return;
		
		//parse loh log file
		String[] lohResults = IO.loadFile(lohFile);
		for (String r: lohResults) {
			if (r.contains("significant increase")) {
				lohLine = r;
				break;
			}
		}
		if (lohLine == null) return;
		
		//0.08659602 (5295/61146) : Fraction heterozygous germline snvs and indels with a significant increase in their allele fraction in the tumor
		String[] tokens = Misc.WHITESPACE.split(lohLine);
		loh = tokens[0].trim();;
		
	}

	private void parseNumGermVars() {
		//any GermlineVariantCalling
		if (dirNameFile.containsKey(germDir)==false) return;

		//pull Illumina Anno dir
		File[] dirs = IO.extractOnlyDirectories(dirNameFile.get(germDir));
		File illDir = null;
		for (File f: dirs) {
			if (f.getName().endsWith(illAnnDir)) {
				illDir = f;
				break;
			}
		}
		if (illDir == null) return;

		//pull Log dir
		File logDir = new File(illDir, "Logs");
		if (logDir.exists() == false) return;

		//pull anno parser results
		File[] logFiles = IO.extractOnlyFiles(logDir);
		File log = null;
		for (File f: logFiles) {
			if (f.getName().endsWith(annVcfParser)) {
				log = f;
				break;
			}
			
		}
		if (log != null) {
			//parse  file
			String[] results = IO.loadFile(log);
			String germlineLine = null;
			for (String r: results) {
				if (r.contains("Number of VCF Records")) {
					germlineLine = r.trim();
					break;
				}
			}
			if (germlineLine != null) {
				//113988	Number of VCF Records
				String[] tokens = Misc.TAB.split(germlineLine);
				numGermVars = tokens[0];
			}
		}
	}

	
	
	
	private void parseTmbSomVar() {
		//any SomaticVariantCalls
		if (dirNameFile.containsKey(svcDir)==false) return;

		//pull Illumina dir
		File[] svcDirs = IO.extractOnlyDirectories(dirNameFile.get(svcDir));
		File illDir = null;
		for (File f: svcDirs) {
			if (f.getName().endsWith(illExt)) {
				illDir = f;
				break;
			}
		}
		if (illDir == null) return;

		//pull Log dir
		File logDir = new File(illDir, "Logs");
		if (logDir.exists() == false) return;

		//pull tmb results
		File[] logFiles = IO.extractOnlyFiles(logDir);
		File tmbFile = null;
		File bkzFile = null;
		for (File f: logFiles) {
			if (f.getName().endsWith(tmbExt)) tmbFile = f;
			else if (f.getName().endsWith(bkzExt)) bkzFile = f;
		}
		if (tmbFile != null) {
			//parse tmb file
			String[] tmbResults = IO.loadFile(tmbFile);
			tmbLine = null;
			for (String r: tmbResults) {
				if (r.contains("tmb")) {
					tmbLine = r;
					break;
				}
			}
			if (tmbLine != null) {
				//tmb 96.6546
				String[] tokens = Misc.WHITESPACE.split(tmbLine);
				tmb = tokens[1];
				//reset tmbLine, want merge
				tmbLine = Misc.stringArrayToString(tmbResults, ", ");
			}
		}
		
		if (bkzFile != null) {
			//parse bkzFile file
			String[] results = IO.loadFile(bkzFile);
			String toParse = null;
			for (int i=0; i< results.length; i++) {
				if (results[i].contains("FailingBKZScore")) {
					toParse = results[i+1];
					break;
				}
			}
			if (toParse != null) {
				//Name	Records	Saved	NotScored	FailingBKZScore	FailingBKAFFilter
				//ZZ6RIP0UFE_TWSv2_FT-SA160156_FT-SA160113D_FT-SA160113R_Illumina_Hg38_Strelka.raw_Filtered.vcf.gz	9770	8297	172	192	1281
				//											0														  1		  2
				String[] tokens = Misc.WHITESPACE.split(toParse);
				numSomVars = tokens[2];
			}
		}
	}


	private void parseMsi() {
		//any SomaticVariantCalls
		if (dirNameFile.containsKey(svcDir)==false) return;
		
		//pull msi dir
		File[] svcDirs = IO.extractOnlyDirectories(dirNameFile.get(svcDir));
		File msiFile = null;
		for (File f: svcDirs) {
			if (f.getName().endsWith(msiExt)) {
				msiFile = f;
				break;
			}
		}
		if (msiFile == null) return;
		
		//pull msi results
		File mantis = null;
		File[] msiFiles = IO.extractOnlyFiles(msiFile);
		for (File f: msiFiles) {
			if (f.getName().endsWith(mantisExt)) {
				mantis = f;
				break;
			}
		}
		if (mantis == null) return;
		
		//parse mantis file, looking for 'Step-Wise Difference (DIF)	0.3702	0.4000   	Stable'
		String[] mantisResults = IO.loadFile(mantis);
		for (String r: mantisResults) {
			if (r.contains("(DIF)")) {
				msiLine = r;
				break;
			}
		}
		if (msiLine == null) return;
		//Step-Wise Difference (DIF)	0.3702	0.4000   	Stable
		String[] tokens = Misc.TAB.split(msiLine);
		msi = tokens[tokens.length-1].trim();
		
	}
	
	public static void main (String[] args) {
		File dir = new File ("/Users/u0028003/Downloads/JobDirDelme/ZZ6RIP0UFE_TWSv2_FT-SA160156_FT-SA160113D_FT-SA160113R");
		TumorInfoSample tis = new TumorInfoSample(dir);
		IO.pl("Name\tMSI\tTMB\tFracCopyRatio\tFracLoH\t#SomCalls\t#GermCalls");

		IO.pl(tis);
	}

}
