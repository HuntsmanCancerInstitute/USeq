package edu.utah.seq.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;

public class VariantInfoSample {
	String id = null;
	
	//germline
	private String germlineLabel = "Normal";
	private String germlineVcfMatch = "_JointGenotyped.vcf.gz";
	private String germlineVcfMatchROI = "_Hg38.anno.filt.roi.vcf.gz";
	private String germlineCoverageBedMatch = "_Pass.bed.gz";
	
	//somatic
	private String somaticLabel = "Tumor";
	private String somaticVcfMatch = "_final.vcf.gz";
	private String somaticCoverageBedMatch = "_CoveredRegion.bed.gz";
	
	//for both, use label to figure out which it is
	private String annoFiltVcfMatch = ".anno.filt.vcf.gz";
	
	//json info file
	private String jsonInfoMatch = "_AvatarInfo.json.gz";
	
	//files
	private File jsonInfo = null;
	private File germlineVcf = null;
	private File germlineBed = null;
	private File germlineAnnoVcf = null;
	private File germlineAnnoVcfROI = null;
	private File somaticVcf = null;
	private File somaticBed = null;
	private File somaticAnnoVcf = null;
	
	//counts
	private int numGermlineVars = -1;
	private int numGermlineBps = -1;
	private int numAnnoFiltGermlineVars = -1;
	private int numAnnoFiltGermlineVarsROI = -1;
	private int numSomaticVars = -1;
	private int numSomaticBps = -1;
	private int numAnnoFiltSomaticVars = -1;
	private float tmb;
	private float gmb;
	
	private Pattern infoVal = Pattern.compile(".+:\\s+\"(.+)\",*");
	private String diseaseType = null;
	private String analysisId = null;
	
	//gene names
	private TreeMap<String, Integer> germlineGenes = new TreeMap<String, Integer>();
	private TreeMap<String, Integer> germlineGenesROI = new TreeMap<String, Integer>();
	private TreeMap<String, Integer> somaticGenes = new TreeMap<String, Integer>();
	
	public VariantInfoSample(File dir) {
		this.id = dir.getName();
		
		findFiles(dir);
		
		countVcfs();
		
		countBps();
		
		loadGenes();
		
		tmb= 1000000f * (float)numSomaticVars/ (float)numSomaticBps;
		gmb= 1000000f * (float)numGermlineVars/ (float)numGermlineBps;
		
		fetchJsonInfo();
	}
	
	public static final String header = "ID\tDisease\tAnalysis ID\tGermline Bps\tGermline Vars\tFilt Germline Vars\tFilt Germline ROI Vars\tGMB\tFilt Germline Genes: Hits\tFilt Germline ROI Genes: Hits\tSomatic Bps\tSomatic Vars\tFilt Somatic Vars\tTMB\tFilt Somatic Genes: Hits";
	
	public String toString() {
		ArrayList<String> al = new ArrayList<String>();
		//id
		al.add(id);
		//diseaseType
		al.add(diseaseType);
		//analysisId
		al.add(analysisId);
		
		//numGermlineBps
		al.add(new Integer(numGermlineBps).toString());
		//numGermlineVcfs
		al.add(new Integer(numGermlineVars).toString());
		//numFiltGermlineVcfs
		al.add(new Integer(numAnnoFiltGermlineVars).toString());
		//numFiltGermlineVcfsROI
		al.add(new Integer(numAnnoFiltGermlineVarsROI).toString());
		//gmb
		al.add(new Float(gmb).toString());
		//filtGermlineGenes
		al.add(treeMapToString(germlineGenes));
		//filtGermlineGenesROI
		al.add(treeMapToString(germlineGenesROI));
		
		//numSomaticBps
		al.add(new Integer(numSomaticBps).toString());
		//numSomaticVcfs
		al.add(new Integer(numSomaticVars).toString());
		//numFiltSomaticVcfs
		al.add(new Integer(numAnnoFiltSomaticVars).toString());
		//tmb
		al.add(new Float(tmb).toString());
		//filtSomaticGenes
		al.add(treeMapToString(somaticGenes));
		
		return Misc.stringArrayListToString(al, "\t");
	}
	
	private static String treeMapToString(TreeMap<String, Integer> t) {
		if (t.size()==0) return "";
		ArrayList<String> al = new ArrayList<String>();
		for (String key: t.keySet()) {
			Integer i = t.get(key);
			al.add(key+":"+i);
		}
		return Misc.stringArrayListToString(al, " ");
	}
	
	private void fetchJsonInfo() {
		if (jsonInfo == null) return;

		try {
			//looking for
			// "Diagnosis": "HEM-CLL"
			// "AnalysisId": "A1133"
			BufferedReader in = IO.fetchBufferedReader(jsonInfo);
			String line = null;
			while ((line = in.readLine()) != null){
				if (line.contains("Diagnosis")) {					
					Matcher m = infoVal.matcher(line);
					if (m.matches()) diseaseType = m.group(1);
				}
				else if (line.contains("AnalysisId")) {
					Matcher m = infoVal.matcher(line);
					if (m.matches()) analysisId = m.group(1);
				}
			}
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("Problem loading json info "+jsonInfo);
		}
	}

	private void loadGenes() {
		try {
			if (germlineAnnoVcf != null) loadGeneNames(germlineAnnoVcf, germlineGenes);
			if (germlineAnnoVcfROI != null) loadGeneNames(germlineAnnoVcfROI, germlineGenesROI);
			if (somaticAnnoVcf != null) loadGeneNames(somaticAnnoVcf, somaticGenes);
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("Problem loading gene names.");
		}
	}
	
	private void loadGeneNames(File vcf, TreeMap<String, Integer> genes) throws IOException {
		BufferedReader in = IO.fetchBufferedReader(vcf);
		String line = null;
		HashSet<String> lineLevelGenes = new HashSet<String>();
		while ((line = in.readLine()) != null) if (line.trim().length() !=0 && line.startsWith("#") == false) {
			//split on tab
			String[] tabFields = Misc.TAB.split(line);
			String info = tabFields[7];
			//split on ;
			String[] infoFields = Misc.SEMI_COLON.split(info);
			//find ANN=
			String ann = null;
			for (String s: infoFields) {
				if (s.startsWith("ANN=")) {
					ann = s;
					break;
				}
			}
			if (ann == null) continue; //these exist sometimes, throw new IOException("Failed to find a ANN= field in "+line+" from "+vcf);
			//split on comma
			for (String transcript: Misc.COMMA.split(ann)){
				String[] transFields = Misc.PIPE.split(transcript);
				for (int i=0; i< transFields.length; i++) {
					if (transFields[i].startsWith("ENSG")) {
						lineLevelGenes.add(transFields[i-1]);
						break;
					}
				}
			}
			//add in line level to global, increment counter
			for (String lg: lineLevelGenes) {
				Integer count = genes.get(lg);
				if (count == null) count = new Integer(1);
				else count = new Integer(count.intValue()+1);
				genes.put(lg, count);
			}
			lineLevelGenes.clear();
		}
		in.close();
	}

	private void countBps() {
		try {
			if (germlineBed != null) numGermlineBps = countBedRegionBps(germlineBed);
			if (somaticBed != null) numSomaticBps = countBedRegionBps(somaticBed);
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("Problem counting number of bed region bps.");
		}
		
	}

	private void countVcfs() {
		try {
			if (germlineVcf != null) numGermlineVars = countVcfRecords(germlineVcf);
			if (germlineAnnoVcf != null) numAnnoFiltGermlineVars = countVcfRecords(germlineAnnoVcf);
			if (germlineAnnoVcfROI != null) numAnnoFiltGermlineVarsROI = countVcfRecords(germlineAnnoVcfROI);
			if (somaticVcf != null) numSomaticVars = countVcfRecords(somaticVcf);
			if (somaticAnnoVcf != null) numAnnoFiltSomaticVars = countVcfRecords(somaticAnnoVcf);

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("Problem counting number of vcf records.");
		}
	}

	private int countVcfRecords(File vcf) throws IOException {
		int numVars = 0;
		BufferedReader in = IO.fetchBufferedReader(vcf);
		String line = null;
		while ((line = in.readLine()) != null) if (line.trim().length() !=0 && line.startsWith("#") == false) numVars++;
		in.close();
		return numVars;
	}
	
	private int countBedRegionBps(File bed) throws Exception {
		int num = 0;
		BufferedReader in = IO.fetchBufferedReader(bed);
		String line = null;
		while ((line = in.readLine()) != null) if (line.trim().length() !=0 && line.startsWith("#") == false) {
			String[] fields = Misc.TAB.split(line);
			int start = Integer.parseInt(fields[1]);
			int end = Integer.parseInt(fields[2]);
			num += (end-start);
		}
		in.close();
		return num;
	}

	private void findFiles(File dir) {
		//walk all files in the dir looking to assign them
		File[] files = IO.fetchFilesRecursively(dir);
		for (File f: files) {
			String name = f.getName();
			//jsonInfo
			if (name.endsWith(jsonInfoMatch)) jsonInfo = f;
			//germlineVcfMatch
			else if (name.endsWith(germlineVcfMatch)) germlineVcf = f;
			//germlineVcfMatchROI
			else if (name.endsWith(germlineVcfMatchROI) && name.contains(germlineLabel)) germlineAnnoVcfROI = f;
			//germlineCoverageBed
			else if (name.endsWith(germlineCoverageBedMatch) && name.contains(germlineLabel)) germlineBed = f;
			//somaticVcfMatch
			else if (name.endsWith(somaticVcfMatch)) somaticVcf = f;
			//somaticCoverageBed
			else if (name.endsWith(somaticCoverageBedMatch)) somaticBed = f;
			//annoFiltVcfMatch Tumor or Normal
			else if (name.endsWith(annoFiltVcfMatch)) {
				if (name.contains(germlineLabel)) germlineAnnoVcf = f;
				else somaticAnnoVcf = f;
			}
		}
	}

	public String getId() {
		return id;
	}

	public int getNumGermlineVars() {
		return numGermlineVars;
	}

	public int getNumGermlineBps() {
		return numGermlineBps;
	}

	public int getNumAnnoFiltGermlineVars() {
		return numAnnoFiltGermlineVars;
	}

	public int getNumSomaticVars() {
		return numSomaticVars;
	}

	public int getNumSomaticBps() {
		return numSomaticBps;
	}

	public int getNumAnnoFiltSomaticVars() {
		return numAnnoFiltSomaticVars;
	}

	public float getTmb() {
		return tmb;
	}

	public float getGmb() {
		return gmb;
	}

	public String getDiseaseType() {
		return diseaseType;
	}

	public TreeMap<String, Integer> getGermlineGenes() {
		return germlineGenes;
	}

	public TreeMap<String, Integer> getSomaticGenes() {
		return somaticGenes;
	}
}
