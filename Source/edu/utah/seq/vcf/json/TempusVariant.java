package edu.utah.seq.vcf.json;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;

import org.json.JSONException;
import org.json.JSONObject;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Json;
import util.gen.Misc;
import util.gen.Num;

public class TempusVariant{
	
	//This is a composite of all of the variant fields in the 80 json reports we've received. Many of these fields will remain null
	private String variantSource = null;
	private String gene = null;
	private String transcript = null;
	private String mutationEffect = null;
	private String cHGVS = null;
	private String pHGVS = null;
	private String pFullHGVS = null;
	private String referenceGenome = null;
	private String variantType = null;
	private String variantDescription = null;
	private String nucleotideAlteration = null;
	private String chromosome = null;
	private String pos = null;
	private String ref = null;
	private String alt = null;
	private String allelicFraction = null;
	private String coverage = null;
	private String copyNumber = null;
	private String clinicalSignificance = null;
	private String disease = null;
	private String structuralVariant = null;
	private String fusionType = null;
	private String gene3 = null;
	private String gene5 = null;
	private String geneDescription = null;
	private TempusJson2Vcf tempusJson2Vcf = null;
	private String accessionId = null;
	
	//for CNV
	private int start;
	private int end;
	private int svLen; //will be negative if deletion
		
	/**Object to represent a Tempus variant, SNV/ INDEL/ CNV/ Fusion
	 * @throws IOException */
	public TempusVariant(String variantSource, String geneName, JSONObject object, TempusJson2Vcf tempusJson2Vcf) throws JSONException, IOException {
		this.variantSource = variantSource;
		this.tempusJson2Vcf = tempusJson2Vcf;
		//attempt to parse gene, if null pull from constructor
		gene = Json.forceGetString(object, "gene");
		if (gene == null) gene = geneName;

		mutationEffect = Json.forceGetString(object, "mutationEffect");
		transcript = Json.forceGetString(object, "transcript");
		cHGVS = Json.forceGetString(object, "HGVS.c");
		if (cHGVS == null) cHGVS = Json.forceGetString(object, "HGVS_c");
		pHGVS = Json.forceGetString(object, "HGVS.p");
		pFullHGVS = Json.forceGetString(object, "HGVS.pFull");
		referenceGenome = Json.forceGetString(object, "referenceGenome");

		variantType = Json.forceGetString(object, "variantType");
		variantDescription = Json.forceGetString(object, "variantDescription");

		nucleotideAlteration = Json.forceGetString(object, "nucleotideAlteration");
		chromosome = Json.forceGetString(object, "chromosome");
		pos = Json.forceGetString(object, "pos");
		ref = Json.forceGetString(object, "ref");
		alt = Json.forceGetString(object, "alt");

		allelicFraction = Json.forceGetString(object, "allelicFraction");
		coverage = Json.forceGetString(object, "coverage");

		copyNumber = Json.forceGetString(object, "copyNumber");

		clinicalSignificance = Json.forceGetString(object, "clinicalSignificance");
		disease = Json.forceGetString(object, "disease");

		// structural vars are all over the place so lots of edge cases
		structuralVariant = Json.forceGetString(object, "structuralVariant");
		fusionType = Json.forceGetString(object, "fusionType");
		gene3 = Json.forceGetString(object, "gene3");
		gene5 = Json.forceGetString(object, "gene5");
		if (gene3 != null && gene5 != null) variantType = "fusion";
		
		if (variantType == null && variantSource.toLowerCase().contains("fusion")) variantType = "fusion";
		
		geneDescription = Json.forceGetString(object, "geneDescription");
		
		//set coordinates from nucleotideAlteration?
		if (nucleotideAlteration != null){
			String[] t = Misc.COLON.split(nucleotideAlteration);
			if (t.length == 4 ) {
				chromosome = t[0];
				pos = t[1];
				ref = t[2];
				alt = t[3];
			}
		}
		
		//watch out for inherited variants
		if (variantType == null) {
			if (cHGVS != null && transcript != null && gene != null) variantType = "SNV";
			else throw new IOException("\nNull Variant Type\n"+this.toString());
		}
		
		//increment counters
		incrementCounters(tempusJson2Vcf);
		
	}
	

	public boolean addGenomicCoordinateInfo(String[] vcfLinesFirst, String[] vcfLinesSecond) throws IOException {
		if (gene == null || cHGVS == null || transcript == null) {
			Misc.printErrAndExit("\nMissing SNV Info for genomic coordinate parsing\n"+toString());
		}
		// Watch out for gene aliases, Tempus is switching these in the report from what is in the vcf annotations
		// They are also upper casing the gene names in the vcf
		String[] geneAliases = tempusJson2Vcf.getGeneAliases().get(gene);
		if (geneAliases == null) {
			//try just adding in the gene name from the json
			geneAliases = new String[]{gene};
			//throw new IOException("No gene aliases for "+gene + " from \n"+ this.toString());
			IO.el("\nWARNING: No gene aliases for "+gene + " from \n"+ this.toString());
		}
	
		Pattern[] genePats = new Pattern[geneAliases.length];
		for (int i=0; i< genePats.length; i++) genePats[i] = Pattern.compile(".*"+ geneAliases[i]+ ".*", Pattern.CASE_INSENSITIVE);
		
		// Watch out for + in the string
		String cleanCDot = cHGVS.replace("+", "\\+");
		
		// Major issues with the reported cHGVS
		//Tempus is changing c.1442_1443delGCinsCT to c.1442_1443delinsCT in the reports from what is in the vcf annotations
		if (cleanCDot.contains("delins")) {
			cleanCDot = cleanCDot.substring(0, cleanCDot.indexOf("delins")+3);
			//Misc.printErrAndExit(cHGVS+" -> "+cleanCDot);
		}
		
		//Tempus is changing c.10560_10561insAGTGGCGGC to c.10560_10561ins(9), why? Are they deliberately blocking the conversion to genomic coordinates
		else if (cleanCDot.contains("ins(")) {
			cleanCDot = cleanCDot.substring(0, cleanCDot.indexOf("ins(")+3);
			//Misc.printErrAndExit(cHGVS+" -> "+cleanCDot);
		}
		
		//Why are all of the AF's in the vcfs 0.167 ?  Yet in the json reports they are something else. Again, is Tempus deliberated obscuring these data?
		
		Pattern cdotPat = Pattern.compile(".*"+ cleanCDot+ ".*");
		Pattern tranPat = Pattern.compile(".*"+ transcript+ ".*");
		
		ArrayList<String> cDotTransMatches = new ArrayList<String>();
		if (vcfLinesFirst != null) {
			for (String line: vcfLinesFirst) {
				if (cdotPat.matcher(line).matches() && tranPat.matcher(line).matches()) {
					for (Pattern gene: genePats) {
						if (gene.matcher(line).matches()) {
							//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	...
							//   0       1   2   3   4
							String[] f = Misc.TAB.split(line);
							chromosome = f[0];
							pos = f[1];
							ref = f[3];
							alt = f[4];
							return true;
						}
					}
				}
			}
		}

		if (vcfLinesSecond != null) {
			for (String line: vcfLinesSecond) {
				if (cdotPat.matcher(line).matches() && tranPat.matcher(line).matches()) {
					for (Pattern gene: genePats) {
						if (gene.matcher(line).matches()) {
							//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	...
							//   0       1   2   3   4
							String[] f = Misc.TAB.split(line);
							chromosome = f[0];
							pos = f[1];
							ref = f[3];
							alt = f[4];
							return true;
						}
					}
				}
			}
		}
		
		//I really don't want to do this. Tempus is using some non existent gene aliases. E.g. LOC100129520 vcf -> TEX13C json report
		//So if only one match based on cDot and Transcript, then assume those coordinates.
		if (cDotTransMatches.size()==1) {
			String[] f = Misc.TAB.split(cDotTransMatches.get(0));
			chromosome = f[0];
			pos = f[1];
			ref = f[3];
			alt = f[4];
			IO.el("\nWARNING: Assigning genomic coordinates from just the cDot and partial transcript ID in: \n"+ this.toString());
			return true;
		}

		IO.el("\nWARNING: Failed to find genomic coordinate info, skippping:\n"+toString());
		return false;
		
	}
	
	public void setSnvVariantType() {
		//check variant type, Tempus groups all short vars as SNV
		if (variantType.equals("SNV") && ref!= null && alt!=null) {
			if (ref.length()< alt.length()) variantType = "INS";
			else if (ref.length()> alt.length()) variantType = "DEL";
		}
	}
	
	private void incrementCounters(TempusJson2Vcf t) {
		//these only get added if not null
		TempusJson2Vcf.add(gene, t.genes);
		TempusJson2Vcf.add(referenceGenome, t.referenceGenome);
		TempusJson2Vcf.add(variantType, t.variantType);
		TempusJson2Vcf.add(variantDescription, t.variantDescription);
		//histograms
		if (allelicFraction !=null) t.tumorAF.count(Double.parseDouble(allelicFraction));
		if (coverage !=null) t.tumorDP.count(Double.parseDouble(coverage));
	}
	
	/** Returns CHROM POS ID REF ALT QUAL FILTER INFO (EG FE ST PE DP AF)
	 * @throws IOException */
	public String toVcf(int id) throws IOException{
		if (gene5 != null && gene3 !=null && tempusJson2Vcf.getCnvGeneNameBed() != null) return fusionVariantToVcf(id);
		if (chromosome==null || pos==null || ref==null || alt==null) return null;
		StringBuilder sb = new StringBuilder();
		//CHROM
		sb.append(chromosome); sb.append("\t");
		//POS
		sb.append(pos); sb.append("\tTempus_");
		//ID
		sb.append(id);
		sb.append("\t");
		//REF
		sb.append(ref); sb.append("\t");
		//ALT
		sb.append(alt); sb.append("\t");
		//QUAL
		sb.append(".\t");
		//FILTER
		sb.append(".\t");
		//INFO (EG CL FE DP AF)
		sb.append("EG=");
		sb.append(gene);
		sb.append(";CL=");
		sb.append(variantSource);
		if (variantDescription != null) {
			sb.append(";FE=");
			sb.append(variantDescription.replaceAll(" ", "_"));
		}
		if (coverage !=null) {
			sb.append(";DP=");
			sb.append(coverage);
		}
		if (allelicFraction !=null) {
			double ap = Double.parseDouble(allelicFraction);
			double af = ap/100.0;
			sb.append(";AF=");
			sb.append(Num.formatNumber(af, 2));
		}
		//cnv?
		if (end !=0) {
			sb.append(";IMPRECISE;SVTYPE=CNV");
			if (copyNumber != null) {
				sb.append(";CN=");
				sb.append(copyNumber);
			}
			sb.append(";END=");
			sb.append(end);
			sb.append(";SVLEN=");
			sb.append(svLen);
		}
		return sb.toString();
	}
	
	private String fusionVariantToVcf(int id) throws IOException {
		//fetch bed coordinates
		Bed bed5 = tempusJson2Vcf.getCnvGeneNameBed().get(gene5);
		Bed bed3 = tempusJson2Vcf.getCnvGeneNameBed().get(gene3);
		if (bed5 == null || bed3 == null) throw new IOException("Failed to find fusion gene info in the lookup hash for "+ gene5 +" and/or "+gene3);
		
		//create common info 
		StringBuilder sb = new StringBuilder();
		sb.append("EG="); sb.append(gene5); sb.append(","); sb.append(gene3);
		sb.append(";CL="); sb.append(variantSource);
		sb.append(";FE="); sb.append(variantDescription.replaceAll(" ", "_"));
		//don't do this, lots of odd chars not permitted in INFO
		//if (geneDescription != null) {
			//sb.append(";DESC="); sb.append(geneDescription.replaceAll(" ", "_"));
		//}
		sb.append(";IMPRECISE;SVTYPE=BND");
		String commonInfo = sb.toString();
		
		String id5 = "Tempus_"+id+"_5";
		String id3 = "Tempus_"+id+"_3";
		
		String vcf5 = fetchFusionVcf(bed5, commonInfo, id5, id3);
		String vcf3 = fetchFusionVcf(bed3, commonInfo, id3, id5);
		
		return vcf5+"\n"+vcf3;
	}

	private String fetchFusionVcf(Bed bed, String commonInfo, String idThis, String idMate) {
		StringBuilder sb = new StringBuilder();
		//CHROM
		sb.append(bed.getChromosome()); sb.append("\t");
		//POS, bed is zero based
		int pos = bed.getStart()+1;
		sb.append(pos); sb.append("\t");
		//ID
		sb.append(idThis);
		sb.append("\t");
		//REF
		ReferenceSequence rs = tempusJson2Vcf.getFasta().getSubsequenceAt(bed.getChromosome(), pos, pos);
		sb.append(new String(rs.getBases())); sb.append("\t");
		//ALT
		sb.append("<BND>"); sb.append("\t");
		//QUAL and FILTER
		sb.append(".\t.\t");

		//INFO
		sb.append(commonInfo);
		sb.append(";END=");
		sb.append(bed.getStop());
		sb.append(";MATEID=");
		sb.append(idMate);
		return sb.toString();
	}

	/**Returns all non null values*/
	public String toString(){
		StringBuilder sb = new StringBuilder();
		if (accessionId!=null) sb.append("accessionId\t"+ accessionId+ "\n");
		if (variantSource != null) sb.append("variantSource\t"+variantSource+ "\n");
		if (gene != null) sb.append("gene\t"+gene+ "\n");
		if (transcript != null) sb.append("transcript\t"+transcript+ "\n"); 
		if (mutationEffect != null) sb.append("mutationEffect\t"+mutationEffect+ "\n"); 
		if (cHGVS != null) sb.append("cHGVS\t"+cHGVS+ "\n"); 
		if (pHGVS != null) sb.append("pHGVS\t"+pHGVS+ "\n"); 
		if (pFullHGVS != null) sb.append("pFullHGVS\t"+pFullHGVS+ "\n"); 
		if (referenceGenome != null) sb.append("referenceGenome\t"+referenceGenome+ "\n"); 
		if (variantType != null) sb.append("variantType\t"+variantType+ "\n"); 
		if (variantDescription != null) sb.append("variantDescription\t"+variantDescription+ "\n"); 
		if (nucleotideAlteration != null) sb.append("nucleotideAlteration\t"+nucleotideAlteration+ "\n"); 
		if (chromosome != null) sb.append("chromosome\t"+chromosome+ "\n"); 
		if (pos != null) sb.append("pos\t"+pos+ "\n"); 
		if (ref != null) sb.append("ref\t"+ref+ "\n"); 
		if (alt != null) sb.append("alt\t"+alt+ "\n"); 
		if (allelicFraction != null) sb.append("allelicFraction\t"+allelicFraction+ "\n"); 
		if (coverage != null) sb.append("coverage\t"+coverage+ "\n"); 
		if (copyNumber != null) sb.append("copyNumber\t"+copyNumber+ "\n"); 
		if (clinicalSignificance != null) sb.append("clinicalSignificance\t"+clinicalSignificance+ "\n"); 
		if (disease != null) sb.append("disease\t"+disease+ "\n"); 
		if (structuralVariant != null) sb.append("structuralVariant\t"+structuralVariant+ "\n"); 
		if (fusionType != null) sb.append("fusionType\t"+fusionType+ "\n"); 
		if (gene3 != null) sb.append("gene3\t"+gene3+ "\n"); 
		if (gene5 != null) sb.append("gene5\t"+gene5+ "\n"); 
		return sb.toString();
	}

	public void addCnvInfo(HashMap<String, Bed> cnvGeneNameBed, IndexedFastaSequenceFile fasta) throws IOException {
		Bed b = cnvGeneNameBed.get(gene);
		if (b == null) {
			//ugg Tempus!
			if (gene.equals("PD-L1")) b = cnvGeneNameBed.get("CD274");
			else if (gene.equals("PD-L2")) b = cnvGeneNameBed.get("PDCD1LG2");
			if (b == null) throw new IOException("Failed to find gene info in CNV lookup hash for "+ toString());
		}
		chromosome = b.getChromosome();
		start = b.getStart();
		pos = new Integer(start).toString();
		end = b.getStop();
		alt = "<CNV>";
		svLen = b.getLength();
		if ((variantType!=null && variantType.contains("deletion")) || (variantDescription!=null && variantDescription.toLowerCase().contains("loss"))) {
			svLen = -1*svLen;
			variantDescription = "CNV_LOSS";
		}
		else variantDescription = "CNV_GAIN";
		//find ref
		ReferenceSequence rs = fasta.getSubsequenceAt(chromosome, start, start);
		ref = new String(rs.getBases());
	}

	public String getVariantSource() {
		return variantSource;
	}

	public String getGene() {
		return gene;
	}

	public String getMutationEffect() {
		return mutationEffect;
	}

	public String getTranscript() {
		return transcript;
	}

	public String getcHGVS() {
		return cHGVS;
	}

	public String getpHGVS() {
		return pHGVS;
	}

	public String getpFullHGVS() {
		return pFullHGVS;
	}

	public String getReferenceGenome() {
		return referenceGenome;
	}

	public String getVariantType() {
		return variantType;
	}

	public String getVariantDescription() {
		return variantDescription;
	}

	public String getNucleotideAlteration() {
		return nucleotideAlteration;
	}

	public String getChromosome() {
		return chromosome;
	}

	public String getPos() {
		return pos;
	}

	public String getRef() {
		return ref;
	}

	public String getAlt() {
		return alt;
	}

	public String getAllelicFraction() {
		return allelicFraction;
	}

	public String getCoverage() {
		return coverage;
	}

	public String getCopyNumber() {
		return copyNumber;
	}

	public String getClinicalSignificance() {
		return clinicalSignificance;
	}

	public String getDisease() {
		return disease;
	}

	public String getStructuralVariant() {
		return structuralVariant;
	}

	public String getFusionType() {
		return fusionType;
	}

	public String getGene3() {
		return gene3;
	}

	public String getGene5() {
		return gene5;
	}


	public void setAccessionId(String accessionId) {
		this.accessionId = accessionId;
		
	}
}