package edu.utah.seq.vcf.json;

import java.io.IOException;
import java.util.HashMap;
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
	
	//for CNV
	private int start;
	private int end;
	private int svLen; //will be negative if deletion
		
	/**Object to represent a Tempus variant, SNV/ INDEL/ CNV/ Fusion*/
	public TempusVariant(String variantSource, String geneName, JSONObject object, TempusJson2Vcf tempusJson2Vcf) throws JSONException {
		this.variantSource = variantSource;
		
		//attempt to parse gene, if null pull from constructor
		gene = Json.forceGetString(object, "gene");
		if (gene == null) gene = geneName;

		mutationEffect = Json.forceGetString(object, "mutationEffect");
		transcript = Json.forceGetString(object, "transcript");
		cHGVS = Json.forceGetString(object, "HGVS.c");
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

		structuralVariant = Json.forceGetString(object, "structuralVariant");
		fusionType = Json.forceGetString(object, "fusionType");
		gene3 = Json.forceGetString(object, "gene3");
		gene5 = Json.forceGetString(object, "gene5");
		
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
		
		//increment counters
		incrementCounters(tempusJson2Vcf);
		
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
	
	/** Returns CHROM POS ID REF ALT QUAL FILTER INFO (EG FE ST PE DP AF)*/
	public String toVcf(int id){
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
	
	/**Returns all non null values*/
	public String toString(){
		StringBuilder sb = new StringBuilder();
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
		if (b == null) throw new IOException("Failed to find gene info in CNV lookup hash for "+ toString());
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
}