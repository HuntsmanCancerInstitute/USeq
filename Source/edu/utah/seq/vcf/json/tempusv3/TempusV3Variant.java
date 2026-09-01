package edu.utah.seq.vcf.json.tempusv3;

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

public class TempusV3Variant{

	//This is a composite of all of the variant fields in the 10 v3.30 json reports we've received. Still missing fusion info.
	//New
	private String variantType = null; //snvIndels, cnvs, isoforms
	private String variantSource  = null; //potentiallyActionable, biologicallyRelevant, likelyPathogenic, riskAllele, unknownSignificance
	private String allelicFraction = null;
	private String cVar = null;
	private String pVar = null;
	private String gene = null;
	private String fusedGene = null;
	private String gene3 = null;
	private String gene5 = null;
	private String geneDescription = null;
	private String genomicSourceClass = null; //somatic, genomic variant,
	private String transcriptId = null; //not in v3.4
	private String variantDescription = null; //Copy number gain, Copy number loss, Frameshift
	private String reportChr = null; //in v3.4+
	private String reportPos = null; //in v3.4+

	//Old
	private String chromosome = null;
	private String pos = null;
	private String ref = null;
	private String alt = null;

	private TempusV3Json2Vcf tempusJson2Vcf = null;
	private String accessionId = null;

	//for CNV
	private int start;
	private int end;
	private int svLen; //will be negative if deletion

	/**Object to represent a Tempus variant, SNV/ INDEL/ CNV/ Fusion
	 * @throws IOException */
	public TempusV3Variant(String variantSource, String variantType, JSONObject object, TempusV3Json2Vcf tempusJson2Vcf) throws JSONException, IOException {
		this.variantSource = variantSource;
		this.tempusJson2Vcf = tempusJson2Vcf;
		this.variantType = variantType;		

		allelicFraction = Json.forceGetString(object, "allelicFraction");
		cVar = Json.forceGetString(object, "cVar");
		fusedGene  = Json.forceGetString(object, "fusedGene");
		gene = Json.forceGetString(object, "gene");
		gene3 = Json.forceGetString(object, "gene3");
		gene5 = Json.forceGetString(object, "gene5");
		geneDescription = Json.forceGetString(object, "geneDescription");
		genomicSourceClass  = Json.forceGetString(object, "genomicSourceClass");
		pVar = Json.forceGetString(object, "pVar");
		transcriptId  = Json.forceGetString(object, "transcriptId");
		variantDescription = Json.forceGetString(object, "variantDescription");
		
		//v3.4 contains chr and pos, hg19
		attemptToParseChrPos(object);

		//increment counters
		incrementCounters(tempusJson2Vcf);
	}

	/*
	 * "chromosome": {
                              "code": "48000-4",
                              "system": "LNC",
                              "value": {
                                   "code": "LA21267-2",
                                   "system": "LNC",
                                   "display": "14"
                              }
                         },

                              "position": {
                              "code": "81254-5",
                              "system": "LNC",
                              "display": "95590735"
                         }
	 */
	private void attemptToParseChrPos(JSONObject object) {
		//attempt to parse chromosome
		String chr = null;
		if (object.has("chromosome")) {
			JSONObject c = TempusV3Json2Vcf.getJsonObjectOrNull("chromosome", object);
			if (c!=null && c.has("value")) {
				JSONObject v= TempusV3Json2Vcf.getJsonObjectOrNull("value",c);
				if (v!=null && v.has("display")) chr = v.getString("display");
			}
		}
		String pos = null;
		if (object.has("position")) {
			JSONObject p= TempusV3Json2Vcf.getJsonObjectOrNull("position",object);
			if (p!=null && p.has("display")) pos = p.getString("display");
		}

		if (chr !=null && pos != null) {
			reportChr = chr;
			reportPos = pos;
		}

	}


	public boolean addGenomicCoordinateInfo(String[] vcfLinesFirst, String[] vcfLinesSecond) throws IOException {
		if (gene == null || ( cVar == null && reportPos == null)) {
			throw new IOException("\nMissing SNV Info for genomic coordinate parsing\n"+toString());
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

		//any cVar?
		Pattern cdotPat = null;
		
		if (cVar != null) {
			// Watch out for + in the string
			String cleanCDot = cVar.replace("+", "\\+");
			// Major issues with the reported cHGVS
			//Tempus is changing c.1442_1443delGCinsCT to c.1442_1443delinsCT in the reports from what is in the vcf annotations
			if (cleanCDot.contains("delins")) {
				cleanCDot = cleanCDot.substring(0, cleanCDot.indexOf("delins")+3);
			}
			//Tempus is changing c.10560_10561insAGTGGCGGC to c.10560_10561ins(9), why? Are they deliberately blocking the conversion to genomic coordinates
			else if (cleanCDot.contains("ins(")) {
				cleanCDot = cleanCDot.substring(0, cleanCDot.indexOf("ins(")+3);
			}
			//Why are all of the AF's in the vcfs 0.167 ?  Yet in the json reports they are something else. Again, is Tempus deliberated obscuring these data?
			cdotPat = Pattern.compile(".*"+ cleanCDot+ ".*");
		}
		
		Pattern tranPat = null;
		Pattern posPat = null;
		if (transcriptId !=null) tranPat = Pattern.compile(".*"+ transcriptId+ ".*");
		if (reportPos !=null) posPat = Pattern.compile(".+"+ reportPos+ ".+");
		
		//Tempus is now using ensembl transcript IDs in the json to refseq in the vcf!!!!! So the transcript matching won't work, e.g. ENST00000288135.5 in the report to NM_000222.2 in the vcf 
		// The only way to get ref and alt is from the vcf, these are not in the json
		if (vcfLinesFirst != null) {
			for (String line: vcfLinesFirst) {
				boolean itMatches = vcfMatches(line, cdotPat, tranPat, posPat, genePats);
				if (itMatches) {
					String[] f = Misc.TAB.split(line);
					chromosome = f[0];
					pos = f[1];
					ref = f[3];
					alt = f[4];
					return true;
				}
			}
		}
	
		if (vcfLinesSecond != null) {
			for (String line: vcfLinesSecond) {
				boolean itMatches = vcfMatches(line, cdotPat, tranPat, posPat, genePats);
				if (itMatches) {
					String[] f = Misc.TAB.split(line);
					chromosome = f[0];
					pos = f[1];
					ref = f[3];
					alt = f[4];
					return true;
				}
			}
		}
		IO.el("\nWARNING: Failed to find genomic coordinate info, skippping:\n"+toString());
		return false;
	}

	private boolean vcfMatches(String line, Pattern cdotPat, Pattern tranPat, Pattern posPat, Pattern[] genePats) {
		/*IO.pl("\nTesting\t"+line);
		IO.pl("\tcDotPat\t"+cdotPat);
		IO.pl("\ttranPat\t"+tranPat);
		IO.pl("\tposPat\t"+posPat);
		IO.pl("\tgenePats\t"+genePats[0]);*/
		
		//one of these must match
		boolean ok = false;
		if (cdotPat!=null && cdotPat.matcher(line).matches()) ok= true;
		if (tranPat!=null && tranPat.matcher(line).matches()) ok= true;
		if (posPat!=null && posPat.matcher(line).matches()) ok= true;
		if (ok == false) return false;
		
		//likewise one of the genes must match
		for (Pattern gene: genePats) if (gene.matcher(line).matches()) return true;
		return false;
	}

	public void setSnvVariantType() {
		//check variant type, Tempus groups all short vars as SNV
		if (variantType.equals("snvIndels") && ref!= null && alt!=null) {
			if (ref.length()< alt.length()) variantType = "INS";
			else if (ref.length()> alt.length()) variantType = "DEL";
		}
	}

	private void incrementCounters(TempusV3Json2Vcf t) {
		//these only get added if not null
		TempusV3Json2Vcf.add(gene, t.genes);
		TempusV3Json2Vcf.add(genomicSourceClass, t.genomicSourceClass);
		TempusV3Json2Vcf.add(variantType, t.variantType);
		TempusV3Json2Vcf.add(variantDescription, t.variantDescription);
		//histograms
		if (allelicFraction !=null) t.tumorAF.count(Double.parseDouble(allelicFraction));
	}

	/** Returns CHROM POS ID REF ALT QUAL FILTER INFO (EG FE ST PE DP AF)
	 * @throws IOException */
	public String toVcf(int id) throws IOException{

		if (variantType.equals("fusions")) return fusionVariantToVcf(id);
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
		if (allelicFraction !=null) {
			double ap = Double.parseDouble(allelicFraction);
			double af = ap/100.0;
			sb.append(";AF=");
			sb.append(Num.formatNumber(af, 2));
		}
		//cnv?
		if (end !=0) {
			sb.append(";IMPRECISE;SVTYPE=CNV");
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
		if (bed5 == null && bed3 == null) {
			throw new IOException("Failed to find fusion gene info in the lookup hash for "+ gene5 +" and "+gene3);
		}
		// some genes have internal alterations, like EGFRvIII, and Tempus doesn't report the gene3
		if (bed5 == null) bed5 = bed3;
		else if (bed3 == null) bed3 = bed5;

		//create common info 
		StringBuilder sb = new StringBuilder();
		sb.append("EG="); sb.append(gene5); sb.append(","); sb.append(gene3);
		sb.append(";CL="); sb.append(variantSource);
		sb.append(";FE="); sb.append(variantDescription.replaceAll(" ", "_"));
		//don't do this, lots of odd chars not permitted in INFO
		//if (geneDescription != null) {
		//sb.append(";DESC="); sb.append(geneDescription.replaceAll(" ", "_"));
		//}
		if (fusedGene!=null) sb.append(";FG="+fusedGene);

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
		if (fusedGene != null) sb.append("fusedGene\t"+fusedGene+ "\n");
		if (transcriptId != null) sb.append("transcriptId\t"+transcriptId+ "\n"); 
		if (variantDescription != null) sb.append("variantDescription\t"+variantDescription+ "\n"); 
		if (cVar != null) sb.append("cVar\t"+cVar+ "\n"); 
		if (pVar != null) sb.append("pVar\t"+pVar+ "\n"); 
		if (variantType != null) sb.append("variantType\t"+variantType+ "\n"); 
		if (variantDescription != null) sb.append("variantDescription\t"+variantDescription+ "\n"); 
		if (geneDescription != null) sb.append("geneDescription\t"+geneDescription+ "\n"); 
		if (genomicSourceClass != null) sb.append("genomicSourceClass\t"+genomicSourceClass+ "\n"); 
		if (chromosome != null) sb.append("chromosome\t"+chromosome+ "\n"); 
		if (pos != null) sb.append("pos\t"+pos+ "\n"); 
		if (ref != null) sb.append("ref\t"+ref+ "\n"); 
		if (alt != null) sb.append("alt\t"+alt+ "\n"); 
		if (reportChr != null) sb.append("reportChr\t"+reportChr+ "\n"); 
		if (reportPos != null) sb.append("reportPos\t"+reportPos+ "\n"); 
		if (allelicFraction != null) sb.append("allelicFraction\t"+allelicFraction+ "\n"); 
		if (gene3 != null) sb.append("gene3\t"+gene3+ "\n"); 
		if (gene5 != null) sb.append("gene5\t"+gene5+ "\n"); 
		return sb.toString();
	}

	public void addCnvInfo(HashMap<String, Bed> cnvGeneNameBed, IndexedFastaSequenceFile fasta) throws IOException {
		Bed b = cnvGeneNameBed.get(gene);
		if (b == null) {
			//ugg Tempus! better to just add these genes to the bed file
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
		if (variantType!=null && variantType.contains("loss")) {
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

	public String getVariantType() {
		return variantType;
	}

	public String getVariantDescription() {
		return variantDescription;
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