package edu.utah.seq.vcf.json.tempusv3;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import util.bio.annotation.Bed;
import util.gen.Json;
import util.gen.Misc;

public class TempusV3GenomicVariants {
	
	private String tumorMutationalBurden = null;
	private String tumorMutationBurdenPercentile = null;
	private String bloodTumorMutationalBurden = null;
	private String msiStatus = null;
	private int numberSnvIndels = 0;
	private int numberCnvs = 0;
	private int numberFusions = 0;
	private ArrayList<TempusV3Variant> variants = new ArrayList<TempusV3Variant>();
	private TempusV3Json2Vcf tempusJson2Vcf = null;
	private String[] diffExpressedGenes = null;
	private TempusV3Ihc ihc = null;
	private String ctDnaFraction = null;
	
	public TempusV3GenomicVariants(JSONObject object, TempusV3Json2Vcf tempusJson2Vcf) throws JSONException, IOException {	
		this.tempusJson2Vcf = tempusJson2Vcf;
		
		//msi
		JSONObject msiObj = object.getJSONObject("microsatelliteInstability");
		msiStatus = Json.forceGetString(msiObj, "status");

		JSONObject genomicVariants = object.getJSONObject("genomicVariants");
		
		//tumorMutationBurden
		JSONObject tmb = genomicVariants.getJSONObject("tumorMutationalBurden");
		tumorMutationalBurden = Json.forceGetString(tmb, "tmb");
		tumorMutationBurdenPercentile = Json.forceGetString(tmb, "tmbPercentile");
		
		//bloodTumorMutationalBurden
		JSONObject bloodTbm = genomicVariants.getJSONObject("bloodTumorMutationalBurden");
		bloodTumorMutationalBurden = Json.forceGetString(bloodTbm, "btmb");
		
		//variants
		tempusJson2Vcf.setWorkingNumPotentiallyActionable( loadVariants("potentiallyActionable", genomicVariants));
		tempusJson2Vcf.setWorkingNumBiologicallyRelevant( loadVariants("biologicallyRelevant", genomicVariants));
		tempusJson2Vcf.setWorkingNumPathogenic( loadVariants("pathogenic", genomicVariants));
		tempusJson2Vcf.setWorkingNumLikelyPathogenic( loadVariants("likelyPathogenic", genomicVariants));
		tempusJson2Vcf.setWorkingNumRiskAllele( loadVariants("riskAllele", genomicVariants));
		tempusJson2Vcf.setWorkingNumUnknownSignificance( loadVariants("unknownSignificance", genomicVariants));

		//cnvs processing, gene fusions are handled at the time of vcf generation
		addCoordinatesToCNVs(tempusJson2Vcf.getCnvGeneNameBed(), tempusJson2Vcf.getFasta());
		
		//variant processing
		addCoordinatesToVariants(tempusJson2Vcf.getWorkingSomVcfLines(), tempusJson2Vcf.getWorkingGermVcfLines(), tempusJson2Vcf.getFailedToFindCooridinates());

		//rna
		addRna(object);
		
		//IHC pdl1 and mmr, might be more
		ihc = new TempusV3Ihc(object, tempusJson2Vcf);
		
		//ctDNA fraction
		addCtDna(object);
		
	}
	
	private void addRna(JSONObject object) {
		if (object.isNull("rna")) return;
		JSONObject rnaObj = object.getJSONObject("rna");
		JSONArray rnaGenes = rnaObj.getJSONArray("variants");
		int num = rnaGenes.length();
		if (num !=0) {
			diffExpressedGenes = new String[num];
			for (int i=0; i<num; i++) {
				JSONObject gene = rnaGenes.getJSONObject(i);
				String geneName = gene.getString("gene");
				String mechanism = gene.getString("mechanism");
				diffExpressedGenes[i]= geneName+ "_"+ mechanism;
			}
		}
	}
	
	/*
	"ctDNA": {
    "disclaimer": "ctDNA tumor fraction is a quantitative measure of circulating tumor DNA.",
    "circulatingTumorCells": {
      "code": "C177519",
      "system": "NCI",
      "quantity": {
        "code": "246205007",
        "amount": 10.0,
        "system": "SNOMEDCT_US",
        "comparator": null,
        "unitOfMeasurement": {
          "code": "118582008",
          "system": "SNOMEDCT_US",
          "display": "%"
          
          Really care only about "amount"
	 */
	private void addCtDna(JSONObject object) {
		if (object.isNull("ctDNA")) return;
		JSONObject ctObj = object.getJSONObject("ctDNA");
		if (ctObj.isNull("circulatingTumorCells")) return;
		JSONObject ctCells = ctObj.getJSONObject("circulatingTumorCells");
		JSONObject quantity = ctCells.getJSONObject("quantity");
		ctDnaFraction = Json.forceGetString(quantity, "amount");
	}

	private int loadVariants(String variantSource, JSONObject genomicVariants) throws JSONException, IOException {
		int varsLoaded = 0;
		JSONObject pa = genomicVariants.getJSONObject(variantSource);
		JSONArray varTypes = pa.getJSONArray("variants");
		int numVars = varTypes.length();
		for (int i=0; i< numVars; i++) {
			JSONObject jo = varTypes.getJSONObject(i);
			String variantType = jo.getString("variantType"); //snvIndels, cnvs, and fusions
			JSONArray vd = jo.getJSONArray("variantDetails");
			int numVarDet = vd.length();
			varsLoaded += numVarDet;
			for (int j=0; j< numVarDet; j++) {
				JSONObject variant = vd.getJSONObject(j);
				variants.add( new TempusV3Variant(variantSource, variantType, variant, tempusJson2Vcf) );
			}
		}
		return varsLoaded;
		
	}

	private void addCoordinatesToVariants(String[] somVcfLines, String[] germVcfLines, ArrayList<TempusV3Variant> failedToFindCooridinates) throws IOException {
		if (somVcfLines == null && germVcfLines == null) return;
		
		for (TempusV3Variant tv: variants) {
			if (tv.getVariantType().contains("snv") == false) continue;
			
			numberSnvIndels++;
			
			// it's not clear which of these is somatic and which are inherited
			// potentiallyActionable, biologicallyRelevant, likelyPathogenic, riskAllele, unknownSignificance
			boolean found = false;
			//somatic 
			if (tv.getVariantSource().equals("potentiallyActionable") || tv.getVariantSource().equals("biologicallyRelevant") ) {
				found = tv.addGenomicCoordinateInfo(somVcfLines, germVcfLines);
				if (found == false) {
					failedToFindCooridinates.add(tv);
					tv.setAccessionId(tempusJson2Vcf.getWorkingOrder().getAccessionId());
				}
			}
			// others 
			else {
				found = tv.addGenomicCoordinateInfo(germVcfLines, somVcfLines);
				if (found == false) {
					failedToFindCooridinates.add(tv);
					tv.setAccessionId(tempusJson2Vcf.getWorkingOrder().getAccessionId());
				}
			}
			// set the snv type where relevant
			if (found) tv.setSnvVariantType();
		}
	}


	private void addCoordinatesToCNVs(HashMap<String, Bed> cnvGeneNameBed, IndexedFastaSequenceFile fasta) throws IOException {
		if (cnvGeneNameBed == null) return;
		for (TempusV3Variant tv: variants) {
			if (tv.getVariantType().equals("cnvs")) {
				tv.addCnvInfo(cnvGeneNameBed, fasta);
				numberCnvs++;
			}
			else if (tv.getVariantType().equals("fusions")) numberFusions++;
		}
	}

	/**Added to VCF Header*/
	public void addMetaData(LinkedHashMap<String, String> meta) {
		Misc.addConcatinatedValue (meta, "tempusTumorMutationalBurden", tumorMutationalBurden, "; ");
		Misc.addConcatinatedValue (meta, "tempusTumorMutationBurdenPercentile", tumorMutationBurdenPercentile, "; ");
		Misc.addConcatinatedValue (meta, "tempusMsiStatus", msiStatus, "; ");
		Misc.addConcatinatedValue (meta, "tempusBloodTumorMutationalBurden", bloodTumorMutationalBurden, "; ");
		//diff exp genes
		if (diffExpressedGenes!=null) Misc.addConcatinatedValue (meta, "tempusDiffExpGenes", Misc.stringArrayToString(diffExpressedGenes, ";"), "; ");
		//ihc
		Misc.addConcatinatedValue (meta, "PDL-1", ihc.getPdl1Result(), "; ");
		Misc.addConcatinatedValue (meta, "MMR", ihc.getMmrResult(), "; ");	
		//ctDNA?
		Misc.addConcatinatedValue (meta, "tempusCtDnaFraction", ctDnaFraction, "; ");
		
	}
	
	/**Added to spreadsheet, by json file*/
	public void addAttributes(LinkedHashMap<String, String> meta) {
		if (tumorMutationalBurden != null) meta.put("tumorMutationalBurden", tumorMutationalBurden);
		if (tumorMutationBurdenPercentile != null) meta.put("tumorMutationBurdenPercentile", tumorMutationBurdenPercentile);
		if (msiStatus != null) meta.put("msiStatus", msiStatus);
		if (bloodTumorMutationalBurden != null) meta.put("bloodTumorMutationalBurden", bloodTumorMutationalBurden);
		meta.put("numberReportedSnvIndels", ""+numberSnvIndels);
		meta.put("numberReportedCnvs", ""+numberCnvs);
		meta.put("numberReportedFusions", ""+numberFusions);
		//ihc
		if (ihc.getPdl1Result()!= null) {
			meta.put("PDL-1Result", ihc.getPdl1Result());
			if (ihc.getPdl1CombinedPositiveScore()!=null) meta.put("PDL-1CombinedPositiveScore", ihc.getPdl1CombinedPositiveScore());
			if (ihc.getPdl1TumorProportionScore()!=null) meta.put("PDL-1TumorProportionScore", ihc.getPdl1TumorProportionScore());
			if (ihc.getPdl1PercentImmuneCellStaining()!=null) meta.put("PDL-1PercentImmuneCellStaining", ihc.getPdl1PercentImmuneCellStaining());
			if (ihc.getPdl1PercentTumorCellStaining()!=null) meta.put("PDL-1PercentTumorCellStaining", ihc.getPdl1PercentTumorCellStaining());
		}
		if (ihc.getMmrResult()!=null) meta.put("MMRResult", ihc.getMmrResult());
		//diff exp genes
		if (diffExpressedGenes != null) meta.put("DiffExpGenes", Misc.stringArrayToString(diffExpressedGenes, ","));
		//ctDNA?
		if (ctDnaFraction != null) meta.put("ctDnaFraction", ctDnaFraction);
	}

	public String getTumorMutationalBurden() {
		return tumorMutationalBurden;
	}

	public String getTumorMutationBurdenPercentile() {
		return tumorMutationBurdenPercentile;
	}

	public String getMsiStatus() {
		return msiStatus;
	}

	public ArrayList<TempusV3Variant> getVariants() {
		return variants;
	}

	public TempusV3Ihc getIhc() {
		return ihc;
	}


}
