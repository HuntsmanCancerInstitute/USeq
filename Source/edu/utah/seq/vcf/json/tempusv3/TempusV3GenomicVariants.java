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

public class TempusV3GenomicVariants {
	
	private String tumorMutationalBurden = null;
	private String tumorMutationBurdenPercentile = null;
	private String msiStatus = null;
	private int numberSnvIndels = 0;
	private int numberCnvs = 0;
	private ArrayList<TempusV3Variant> variants = new ArrayList<TempusV3Variant>();
	private TempusV3Json2Vcf tempusJson2Vcf = null;
	
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
		
		//variants
		tempusJson2Vcf.setWorkingNumPotentiallyActionable( loadVariants("potentiallyActionable", genomicVariants));
		tempusJson2Vcf.setWorkingNumBiologicallyRelevant( loadVariants("biologicallyRelevant", genomicVariants));
		tempusJson2Vcf.setWorkingNumLikelyPathogenic( loadVariants("likelyPathogenic", genomicVariants));
		tempusJson2Vcf.setWorkingNumRiskAllele( loadVariants("riskAllele", genomicVariants));
		tempusJson2Vcf.setWorkingNumUnknownSignificance( loadVariants("unknownSignificance", genomicVariants));

		//cnvs processing
		addCoordinatesToCNVs(tempusJson2Vcf.getCnvGeneNameBed(), tempusJson2Vcf.getFasta());
		
		//variant processing
		addCoordinatesToVariants(tempusJson2Vcf.getWorkingSomVcfLines(), tempusJson2Vcf.getWorkingGermVcfLines(), tempusJson2Vcf.getFailedToFindCooridinates());
	}
	
	private int loadVariants(String variantSource, JSONObject genomicVariants) throws JSONException, IOException {
		int varsLoaded = 0;
		JSONObject pa = genomicVariants.getJSONObject(variantSource);
		JSONArray varTypes = pa.getJSONArray("variants");
		int numVars = varTypes.length();
		for (int i=0; i< numVars; i++) {
			JSONObject jo = varTypes.getJSONObject(i);
			String variantType = jo.getString("variantType"); //snvIndels cnvs fusions? not seeing any
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
		}
	}

	/**Added to VCF Header*/
	public void addMetaData(LinkedHashMap<String, String> meta) {
		if (tumorMutationalBurden != null) meta.put("tempusTumorMutationalBurden", tumorMutationalBurden.toString());
		if (tumorMutationBurdenPercentile != null) meta.put("tempusTumorMutationBurdenPercentile", tumorMutationBurdenPercentile.toString());
		if (msiStatus != null) meta.put("tempusMsiStatus", msiStatus);
	}
	
	/**Added to spreadsheet*/
	public void addAttributes(LinkedHashMap<String, String> meta) {
		if (tumorMutationalBurden != null) meta.put("tumorMutationalBurden", tumorMutationalBurden.toString());
		if (tumorMutationBurdenPercentile != null) meta.put("tumorMutationBurdenPercentile", tumorMutationBurdenPercentile.toString());
		if (msiStatus != null) meta.put("msiStatus", msiStatus);
		meta.put("numberSnvIndels", ""+numberSnvIndels);
		meta.put("numberCnvs", ""+numberCnvs);
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


}
