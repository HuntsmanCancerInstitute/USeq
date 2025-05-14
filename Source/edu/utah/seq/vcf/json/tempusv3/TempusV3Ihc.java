package edu.utah.seq.vcf.json.tempusv3;

import java.io.IOException;
import org.json.JSONObject;
import util.gen.Json;


/**Seems like there are only two IHC tests, pdl1 and mmr.  More might follow. This will look for new ones and error out.*/
public class TempusV3Ihc {
	
	private String pdl1Result = null;	//"positive - xxx", "negative", "null"
	private String pdl1TumorProportionScore = null;
	private String pdl1CombinedPositiveScore = null;
	private String pdl1PercentTumorCellStaining = null;
	private String pdl1PercentImmuneCellStaining = null;
	
	private String mmrResult = null;	//"abnormal", "?"
	
	public TempusV3Ihc (JSONObject obj, TempusV3Json2Vcf parser) throws IOException {
		if (obj.isNull("ihc")) return;
		
		JSONObject ihcObj = obj.getJSONObject("ihc");
		
		//look for pdl1, mmr, and anything else
		boolean mmrFound = false;
		boolean pdlFound = false;
		boolean otherFound = false;
				
		for (String keyName: ihcObj.getNames(ihcObj)) {
			if (keyName.equals("pdl1")) pdlFound = true;
			else if (keyName.equals("mmr")) mmrFound = true;
			else otherFound = true;
		}
		if (otherFound) throw new IOException("ERROR: new IHC result found in "+parser.getWorkingJsonFile());
		if (mmrFound==false || pdlFound==false) throw new IOException("ERROR: failed to find mmr or pdl1 IHC results in "+parser.getWorkingJsonFile());
		
		if (ihcObj.isNull("pdl1") == false) {
			JSONObject pdl1Obj = ihcObj.getJSONObject("pdl1");
			if (pdl1Obj.has("interpretation")) {
				JSONObject interp = pdl1Obj.getJSONObject("interpretation");
				pdl1Result = Json.forceGetString(interp, "result");
				pdl1TumorProportionScore = Json.forceGetString(interp, "tumorProportionScore");
				pdl1CombinedPositiveScore = Json.forceGetString(interp, "combinedPositiveScore");
				pdl1PercentTumorCellStaining = Json.forceGetString(interp, "percentTumorCellStaining");
				pdl1PercentImmuneCellStaining = Json.forceGetString(interp, "percentImmuneCellStaining");
			}
		}
		
		if (ihcObj.isNull("mmr") == false) {
			JSONObject mmrObj = ihcObj.getJSONObject("mmr");
			mmrResult = Json.forceGetString(mmrObj, "interpretation");
		}
	}

	public String getPdl1Result() {
		return pdl1Result;
	}

	public String getPdl1TumorProportionScore() {
		return pdl1TumorProportionScore;
	}

	public String getPdl1CombinedPositiveScore() {
		return pdl1CombinedPositiveScore;
	}

	public String getPdl1PercentTumorCellStaining() {
		return pdl1PercentTumorCellStaining;
	}

	public String getPdl1PercentImmuneCellStaining() {
		return pdl1PercentImmuneCellStaining;
	}

	public String getMmrResult() {
		return mmrResult;
	}
	
	/*
"ihc": {
    "pdl1": {
      "interpretation": {
        "result": "positive - high",   or "negative",   or null
        "tumorProportionScore": "90%",
        "combinedPositiveScore": "92",
        "percentTumorCellStaining": null,
        "percentImmuneCellStaining": null
      }
    },
    "mmr": {
      "details": [
        {
          "result": "absent",
          "protein": "mlh1"
        },
        {
          "result": "present",
          "protein": "msh2"
        },
        {
          "result": "present",
          "protein": "msh6"
        },
        {
          "result": "absent",
          "protein": "pms2"
        }
      ],
      "interpretation": "abnormal"
    }
  },	 * */
}
