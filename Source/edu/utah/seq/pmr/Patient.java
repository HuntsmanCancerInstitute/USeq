package edu.utah.seq.pmr;
import java.util.HashMap;

public class Patient {
	
	//fields
	private String coreId = null;
	private HashMap<String, Dataset> idDataSets = new HashMap<String, Dataset>();

	public Patient(String coreId) {
		this.coreId = coreId;
	}

	public void addDataFile(String[] keys) {
		//Patients   AA2mF6Vy   Avatar   A032049_SL419345_SL419548_SL420681     ClinicalReport    A032049_SL419345_SL419548_SL420681_IDTv1_SAR_F.json
		//   0          1         2                 3                                  4                5+
		String source = keys[2];
		String datasetId = keys[3];
		Dataset d = idDataSets.get(datasetId);
		if (d == null) {
			d = new Dataset(source, datasetId);
			idDataSets.put(datasetId, d);
		}
		StringBuilder sb = new StringBuilder(keys[4]);
		for (int i=5; i< keys.length; i++) {
			sb.append("/");
			sb.append(keys[i]);
		}
		d.getPartialPaths().add(sb.toString());
	}

	public String getCoreId() {
		return coreId;
	}

	public HashMap<String, Dataset> getIdDataSets() {
		return idDataSets;
	}


}
