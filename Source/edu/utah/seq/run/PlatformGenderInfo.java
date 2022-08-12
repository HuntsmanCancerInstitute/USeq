package edu.utah.seq.run;

import util.gen.Misc;

public class PlatformGenderInfo {
	//fields
	private boolean parsed = false;
	private String platform = null;
	private String panel = null;
	private String gender = null;
	private String originalName = null;
	
	public PlatformGenderInfo(String clinInfoFile) {
		originalName = clinInfoFile;
		//TL-18-843E9B_XT.V1_2018-10-26_deid_Neeraj_Agarwal_M.json
		//     0         1        2       3     4      5    6
		String[] fields = Misc.UNDERSCORE.split(clinInfoFile);
		if (fields.length == 7) {
			if (fields[3].equals("deid") && fields[6].startsWith("M.") || fields[6].startsWith("F.")) {
				//platform
				if (fields[0].startsWith("TL-")) {
					platform = "Tempus";
					//panel
					panel = fields[1];
					//skip RNA only reports, RS.vxxx
					if (panel.toLowerCase().contains("rs.v")) return;
					//gender
					gender = fields[6].substring(0, 1);
					parsed = true;
				}
			}
		}
	}

	public boolean isParsed() {
		return parsed;
	}

	public String getPlatform() {
		return platform;
	}

	public String getPanel() {
		return panel;
	}

	public String getGender() {
		return gender;
	}

	public String getOriginalName() {
		return originalName;
	}
}
