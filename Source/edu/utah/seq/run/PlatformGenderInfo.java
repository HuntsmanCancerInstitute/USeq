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
		if (fields.length == 7 && (clinInfoFile.endsWith("_M.json") || clinInfoFile.endsWith("_F.json"))) {
			if (fields[3].equals("deid")) {
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
			//New Avatar Format, AvatarID_NE_TE_TT_platform_submittingGroup_Gender.json
			//                      0      1  2  3     4         5            6
			//                   A028564_A61031_A60999_NA_IDTv1_GU_M.json 
			//                       0      1      2    3   4   5  6 
			else {
				platform = "Avatar";
				panel = fields[4];
				gender = fields[6].substring(0, 1);
				parsed = true;
			}
		}
		//Old Avatar format, AvatarID_Panel_Gender.json 
		//ZBG7MKFOZ1_IDT.V1_M.json
		//YLD7EH98AF_NIM.V1_M.json
		else if (fields.length == 3 && (clinInfoFile.endsWith("_M.json") || clinInfoFile.endsWith("_F.json"))) {
			platform = "Avatar";
			panel = fields[1];
			gender = fields[2].substring(0, 1);
			parsed = true;
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
