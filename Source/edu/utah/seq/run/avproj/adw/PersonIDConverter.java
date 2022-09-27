package edu.utah.seq.run.avproj.adw;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;

public class PersonIDConverter {
	
	//fields
	private boolean parsed = true;
	private String header = null;
	private HashMap<String, Integer> headerKeyIndex = new HashMap<String, Integer>();
	private HashMap<String, String[]> personIdDataLine = new HashMap<String, String[]>();
	private HashMap<String, String> avatarIDGender = null;
	private static final Pattern genderPat = Pattern.compile("[MF]");
	private int lastNameIndex = -1;
	private int firstNameIndex = -1;
	private int dobIndex = -1;
	private int genderIndex = -1;
	private int mrnIndex = -1;
	private int avatarIdIndex = -1;
	private int hciIdIndex = -1;

	public PersonIDConverter( File csv) throws IOException {
		
		parseIt(csv);
		
		makeHeaderLookup();
		
		findIndexesForSMM();

	}
	
	private void findIndexesForSMM() throws IOException {
		StringBuilder sb = new StringBuilder();
		lastNameIndex = findIndex("lastName", sb);
		firstNameIndex = findIndex("firstName", sb);
		dobIndex = findIndex("birthDate", sb);
		genderIndex = findIndex("gender", sb);
		mrnIndex = findIndex("mrn", sb);
		avatarIdIndex = findIndex("avatarId1", sb);
		hciIdIndex = findIndex("hciPersonId", sb);
		if (sb.length() != 0) throw new IOException("Problems finding particular columns in the person ID conversion PHI file.\n"+sb);
	}
	
	private int findIndex(String name, StringBuilder sb) {
		if (headerKeyIndex.containsKey(name)) return headerKeyIndex.get(name);
		else {
			sb.append("\tFailed to find "+name+" in "+header+"\n");
			return -1;
		}
	}

	/**lastName firstName dobMonth(1-12) dobDay(1-31) dobYear(1900-2050) gender(M|F) mrn */
	public String fetchSubjectMatchMakerLine(String[] l) {
		//Textbox24,Textbox26,avatarId1,mrn,hciPersonId,firstName,middleName,lastName,birthDate,gender
		//HCI Person ID,MRN,A037950,12869785,00297846,KENNITH,M,JOHNSON,1/19/1947,M - mock entry
		//HCI Person ID,MRN,A034305,12117335,00862698,SUSAN,P,SMITH,6/31/1955,F - mock entry
		
		StringBuilder sb = new StringBuilder();
		
		sb.append(l[lastNameIndex]);
		sb.append("\t");
		
		sb.append(l[firstNameIndex]);
		sb.append("\t");

		String[] mdy = Misc.FORWARD_SLASH.split(l[dobIndex]);
		if (mdy.length == 3) {
			Integer dobMonth = Integer.parseInt(mdy[0]);
			sb.append(dobMonth.toString());
			sb.append("\t");
			Integer dobDay = Integer.parseInt(mdy[1]);
			sb.append(dobDay.toString());
			sb.append("\t");
			Integer dobYear = Integer.parseInt(mdy[2]);
			sb.append(dobYear.toString());
			sb.append("\t");
		}
		else sb.append("\t\t\t");
		
		sb.append(l[genderIndex]);
		sb.append("\t");
		
		sb.append(l[mrnIndex]);
		sb.append("\t\t");
		
		sb.append(l[avatarIdIndex]);
		sb.append(";");
		
		sb.append(l[hciIdIndex]);
		return sb.toString();
	}
	
	private void writeOutSubjectMatchMakerFile( File txtFile) throws IOException {
		PrintWriter out = new PrintWriter( new FileWriter(txtFile));
		out.println("#lastName\tfirstName\tdobMonth(1-12)\tdobDay(1-31)\tdobYear(1900-2050)\tgender(M|F)\tmrn\tavatarId;hciId");
		for (String[] pId: personIdDataLine.values()) out.println(fetchSubjectMatchMakerLine(pId));
		out.close();
	}


	private void makeGenderLookup() throws IOException {
		avatarIDGender = new HashMap<String, String>();
		for (String avatarID: personIdDataLine.keySet()) {
			String gender = personIdDataLine.get(avatarID)[genderIndex];
			if (genderPat.matcher(gender).matches() == false) throw new IOException("Failed to find a M or F for gender for "+avatarID);
			else avatarIDGender.put(avatarID, gender);
		}
	}


	private void makeHeaderLookup() {
		String[] f = Misc.COMMA.split(header);
		for (int i=0; i<f.length; i++) headerKeyIndex.put(f[i], i);
	}

	private void parseIt(File cvs) {
		String line = null;
		try {
			BufferedReader in = IO.fetchBufferedReader(cvs);
			while ((line = in.readLine()) !=null) {
				if (line.contains("hciPersonId")) {
					header = line.trim();
					makeHeaderLookup();
					int avaIdIndex = headerKeyIndex.get("avatarId1");
					//read in data lines
					while ((line = in.readLine()) !=null) {
						line = line.trim();
						if (line.length() == 0) continue;
						//split on ,
						String[] f =  Misc.COMMA.split(line);
						if (personIdDataLine.containsKey(f[avaIdIndex])) throw new IOException("Found a duplicate avatarId -> "+f[avaIdIndex]);
						else personIdDataLine.put(f[avaIdIndex], f);
					}
				}
			}
			in.close();
			parsed = true;
		} catch (IOException e) {
			IO.el("Failed to parse "+ line + " in "+cvs);
			e.printStackTrace();
		}
	}

	public static void main (String[] args) throws IOException {
		File test = new File ("/Users/u0028003/HCI/AvatarORIEN/AutoAvatar/ResourceFiles/19Sept2022_PatientPHI.csv");
		PersonIDConverter cml = new PersonIDConverter(test);
		IO.pl(cml.getHeaderKeyIndex());
		IO.pl(cml.getPersonIDDataLine());
		cml.writeOutSubjectMatchMakerFile(new File("/Users/u0028003/HCI/AvatarORIEN/AutoAvatar/ResourceFiles/smm_PatientPHI.txt"));
	}

	public boolean isParsed() {
		return parsed;
	}

	public String getHeader() {
		return header;
	}


	public HashMap<String, Integer> getHeaderKeyIndex() {
		return headerKeyIndex;
	}


	public HashMap<String, String[]> getPersonIDDataLine() {
		return personIdDataLine;
	}


	public HashMap<String, String> getAvatarIDGender() {
		return avatarIDGender;
	}

	
}
