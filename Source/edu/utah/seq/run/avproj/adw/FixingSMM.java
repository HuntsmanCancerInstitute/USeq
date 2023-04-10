package edu.utah.seq.run.avproj.adw;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import util.gen.IO;
import util.gen.Misc;

public class FixingSMM {
	
	File corrData = new File("/Users/u0028003/HCI/AvatarORIEN/AutoAvatar/MoveOld2AWS/AvatarDataWranglerDoOver/CorrectingSMMRegistry/GetPersonIDByMRN.csv");
	HashMap<String, String[]> avatarIdCorrData = null;
	HashMap<String, String[]> mrnCorrData = null;
	HashMap<String, ArrayList<String>> duplicateMrns = new HashMap<String, ArrayList<String>>();
	HashMap<String, ArrayList<String>> duplicateDoBs = new HashMap<String, ArrayList<String>>();
	
	File smmFileToFix = new File("/Users/u0028003/HCI/AvatarORIEN/AutoAvatar/MoveOld2AWS/AvatarDataWranglerDoOver/CorrectingSMMRegistry/oldRegistry_1676559934551_PHI.txt");
	File smmFileFixed = new File("/Users/u0028003/HCI/AvatarORIEN/AutoAvatar/MoveOld2AWS/AvatarDataWranglerDoOver/CorrectingSMMRegistry/oldRegistry_1676559934551_PHI.fixed.txt");

	public FixingSMM() throws IOException {
		parseCorrectData();
		walkSmmToFix();
		checkForDuplicates();
		IO.pl("\nDONE");
	}
	
	

	
	private void checkForDuplicates() {
		IO.pl("\nChecking for duplicate MRNs in fixed repo...");
		for (String mrn: duplicateMrns.keySet()) {
			ArrayList<String> al = duplicateMrns.get(mrn);
			if (al.size()>1) {
				for (String l: al) {
					IO.pl(l+"\txxxxxx");
				}
				IO.pl();
			}
		}
		
		IO.pl("\nChecking for duplicate DoBs in fixed repo...");
		for (String mrn: duplicateDoBs.keySet()) {
			ArrayList<String> al = duplicateDoBs.get(mrn);
			if (al.size()>1) {
				for (String l: al) {
					IO.pl(l+"\tyyyyyy");
				}
				IO.pl();
			}
		}
		
	}




	private void walkSmmToFix() throws IOException {
		
		BufferedReader in = IO.fetchBufferedReader(smmFileToFix);
		PrintWriter out = new PrintWriter( new FileWriter(smmFileFixed));
		String[] fields = null;
		String line = null;
		while ((line = in.readLine())!=null) {
			line = line.trim();
			if(line.length()==0) continue;
			if (line.startsWith("#LastName")) out.println(line);
			else {
				//#LastName	FirstName	DoBMonth(1-12)	DoBDay(1-31)	DoBYear(1900-2050)	Gender(M|F)	MRN	CoreId	OtherIds(AvatarId;HciId) - what it should be
				//#LastName	FirstName	DoBMonth(1-12)	DoBDay(1-31)	DoBYear(1900-2050)	Gender(M|F)	HciID	CoreId	OtherIds(AvatarId;MRN)	- what it is for Avatar cases
				//    0         1            2             3                  4                 5         6       7            8
				//so parse the OtherIds, pull the AvatarId, swap the HciID for the MRN, check with the corrData, write out corrected

				fields = Misc.TAB.split(line);
				
				//does it have a otherIds?
				if (fields.length == 8) {
					IO.pl("\nNo OtherIds, just saving -> "+line);
					out.println(line+"\t");
					addToHashes(fields);
				}
				else if (fields.length == 9) {
					//does it have a 'A022289;10400893'  AvatarId;MRN?
					String[] otherIds = Misc.SEMI_COLON.split(fields[8]);
					if (otherIds.length == 2 && otherIds[0].startsWith("A0")) {
						IO.pl("\nAvatarId found, swapping\nbad: "+line);
						String mrn = otherIds[1];
						String hciId = fields[6];
						
						//check against correct
						//avatarId1,mrn,hciPersonId,firstName,middleName,lastName,birthDate,gender
						//   0       1      2         
						String[] corrFields = avatarIdCorrData.get(otherIds[0]);
						if (mrn.equals(corrFields[1]) == false || hciId.equals(corrFields[2]) == false) throw new IOException("\nMrn or HciId mismatch\nSmm: "+line+"\nCorr: "+Misc.stringArrayToString(corrFields, "\t"));
						
						fields[6] = mrn;
						otherIds[1]=hciId;
						fields[8] = otherIds[0]+";"+otherIds[1];
						String fixedLine = Misc.stringArrayToString(fields, "\t");
						IO.pl("fix: "+fixedLine);
						out.println(fixedLine);
						addToHashes(fields);
					}
					else throw new IOException("\nOtherIds but no match "+line);
				}
				else throw new IOException("Odd smm line? "+line+" length: "+fields.length);
			}
			
			
		}
		in.close();
		out.close();
		
	}

	private void addToHashes(String[] f) {
		//#LastName	FirstName	DoBMonth(1-12)	DoBDay(1-31)	DoBYear(1900-2050)	Gender(M|F)	MRN	CoreId	OtherIds(AvatarId;HciId) - what it should be
		//    0         1            2             3                  4                 5        6     7            8
		String mrn = f[6];
		String dob = f[2]+"/"+f[3]+"/"+f[4];
		
		ArrayList<String> alMrns =  duplicateMrns.get(mrn);
		if (alMrns == null) {
			alMrns = new ArrayList<String>();
			duplicateMrns.put(mrn, alMrns);
		}
		alMrns.add(Misc.stringArrayToString(f, "\t"));
		
		ArrayList<String> alDoB =  duplicateDoBs.get(dob);
		if (alDoB == null) {
			alDoB = new ArrayList<String>();
			duplicateDoBs.put(dob, alDoB);
		}
		alDoB.add(Misc.stringArrayToString(f, "\t"));
	}




	private void parseCorrectData() throws IOException {
		avatarIdCorrData = new HashMap<String, String[]>();
		mrnCorrData = new HashMap<String, String[]>();
		BufferedReader in = IO.fetchBufferedReader(corrData);
		String[] fields = null;
		String line = null;
		while ((line = in.readLine())!=null) {
			line = line.trim();
			if (line.startsWith("avatarId1") || line.length()==0) continue;
			//avatarId1,mrn,hciPersonId,firstName,middleName,lastName,birthDate,gender
			fields = Misc.COMMA.split(line);
			String avatarId = fields[0];
			if (avatarIdCorrData.containsKey(avatarId)) throw new IOException("Duplicate AvatarID in the correct data, see "+line);
			avatarIdCorrData.put(avatarId, fields);
			//add via mrn
			if (mrnCorrData.containsKey(fields[1])) throw new IOException("Duplicate MRN in the correct data, see "+line);
			mrnCorrData.put(fields[1], fields);
		}
		in.close();
	}

	public static void main(String[] args) throws IOException {
		new FixingSMM();
		
		

	}

}
