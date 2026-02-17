package edu.utah.seq.pmr;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;

/**Identifies partial Avatar datasets that should be deleted given more recent complete analysis*/
public class AvatarPMRPruner {
	
	//user defined fields
	private File pmrFileList = null;
	private HashSet<String> savedPaths = new HashSet<String>();
	private HashMap<String, ArrayList<PMRDataset>> pmrIdAvatarDatasets = new HashMap<String, ArrayList<PMRDataset>>();

	public AvatarPMRPruner (String[] args) {
		long startTime = System.currentTimeMillis();
		try {
			processArgs(args);

			parsePMRListFile();
			IO.pl(pmrIdAvatarDatasets.size()+"\tDatasets");
			
			walkPatientDatasets();
			
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
			IO.pl("\nDone! "+Math.round(diffTime)+" Sec\n");
		} catch (Exception e) {
			IO.el("\nERROR running the SearchPMR...");
			e.printStackTrace();
			System.exit(1);
		}
	}
	

	private void walkPatientDatasets() {
		for (String pmrId: pmrIdAvatarDatasets.keySet()) {
			ArrayList<PMRDataset> ad = pmrIdAvatarDatasets.get(pmrId);
			if (ad.size() == 1) continue;
			//Any NAs?
			boolean naFound = false;
			for (PMRDataset d: ad) {
				if (d.containsNA()) {
					naFound = true;
					break;
				}
			}
			if (naFound == false) continue;
			
			ArrayList<String> dups = new ArrayList<String>();
			HashSet<String> tes = new HashSet<String>();
			HashSet<String> tts = new HashSet<String>();
			boolean dupsFound = false;
			
			//for each of the datasets
			for (PMRDataset d: ad) {
				dups.add(d.getPath());
				if (d.te.equals("NA")== false) {
					if (tes.contains(d.te)) {
						dupsFound = true;
					}
					else tes.add(d.te);
				}
				if (dupsFound == false) {
					if (d.tt.equals("NA")== false) {
						if (tts.contains(d.tt)) dupsFound = true;
						else tts.add(d.tt);
					}
				}
			}
			if (dupsFound) IO.pl(Misc.stringArrayListToString(dups, "\n")+"\n");
		}
	}


	private void parsePMRListFile() throws IOException {
		BufferedReader in = IO.fetchBufferedReader(pmrFileList);
		String line;
		String[] f;
		String[] u;
		
		while ((line = in.readLine())!=null) {
			line = line.trim();
			if (line.contains("Avatar")) {
				//2023-03-01    09:04:43     2938    Patients/AA2mF6Vy/Avatar/A032049_SL419345_SL419548_SL420681/ClinicalReport/A032049_SL419345_SL419548_SL420681_IDTv1_SAR_F.json
				String patient = line.substring(line.indexOf("Patients/")+9);
				//AA2mF6Vy/Avatar/A032049_SL419345_SL419548_SL420681/ClinicalReport/
				//   0        1               2
				f = Misc.FORWARD_SLASH.split(patient);
				u = Misc.UNDERSCORE.split(f[2]);
				if (u.length != 4) throw new IOException("Failed to parse 4 fields from "+f[2]+" in "+line);
				
				PMRDataset ad = new PMRDataset(f[0], u[0], u[1], u[2], u[3]);
				String path = ad.getPath();
				//already saved?
				if (savedPaths.contains(path) == false) {
					savedPaths.add(path);
					//IO.pl(path);
					//fetch or make ArrayList
					ArrayList<PMRDataset> al = pmrIdAvatarDatasets.get(ad.pmrId);
					if (al == null) {
						al = new ArrayList<PMRDataset>();
						pmrIdAvatarDatasets.put(ad.pmrId, al);
					}
					al.add(ad);
				}
			}
		}
		in.close();
	}
	
	private class PMRDataset {
		
		String pmrId;
		String orienId;
		String ne;
		String te;
		String tt;
		
		public PMRDataset(String pmrId, String orienId, String ne, String te, String tt) {
			this.pmrId = pmrId;
			this.orienId = orienId;
			this.ne = ne;
			this.te = te;
			this.tt = tt;
		}
		
		public boolean containsNA() {
			if (ne.equals("NA") || te.equals("NA") || tt.equals("NA")) return true;
			return false;
		}

		public String getPath() {
			StringBuilder sb = new StringBuilder();
			sb.append(pmrId); sb.append("/");
			sb.append("Avatar/");
			sb.append(orienId); sb.append("_");
			sb.append(ne); sb.append("_");
			sb.append(te); sb.append("_");
			sb.append(tt);
			return sb.toString();
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AvatarPMRPruner(args);
	}		

	/**This method will process each argument and assign new variables
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException {

		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'p': pmrFileList = new File(args[++i]).getCanonicalFile(); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					e.printStackTrace();
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		if (pmrFileList == null || pmrFileList.canRead()== false) {
			Misc.printErrAndExit("Error: failed to find your pmr list file? -p "+pmrFileList);
		}
	}


	public static void printDocs(){
		IO.pl("\n" +
				"xxx\n"+

				"**************************************************************************************\n");
	}
}
