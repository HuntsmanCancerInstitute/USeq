package edu.utah.seq.run.avproj.adw;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;

public class AWSRepoList {
	
	//fields
	public static final Pattern SPACES = Pattern.compile(" +");
	private ArrayList<String> awsObjectNames = new ArrayList<String>();

	//constructor
	public AWSRepoList(File awsRepoListFile) throws IOException {
		//parse the file extracting just the object names
		String line;
		String[] split;
		BufferedReader in = new BufferedReader (new FileReader(awsRepoListFile));
		while ((line = in.readLine())!= null) {
			split = SPACES.split(line);
			//Misc.printExit(line+"  "+split.length);
			//2022-05-31    14:35:19   23107   Patients/AC4Nx5JK/Caris/TN21-177236_2021-12-10/Alignment/TN21-177236_2021-12-10_TumorDNA/Logs.zip
			if (split.length == 4) awsObjectNames.add(split[3]);
		}
		in.close();
	}
	
	public HashSet<String> parseSubDirNamesAfterParent(String parentDirName){
		HashSet<String> subDirs = new HashSet<String>();
		Pattern pat = Pattern.compile("/"+parentDirName+"/([\\w\\.-]+)/");
		Matcher mat;
		for (String name: awsObjectNames) {
			mat = pat.matcher(name);
			if (mat.find()) {
				subDirs.add(mat.group(1));
			}
		}
		return subDirs;
	}
	
	public static void main(String[] args) throws IOException {
		File test = new File("/Users/u0028003/Downloads/PatientMolecularRepoFiles/awsRepoList19Apr2023.txt");
		AWSRepoList arl = new AWSRepoList(test);
		IO.pl("NumObs: "+arl.getAwsObjectNames().size());
		HashSet<String> subDirs = arl.parseSubDirNamesAfterParent("Avatar");
		IO.pl(subDirs);
	}

	public ArrayList<String> getAwsObjectNames() {
		return awsObjectNames;
	}
}
