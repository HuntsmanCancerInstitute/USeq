package util.gen;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
 
public class BucketObjLoadTest {
	
	public static void main(String[] args) throws IOException {
		BufferedReader in = IO.fetchBufferedReader("/Users/u0028003/Downloads/BucketObjTest/hcibioinfo-patient-molecular-repo.list.16Dec2023.txt.gz");
		String line;
		String[] tokens;
		HashMap<String,ArrayList<String>> dirFiles = new HashMap<String,ArrayList<String>>();
		String delimiter = "/";
		long totalSize = 0;
		long numberObjects = 0;
		long startTime = System.currentTimeMillis();
		
		while ((line = in.readLine())!=null) {
			line = line.trim();
			if (line.length()!=0) {
				//this will cut any directory or file name with a space, which is permitted!
				//its always a file, never a directory
				/*
2023-03-01 	09:04:43 	2938 		Patients/AA2mF6Vy/Avatar/A032049_SL419345_SL419548_SL420681/ClinicalReport/A032049_SL419345_SL419548_SL420681_IDTv1_SAR_F.json
2023-01-30 	21:32:56 	7208120597 	Patients/AA2mF6Vy/Avatar/A032049_SL419345_SL419548_SL420681/Fastq/NormalDNA/SL419345_1.fastq.gz
2023-01-30 	21:32:56 	7423738484 	Patients/AA2mF6Vy/Avatar/A032049_SL419345_SL419548_SL420681/Fastq/NormalDNA/SL419345_2.fastq.gz
	0        1             2              3
				*/
				tokens = Misc.WHITESPACE.split(line);
				long size = Long.parseLong(tokens[2]);
				totalSize+= size;
				numberObjects++;
				
				//IO.pl(tokens.length+"\t"+ line);
				
				//create fullPath
				String fullPath = null;
				//no white space
				if (tokens.length == 4) fullPath = tokens[3];
				//with white space in name or directory
				else {
					StringBuilder sb = new StringBuilder(tokens[3]);
					for (int i=4; i< tokens.length; i++) {
						sb.append(" ");
						sb.append(tokens[i]);
					}
					fullPath = sb.toString();
				}
				
				//split full path into dir and name
				String dir = null;
				String fileName = null;
				//is this just in the root dir?
				if (fullPath.contains(delimiter) == false) {
					dir = "ROOT";
					fileName = fullPath;
				}
				else {
					int lastIndex = fullPath.lastIndexOf(delimiter);
					dir = fullPath.substring(0, lastIndex);
					fileName = fullPath.substring(lastIndex+1);
				}
				
				//add to hash
				ArrayList<String> files = dirFiles.get(dir);
				if (files == null) {
					files = new ArrayList<String>();
					dirFiles.put(dir, files);
					if (dir.contains("TL-23-WGGHKD5Z/GermlineVariantCalling")) IO.pl(dir);
				}
				files.add(fileName);
			}
		}
		IO.pl(dirFiles.size());
		IO.pl(numberObjects);
		IO.pl(totalSize);
		
		//Issue is with directories that only contain directories, walk the /x/y/z ?
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("Done! "+Math.round(diffTime)+" sec\n");
		
		IO.pl(dirFiles.get("Patients/ADk3Na5eXy/Tempus/TL-23-WGGHKD5Z/GermlineVariantCalling"));
		
		in.close();
	}
	


	
	
	
	
}
	
	