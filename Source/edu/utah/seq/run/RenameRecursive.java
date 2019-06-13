package edu.utah.seq.run;

import java.io.File;
import java.util.ArrayList;

import util.gen.IO;

public class RenameRecursive {

	public static void main(String[] args) {
		if (args == null || args.length != 1) return;
		
		File toParse = new File(args[0]);
		File[] patientDirs = IO.extractOnlyDirectories(toParse);
		
		//for each patient dir
		for (File patientDir : patientDirs){
			IO.pl("Proc "+patientDir);
			File[] dirs = IO.extractOnlyDirectories(patientDir);
			
			//for each dir
			for (File subDir: dirs){
				//IO.pl("\t"+subDir.getName());
				//skip Fastq
				//if (subDir.getName().equals("Fastq")) continue;
				ArrayList<File> filesToCheck = IO.fetchFilesDirsRecursively(subDir);
				
				//for each File
				for (File file: filesToCheck){
					//look for Exome
					if (file.getName().contains("Exome")){
						String modName = file.getName().replace("Exome", "DNA");
						File newFile = new File(file.getParentFile(), modName);
						//IO.pl("\t\t"+file);
						//IO.pl("\t\t"+newFile);
						file.renameTo(newFile);
					}
					else if (file.getName().contains("Transcriptome")){
						String modName = file.getName().replace("Transcriptome", "RNA");
						File newFile = new File(file.getParentFile(), modName);
						//IO.pl("\t\t"+file);
						//IO.pl("\t\t"+newFile);
						file.renameTo(newFile);
					}
				}
				
			}
		}

	}

}
