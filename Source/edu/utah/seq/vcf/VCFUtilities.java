package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

public class VCFUtilities {

	public VCFUtilities() {
		// TODO Auto-generated constructor stub
	}
	
	public static File unzipTabix(File tabixFile, String pathToTabix) {
		//Make sure the file is the correct format and grab the file prefix
		Pattern fp = Pattern.compile("(.+?.vcf).gz");
		Matcher fm = fp.matcher(tabixFile.getName());
		
		//Exit if not the correct format
		if (!fm.matches()) {
			System.out.println("Tabix-compressed VCF file does not have the correct extension, exiting");
			System.exit(1);
		}
		
		File uncompressed = new File(fm.group(1));
		uncompressed.deleteOnExit();
		
		ProcessBuilder pb = new ProcessBuilder(pathToTabix + "/bgzip","-c","-d",tabixFile.getName());
		
		
		try {
			Process p = pb.start();
			BufferedReader outStream = new BufferedReader(new InputStreamReader(p.getInputStream()));
			BufferedReader errStream = new BufferedReader(new InputStreamReader(p.getErrorStream()));
			StringBuffer sb = new StringBuffer("");
			BufferedWriter bw = new BufferedWriter(new FileWriter(uncompressed));
			
			String line = null;
			
			while ((line = outStream.readLine()) != null) {
				bw.write(line + "\n");
			}
			
			while((line = errStream.readLine()) != null) {
				sb.append(line + "\n");
			}
			
			int retVal = p.waitFor();
			if (retVal != 0) {
				System.out.println("File decompression failed: " +sb.toString());
				System.exit(1);
			}
			
			bw.close();
			
		} catch (IOException ioex ) {
			System.out.println("Failed to read/write file stream : " + ioex.getMessage());
			ioex.printStackTrace();
			System.exit(1);
		} catch (InterruptedException irex) {
			System.out.println("Gunzip was interruped before it was finished: " + irex.getMessage());
			irex.printStackTrace();
			System.exit(1);
		} 
		return uncompressed;
	}
	
	public static File createTabix(File vcfFile, String pathToTabix) {
		//Make sure the file is the correct format and grab the file prefix
		Pattern fp = Pattern.compile(".+?.vcf");
		Matcher fm = fp.matcher(vcfFile.getName());
		
		//Exit if not the correct format
		if (!fm.matches()) {
			System.out.println("VCF file does not have the correct extension, exiting");
			System.exit(1);
		}
		
		//Create file pointing to compressed file
		File compressed = new File(fm.group(0) + ".gz");
		
		ProcessBuilder pb1 = new ProcessBuilder(pathToTabix + "/bgzip",vcfFile.getName());
		ProcessBuilder pb2 = new ProcessBuilder(pathToTabix + "/tabix","-p","vcf",compressed.getName());
		
		ArrayList<ProcessBuilder> pbList = new ArrayList<ProcessBuilder>(); 
		pbList.add(pb1);
		pbList.add(pb2);
		
		try {
			for (ProcessBuilder pb: pbList) {
				pb.redirectErrorStream(true);
				Process p = pb.start();
				
				BufferedReader outStream = new BufferedReader(new InputStreamReader(p.getInputStream()));
				StringBuffer sb = new StringBuffer("");
				
				String line = null;
				
				while ((line = outStream.readLine()) != null) {
					sb.append(line + "\n");
				}
				
				int retVal = p.waitFor();
				if (retVal != 0) {
					System.out.println("File compression failed failed: " + pb.command().toString() + ": " + sb.toString());
					System.exit(1);
				}
			}
		} catch (IOException ioex ) {
			System.out.println("Failed to read/write file stream : " + ioex.getMessage());
			ioex.printStackTrace();
			System.exit(1);
		} catch (InterruptedException irex) {
			System.out.println("Gunzip was interruped before it was finished: " + irex.getMessage());
			irex.printStackTrace();
			System.exit(1);
		} 
		
		//mark VCF file for deletion
		vcfFile.deleteOnExit();
		
		return compressed;
	}

}
