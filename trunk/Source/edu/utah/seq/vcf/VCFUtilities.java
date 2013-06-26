package edu.utah.seq.vcf;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.regex.Pattern;
import java.util.regex.Matcher;


public class VCFUtilities {

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
		
		ProcessBuilder pb = new ProcessBuilder(pathToTabix + "/bgzip","-c","-d",tabixFile.getAbsoluteFile().toString());
		
		
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
	
	public static int readsToChunk = 500000;
	
	public static int countReads(File vcfFile) {
		BufferedReader in=null;
		String line = null;
		int recordCount=0;
		boolean found = false;
		try {
			in = new BufferedReader(new FileReader(vcfFile));
			while ((line=in.readLine()) != null) {
				if (line.startsWith("#CHROM")) {
					found = true;
					break;
				}
			}
			
			if (!found) {
				System.out.println("Could not find the #Chrom line, exiting");
				System.exit(1);
			}
			
			while ((line=in.readLine()) != null){
				recordCount++;
			}
			
			if (recordCount == 0) {
				System.out.println("No VCF records found, exiting");
				System.exit(1);
			}
			
			in.close();
		} catch (FileNotFoundException fnfe) {
			System.out.println("Specified file cannot be found, exiting: " + vcfFile.getAbsolutePath());
			System.exit(1);
		} catch (IOException ioex) {
			System.out.println("Error reading vcf file, exiting: " + vcfFile.getAbsolutePath());
			System.exit(1);
		} catch (Exception e) {
			System.err.println("\nAborting, problem parsing vcf file -> "+vcfFile);
			e.printStackTrace();
			System.exit(1);
		} finally{
			try {
				in.close();
			} catch (IOException e) {}
		}
		return recordCount;
	}
	
	public static void catFiles(ArrayList<File> fileList, File dest) {
		try {
			ArrayList<String> command  = new ArrayList<String>();
			command.add("cat");

			for (File f: fileList) {
				if (!f.exists()) {
					System.out.println("[catFiles] Expected file does not exist: " + f.getAbsolutePath());
					System.exit(1);
				}
				command.add(f.getAbsolutePath());
			}
			
			
			ProcessBuilder pb = new ProcessBuilder(command);
			Process p = pb.start();
			
			BufferedInputStream bis = new BufferedInputStream(p.getInputStream());
			BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(dest));
			
			
			byte[] buffer = new byte[1024*1024*10];
			int n = -1;
			
			while((n = bis.read(buffer))!=-1) {
			  bos.write(buffer,0,n);
			}
		

			int val = p.waitFor();
			bos.close();
			bis.close();

			
			if (val != 0) {
				System.out.println("[catFiles] Error while merging the text files files: " + dest.getAbsolutePath());
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
			
			
			
		} catch (IOException ioex) {
			System.out.println("[mergeVcf] IO Exception while trying to merge the VCF files: " + dest.getAbsolutePath());
			System.exit(1);
		} catch (InterruptedException ieex) {
			System.out.println("[mergeVcf] Process was interrupted while trying to merge the VCF files: " + dest.getAbsolutePath());
			System.exit(1);
		}
	}
	
	public static void mergeVcf(ArrayList<File> vcfList, File dest) {
		try {
			ArrayList<String> command  = new ArrayList<String>();
			command.add("/tomato/app/vcftools/vcf-concat");

			for (File bam: vcfList) {
				if (!bam.exists()) {
					System.out.println("[mergeVcf] Expected file does not exist: " + bam.getAbsolutePath());
					System.exit(1);
				}
				command.add(bam.getAbsolutePath());
			}
			
			
			ProcessBuilder pb = new ProcessBuilder(command);
			Process p = pb.start();
			
			BufferedInputStream bis = new BufferedInputStream(p.getInputStream());
			BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(dest));
			
			
			byte[] buffer = new byte[1024*1024*10];
			int n = -1;
			
			while((n = bis.read(buffer))!=-1) {
			  bos.write(buffer,0,n);
			}
		

			int val = p.waitFor();
			bos.close();
			bis.close();

			
			if (val != 0) {
				System.out.println("[mergeVcf] Error while merging the VCF files: " + dest.getAbsolutePath());
				BufferedReader br2 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
				String line2 = null;
				while((line2 = br2.readLine()) != null) {
					System.out.println(line2);
				}
				System.exit(1);
			}
			
			
			
		} catch (IOException ioex) {
			System.out.println("[mergeVcf] IO Exception while trying to merge the VCF files: " + dest.getAbsolutePath());
			System.exit(1);
		} catch (InterruptedException ieex) {
			System.out.println("[mergeVcf] Process was interrupted while trying to merge the VCF files: " + dest.getAbsolutePath());
			System.exit(1);
		}
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
		
		ProcessBuilder pb1 = new ProcessBuilder(pathToTabix + "/bgzip",vcfFile.getAbsoluteFile().toString());
		ProcessBuilder pb2 = new ProcessBuilder(pathToTabix + "/tabix","-p","vcf",compressed.getAbsoluteFile().toString());
		
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
