package edu.utah.seq.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Num;

public class VCFParser {
	//fields
	private File vcfFile;
	public static final Pattern TAB = Pattern.compile("\\t");
	private ArrayList<String> comments = null;
	private long numberVCFRecords = 0;
	private long numberVCFRecordsFailingFilter = 0;
	private HashMap<String, VCFLookUp> chromVCFRecords = null;
	private boolean excludeFilteredRecords = true;
	
	//fields for ripping vcf records
	int chromosomeIndex= 0;
	int positionIndex = 1;
	int referenceIndex = 3;
	int alternateIndex = 4;
	int qualityIndex = 5;
	int filterIndex = 6;
	int infoIndex = 7;
	int sampleIndex = 9;
	int sampleGenotypeIndex =0;
	int sampleScoreIndex =7;
	int sampleRawReadDepthIndex=2;
	int numberColumnsInVCFRecord = 10;

	public VCFParser(File vcfFile, boolean excludeFilteredRecords){
		this. vcfFile = vcfFile;
		this.excludeFilteredRecords = excludeFilteredRecords;
		parseVCF();
	}

	public void parseVCF() {
		numberVCFRecords = 0;
		numberVCFRecordsFailingFilter = 0;
		chromVCFRecords = new HashMap<String, VCFLookUp>();
		comments = new ArrayList<String>();
		
		BufferedReader in = null;
		String[] fields = null;
		String line = null;
		String oldChrom = null;
		int oldPosition = -1;
		int badCounter = 0;
		
		try {
			in  = IO.fetchBufferedReader(vcfFile);

			ArrayList<VCFRecord> records = new ArrayList<VCFRecord>();
			ArrayList<Integer> positions = new ArrayList<Integer>();
			while ((line=in.readLine()) != null){
				//comments?
				if (line.startsWith("#")) {
					comments.add(line);
					continue;
				}
				//correct number of columns?
				fields = TAB.split(line);
				if (fields.length != numberColumnsInVCFRecord){
					System.err.println("Malformed VCF Record skipping -> "+line);
					if (badCounter++ > 100) throw new Exception("\nToo many malformed VCF Records.\n");
					continue;
				}
				//no info? skip it
				if (fields[sampleIndex].startsWith("./.:.:")) continue;
				
				//get chromosome and position
				String chr = fields[chromosomeIndex];
				VCFRecord vcf = new VCFRecord(fields, this);
				//does it pass filter?
				if (excludeFilteredRecords == true && vcf.getFilter().equals("PASS") == false){
					numberVCFRecordsFailingFilter++;
					continue;
				}
				int position = vcf.getPosition();
				//initialize 
				if (oldChrom == null) {
					oldChrom = chr;					
				}
				else if (chr.equals(oldChrom) == false){
					//new chrom so close old
					setChromData(oldChrom, records, positions);
					oldPosition = -1;
					oldChrom = chr;					
				}
				//check position 
				if (oldPosition >= position) throw new Exception ("\nNew position "+position+" is < or = prior position "+oldPosition+".  Is your vcf file sorted?\n");
				//save
				positions.add(position);
				records.add(vcf);
				numberVCFRecords++;
			}
			//add last
			setChromData(oldChrom, records, positions);
			
		}catch (Exception e) {
			System.err.println("\nAborting, problem parsing vcf record -> "+line);
			e.printStackTrace();
			System.exit(1);
		} finally{
			try {
				in.close();
			} catch (IOException e) {}
		}



	}

	private void setChromData(String chr, ArrayList<VCFRecord> records, ArrayList<Integer> positions) {
		VCFRecord[] v = new VCFRecord[records.size()];
		records.toArray(v);
		records.clear();
		int[] p = Num.arrayListOfIntegerToInts(positions);
		positions.clear();
		chromVCFRecords.put(chr, new VCFLookUp (p, v));

	}

	public File getVcfFile() {
		return vcfFile;
	}

	public ArrayList<String> getComments() {
		return comments;
	}

	public long getNumberVCFRecords() {
		return numberVCFRecords;
	}

	public HashMap<String, VCFLookUp> getChromVCFRecords() {
		return chromVCFRecords;
	}

}
