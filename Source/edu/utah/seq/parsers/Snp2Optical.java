package edu.utah.seq.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Made to parse Illumina snp chip data exported from GenomeStudio into Optical format.  Steve Guthery.
 * @author Nix*/
public class Snp2Optical {

	LinkedHashMap<String, ArrayList<Float>> snpNameSample = new LinkedHashMap<String, ArrayList<Float>>();
	ArrayList<String> sampleNamesAL = new ArrayList<String>();
	String[] sampleNames;
	int numberSamples;
	int twiceNumberSamples;
	LinkedHashSet<String> alleles = new LinkedHashSet<String>();
	
	public Snp2Optical(String[] args) {
		File snpData = new File(args[0]);
		System.out.print("Loading...");
		loadSnpData(snpData);
		sampleNames = Misc.stringArrayListToStringArray(sampleNamesAL);
		numberSamples = sampleNames.length;
		twiceNumberSamples = numberSamples*2;
		System.out.println("\n\t"+snpNameSample.size()+" Snps found");
		System.out.println("\nSample names:");
		Misc.printArray(sampleNamesAL);
		sampleNamesAL = null;
		
		//make hash to see if there is an issue
		HashSet<String> colNames = new HashSet<String>();
		for (String s: sampleNames) colNames.add(s);
		if (colNames.size() != numberSamples) Misc.printErrAndExit("Num of samples differ!");
		System.out.println(numberSamples+" Samples found");
		
		//printdata
		File out = new File(snpData.getParentFile(), Misc.removeExtension(snpData.getName())+"_s2o.txt.gz");
		try {
			Gzipper gz = new Gzipper(out);
			printOpticalFormat(gz);
			gz.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private void printOpticalFormat(Gzipper gz) throws Exception {
		
		//print header
		gz.print("SNP\tCoor\tAlleles");
		for (int i=0; i<numberSamples; i++){
			gz.print("\t");
			gz.print(sampleNames[i]); 
			gz.print("A\t");
			gz.print(sampleNames[i]);
			gz.print("B");
		}
		gz.println();
		
		//for each snp
		for (String snpName : snpNameSample.keySet()){
			ArrayList<Float> samplesAL = snpNameSample.get(snpName);
			int num = samplesAL.size();
			if (num != twiceNumberSamples) Misc.printErrAndExit("\nDifferent number of values than samples!\n");
			
			//print general snp info
			////snp name tab position
			gz.print(snpName); 
			//fudge the alleles
			gz.print("\tXX");
			//for each sample, this is going to throw an out of bounds error if sample number isn't correct
			for (int i=0; i<num; i++){
				gz.print("\t");
				//X
				gz.print(samplesAL.get(i++)); gz.print("\t");
				//Y
				gz.print(samplesAL.get(i));
			}
			gz.println();
		}
		
	}

	private String fetchAlleles(ArrayList<String[]> samplesAL) {
		alleles.clear();
		int num = samplesAL.size();
		for (int i=0; i< num; i++){
			String[] s = samplesAL.get(i);
			alleles.add(s[3]);
			if (alleles.size() == 2) return Misc.hashSetToString(alleles, "");
			alleles.add(s[4]);
			if (alleles.size() == 2) return Misc.hashSetToString(alleles, "");
		}
		return "XX";
	}

	private void loadSnpData(File f) {
		try {
			BufferedReader in = IO.fetchBufferedReader(f);
			String line;
			String[] fields;
			String oldName = "";
			//find header line
			//   0          1          2             3          4       5   6
			//SNPName	Position	SampleID	Allele1-Top	Allele2-Top	X	Y
			while ((line = in.readLine()) != null){
				if (line.contains("Allele1")) break;
			}
			//for each data line
			int counter = 0;
			while ((line = in.readLine()) != null){
				if (counter++ > 1000000){
					System.out.print(".");
					counter = 0;
				}
				fields = Misc.TAB.split(line);
				String snpName = fields[0]+"\t"+fields[1];
				ArrayList<Float> samples = snpNameSample.get(snpName);
				if (samples == null){
					samples = new ArrayList<Float>();
					snpNameSample.put(snpName, samples);
				}
				samples.add(new Float(Float.parseFloat(fields[5])));
				samples.add(new Float(Float.parseFloat(fields[6])));
				//check sample name
				if (fields[2].equals(oldName) == false){
					oldName = fields[2];
					sampleNamesAL.add(oldName);
				}
			}
			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

	public static void main(String[] args){
		new Snp2Optical(args);
	}
	
	
}
