package edu.utah.seq.methylation454;
import java.io.*;
import java.util.*;

public class MethylationParser {
	//fields
	private Sample[] samples;
	//internal fields
	private DataLine[] dataLines;
	private DataLine dataLine;
	private ArrayList sampleAL = new ArrayList();
	private Sample sample = null;
	private ArrayList ampliconAL = new ArrayList();
	private Amplicon amplicon = null;
	private ArrayList cpGAL = new ArrayList();
	private CpG cpG = null;
	private ArrayList readAL = new ArrayList();
	
	//constructor
	public MethylationParser(File dataFile){
		parseIt(dataFile);
	}
	
	//methods
	public void parseIt (File dataFile){
		//parse to DataLine[]
		parseDataLines(dataFile);
		//sort
		Arrays.sort(dataLines);
		//build Sample[]
		makeSamples();
	}
	
	public void makeSamples(){
		//set first
		sample = new Sample(dataLines[0].getSample());
		amplicon = new Amplicon(dataLines[0].getAmplicon());
		cpG = new CpG(dataLines[0].getCpG());
		readAL.add(new BaseRead(dataLines[0]));
		//walk through sorted data lines
		for (int i=1; i< dataLines.length; i++){
			dataLine = dataLines[i];
			//check sample
			if (sample.getId().equals(dataLine.getSample()) == false){
				closeSetReads();
				closeSetCpGs();
				closeSetAmplicons();
				//System.out.println("Sample\t"+dataLine);
			}
			//check amplicon
			else if (amplicon.getId().equals(dataLine.getAmplicon()) == false){
				closeSetReads();
				closeSetCpGs();
				//System.out.println("Amplic\t"+dataLine);
			}
			//check CpG
			else if (cpG.getId().equals(dataLine.getCpG()) == false){
				closeSetReads();
				//System.out.println("CpG   \t"+dataLine);
			}
			//by default must be a new read
			else {
				readAL.add(new BaseRead(dataLine));
				//System.out.println("Read  \t"+dataLine);
			}
		}
		//set last 
		closeSetReads();
		closeSetCpGs();
		closeSetAmplicons();
		//make samples
		samples = new Sample[sampleAL.size()];
		sampleAL.toArray(samples);
	}
	
	public void closeSetAmplicons(){
		//close old amplicons
		Amplicon[] amps = new Amplicon[ampliconAL.size()];
		ampliconAL.toArray(amps);
		ampliconAL.clear();
		//set amplicon
		sample.setAmplicons(amps);
		sampleAL.add(sample);
		//make new
		sample = new Sample(dataLine.getSample());
	}
	
	public void closeSetCpGs(){
		//close old CpGs
		CpG[] cpGs = new CpG[cpGAL.size()];
		cpGAL.toArray(cpGs);
		cpGAL.clear();
		//set cpGs
		amplicon.setCpGs(cpGs);
		ampliconAL.add(amplicon);
		//make new
		amplicon = new Amplicon (dataLine.getAmplicon());
	}
	
	public void closeSetReads(){
		//close old reads
		BaseRead[] reads = new BaseRead[readAL.size()];
		readAL.toArray(reads);
		readAL.clear();
		//set reads
		cpG.setReads(reads);
		cpGAL.add(cpG);
		//make new
		cpG = new CpG(dataLine.getCpG());
		readAL.add(new BaseRead(dataLine));
	}
	
	public void parseDataLines(File dataFile){
		ArrayList al = new ArrayList();
		try {
			BufferedReader in = new BufferedReader ( new FileReader (dataFile));
			String line;
			while ((line = in.readLine()) !=null){
				line = line.trim();
				if (line.startsWith("#") || line.length()==0) continue;
				DataLine dl = new DataLine(line);
				//is the read a T or C?
				String x = dl.getSequence();
				if (x.equals("T") || x.equals("C")) al.add(new DataLine(line));
			}
			in.close();
		} catch (Exception e) {
			System.out.println("\nProblem parsing data lines for "+dataFile);
			e.printStackTrace();
		}
		dataLines = new DataLine[al.size()];
		al.toArray(dataLines);
	}

	public Sample[] getSamples() {
		return samples;
	}
}
