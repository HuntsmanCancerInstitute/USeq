package edu.utah.seq.methylation454;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Num;

public class Reference {

	//fields
	private String name;
	private String sequence;
	private String chromosome = "chr?";
	private int startPosition = 0;
	private int[] cpGPositions;
	private int[] cpGGenomicPositions;
	
	//constructor
	/**Alignment file line for a reference.*/
	public Reference (String line){
		sequence = line.split("\\s+")[1];
		cpGPositions = findCpGs(sequence);
		cpGGenomicPositions = findCpGs(sequence.replaceAll("-", ""));
	}
	
	
	/**Case sensitive.*/
	public static int[] findCpGs(String seq){
		Pattern p = Pattern.compile("C-*G");
		Matcher m = p.matcher(seq);
		ArrayList al = new ArrayList();
		while (m.find()) {
			int base = m.start();
			al.add(new Integer(base));
			//System.out.println(base+" "+seq.charAt(base));
		}   
		return Num.arrayListOfIntegerToInts(al);
	}


	public String getChromosome() {
		return chromosome;
	}


	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}


	public int getStartPosition() {
		return startPosition;
	}


	public void setStartPosition(int startPosition) {
		this.startPosition = startPosition;
	}


	public int[] getCpGGenomicPositions() {
		return cpGGenomicPositions;
	}


	public int[] getCpGPositions() {
		return cpGPositions;
	}


	public String getName() {
		return name;
	}


	public String getSequence() {
		return sequence;
	}


	public void setName(String name) {
		this.name = name;
	}

}
