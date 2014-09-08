package util.bio.annotation;
import java.io.*;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.*;

import util.bio.parsers.*;

public class ExportExons {

	private File ucscTableFile;	
	private UCSCGeneModelTableReader genes;
	private int adder = 0;
	private boolean trimUTRs = false;
	private boolean addExonNumbers = false;
	
	
	
	private interface NameString {
		String createNameString(String name, String strand, int exonNumber);
	}
	
	
	public ExportExons (String[] args){
		//Parse args
		processArgs(args);
		
		//create NameString command basd on addExonNumbers
		NameString nameStringCommand;
		if (this.addExonNumbers) {
			nameStringCommand = new NameString() {
				public String createNameString(String name, String strand, int exonNumber) {
					String nameScoreStrand = "\t" + name + "_" + String.valueOf(exonNumber) +  "\t0\t" + strand;
					return nameScoreStrand;		
				};
			};
		} else {
			nameStringCommand = new NameString() {
				public String createNameString(String name, String strand, int exonNumber) {
					String nameScoreStrand = "\t" + name  +  "\t0\t" + strand;
					return nameScoreStrand;
				};
			};
		}
		
		
		try {
			
			genes = new UCSCGeneModelTableReader(ucscTableFile, 0);
			UCSCGeneLine[] lines = genes.getGeneLines();
			File exonsFile = new File (ucscTableFile.getParentFile(), Misc.removeExtension(ucscTableFile.getName())+"_Exons.bed");
			PrintWriter out = new PrintWriter (new FileWriter(exonsFile));
			if (trimUTRs){
				for (int i=0; i< lines.length; i++){
					ExonIntron[] exons = lines[i].getExons();
					String chr = lines[i].getChrom();
					int startCoding = lines[i].getCdsStart();
					int endCoding = lines[i].getCdsEnd();
					String strand = lines[i].getStrand();
					String name = lines[i].getDisplayNameThenName();
					
					for (int j=0; j< exons.length; j++){
						int start = exons[j].getStart();
						int end = exons[j].getEnd();
						//trim utrs?
						if (start< startCoding) start = startCoding;
						if (end > endCoding) end = endCoding;
						if ((end-start) <=0) continue;
						//adders
						start = start - adder;
						if (start < 0) start =0;
						end = end + adder;
						
						//Determine exon number
						int exon = -1;
						if (lines[i].getStrand().equals("-")) {
							exon = exons.length - j;
						} else if (lines[i].getStrand().equals("+")) {
							exon = j + 1;
						}
						
						//chr start stop name score strand
						String nameScoreStrand = nameStringCommand.createNameString(name, strand, exon);
						out.println(chr+"\t"+start+"\t"+end + nameScoreStrand);
					}
				}
			}
			else {
				for (int i=0; i< lines.length; i++){
					ExonIntron[] exons = lines[i].getExons();
					String chr = lines[i].getChrom();
					String strand = lines[i].getStrand();
					String name = lines[i].getDisplayNameThenName();
					
					
					for (int j=0; j< exons.length; j++){
						int start = exons[j].getStart()-adder;
						if (start < 0) start =0;
						
						//Determine exon number
						int exon = -1;
						if (lines[i].getStrand().equals("-")) {
							exon = exons.length - j;
						} else if (lines[i].getStrand().equals("+")) {
							exon = j + 1;
						}
						
						String nameScoreStrand = nameStringCommand.createNameString(name, strand, exon);
						out.println(chr+"\t"+start+"\t"+(exons[j].getEnd()+adder) + nameScoreStrand);
					}
				}
			}
			out.close();
		}
		catch (Exception e){
			e.printStackTrace();
		}

	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		
		
		
		new ExportExons(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'g': ucscTableFile = new File (args[i+1]); i++; break;
					case 'a': adder = Integer.parseInt(args[++i]); break;
					case 'u': trimUTRs=true; break;
					case 'n': addExonNumbers=true; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		if (ucscTableFile == null || ucscTableFile.canRead() == false) Misc.printExit("\nError: cannot find your UCSC table file!\n");
	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              Export Exons   Sept 2013                            **\n" +
				"**************************************************************************************\n" +
				"EE takes a UCSC Gene table and prints the exons to a bed file.\n\n"+

				"Parameters:\n"+
				"-g Full path file text for the UCSC Gene table.\n"+
				"-a Expand the size of each exon by X bp, defaults to 0\n"+
				"-u Remove UTRs if present, defaults to including\n"+
				"-n Append exon numbers to the gene name field.  This makes the bed file compatible \n" +
				"      with DRDS\n" +

				"\nExample: java -Xmx1000M -jar pathTo/T2/Apps/ExportExons -g /user/Jib/ucscPombe.txt\n" +
				"      -a 50\n"+

		"**************************************************************************************\n");		
	}		



}

