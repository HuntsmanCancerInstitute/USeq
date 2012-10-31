package util.bio.annotation;
import java.io.*;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.*;

import util.bio.parsers.*;


public class ExportTrimmedGenes {

	private File ucscTableFile;	
	private UCSCGeneModelTableReader genes;
	private boolean justUTR = false;
	private boolean justFirstIntron = false;

	public ExportTrimmedGenes (String[] args){
		try {
			//process the args
			processArgs(args);
			
			//load the gene models
			genes = new UCSCGeneModelTableReader(ucscTableFile, 0);
			UCSCGeneLine[] lines = genes.getGeneLines();
			
			String utr = "UTRs+1stCDSIntron";
			if (justUTR) utr = "JustUTRs";
			else if (justFirstIntron) utr = "Just1stCDSIntron";
			File trimmedFile = new File (ucscTableFile.getParentFile(), Misc.removeExtension(ucscTableFile.getName())+"_trimmed_"+ utr +".ucsc");
			Gzipper out = new Gzipper(trimmedFile);
			
			//for each line
			for (int i=0; i< lines.length; i++){
				ExonIntron[] exons = lines[i].getExons();
				//any introns?
				if (exons.length<2) continue;
				
				//make new set of exons, only want those that include or preceed the coding start
				ArrayList<ExonIntron> eiAL = new ArrayList<ExonIntron>();
				
				if (justUTR){
					//plus stranded
					if (lines[i].getStrand().equals("+")){
						int codingStart = lines[i].getCdsStart();
						int end = 0;
						for (ExonIntron ei: exons){
							//before coding sequence?
							if (ei.getEnd() < codingStart){
								eiAL.add(ei);
								end = ei.getEnd();
							}
							//contained?
							else if (ei.contains(codingStart)) {
								eiAL.add(ei);
								end = ei.getEnd();
								break;
							}
							//must be past
							else break;
						}
						//reset end of gene and coding sequence
						if (end !=0){
							lines[i].setCdsEnd(end);
							lines[i].setTxEnd(end);
						}					
					}
					
					//minus strand
					else {
						int codingStart = lines[i].getCdsEnd()-1;
						int end = 0;
						//walk from rear
						for (int j=exons.length-1; j>=0; j--){
							//before?
							if (codingStart < exons[j].getStart()){
								eiAL.add(exons[j]);
								end = exons[j].getStart();
							}
							//contained?
							else if (exons[j].contains(codingStart)) {
								eiAL.add(exons[j]);
								end = exons[j].getStart();
								break;
							}
							//beyond
							else break;
						}
						//reset end of gene and coding sequence
						if (end !=0){
							lines[i].setCdsStart(end);
							lines[i].setTxStart(end);
						}
					}
				}
				
				else if (justFirstIntron){
					//plus stranded
					if (lines[i].getStrand().equals("+")){
						int codingStart = lines[i].getCdsStart();
						int end = 0;
						int numExons = exons.length;
						for (int j=0; j< numExons; j++){
							if (exons[j].contains(codingStart)){
								//does it contain the stop too?
								if (exons[j].contains(lines[i].getCdsEnd())) {
									eiAL=null;
									break;
								}
								eiAL.add(exons[j]);
								lines[i].setCdsStart(exons[j].getStart());
								lines[i].setTxStart(exons[j].getStart());
								
								//add next?
								j++;
								if (j<numExons){
									eiAL.add(exons[j]);
									end = exons[j].getEnd();
									lines[i].setCdsEnd(end);
									lines[i].setTxEnd(end);
									break;
								}
								//no next so no intron, so skip
								else eiAL = null;
								
							}
						}					
					}
					
					//minus strand
					else {
						int codingStart = lines[i].getCdsEnd()-1;
						int end = 0;
						//walk from rear
						for (int j=exons.length-1; j>=0; j--){
							if (exons[j].contains(codingStart)){
								//does it also contain the cds end?
								if (exons[j].contains(lines[i].getCdsStart())){
									eiAL = null;
									break;
								}
								eiAL.add(exons[j]);
								lines[i].setCdsEnd(exons[j].getEnd());
								lines[i].setTxEnd(exons[j].getEnd());
								
								//add next
								j--;
								if (j>=0){
									eiAL.add(exons[j]);
									end = exons[j].getStart();
									lines[i].setCdsStart(end);
									lines[i].setTxStart(end);
									break;
								}
								else eiAL = null;
							}
						}
					}
				}
				
				//utrs and first cs exon
				else {
					//plus stranded
					if (lines[i].getStrand().equals("+")){
						int codingStart = lines[i].getCdsStart();
						int end = 0;
						for (ExonIntron ei: exons){
							//add em, might be the last
							eiAL.add(ei);
							end = ei.getEnd();
							if (ei.contains(codingStart) == false && (ei.getEnd() <= codingStart)== false) break;
						}
						//reset end of gene and coding sequence
						if (end !=0){
							lines[i].setCdsEnd(end);
							lines[i].setTxEnd(end);
						}					
					}
					
					//minus strand
					else {
						int codingStart = lines[i].getCdsEnd()-1;
						int end = 0;
						//walk from rear
						for (int j=exons.length-1; j>=0; j--){
							eiAL.add(exons[j]);
							end = exons[j].getStart();
							if (exons[j].contains(codingStart) == false && (codingStart <= exons[j].getStart()) == false) break;
						}
						//reset end of gene and coding sequence
						if (end !=0){
							lines[i].setCdsStart(end);
							lines[i].setTxStart(end);
						}
					}
				}
				
				
				
				//print it?
				if (eiAL != null && eiAL.size() > 1){
					ExonIntron[] ei = new ExonIntron[eiAL.size()];
					eiAL.toArray(ei);
					lines[i].setExons(ei);
					out.println(lines[i].toString());
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
		new ExportTrimmedGenes(args);
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
					case 'u': justUTR = true; break;
					case 'i': justFirstIntron = true; break;
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
		if (justUTR == true && justFirstIntron == true) Misc.printExit("\nError: must select either -i or -u but not both.  Select neither for default.\n");
	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              Export Trimmed Genes    May 2012                    **\n" +
				"**************************************************************************************\n" +
				"EE takes a UCSC Gene table and clips each gene back to the first intron closed by a\n" +
				"coding sequence exon. Thus these include all of the 5'UTRs. Genes with no introns are\n" +
				"removed.\n\n"+

				"Parameters:\n"+
				"-g Full path file text for the UCSC Gene table.\n"+
				"-u Print just UTRs, defaults to UTRs plus 1st CDS intron with flanking exon.\n"+
				"-i Print just 1st CDS intron with flanking exons.\n"+

				"\nExample: java -Xmx1000M -jar pathTo/T2/Apps/ExportTrimmedGenes -u -g \n" +
				"      /user/Jib/ucscPombe.txt\n"+

		"**************************************************************************************\n");		
	}		



}

