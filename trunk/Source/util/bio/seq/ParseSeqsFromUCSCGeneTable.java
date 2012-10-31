package util.bio.seq;

import java.io.File;
import java.util.*;
import util.bio.parsers.*;
import java.util.regex.*;
import util.gen.*;

public class ParseSeqsFromUCSCGeneTable {

	public static void main(String[] args) {
		//files
		File geneFile = new File ("/Users/nix/HCI/PIs/Thummel/UCSCGeneTables/flybasegenesApr04.xls");//args[0]);
		File fastaDir = new File ("/Users/nix/HCI/PIs/Thummel/Seqs/RepeatMaskedFastas/ExonsMasked");//args[1]);
		File selectList = new File ("/Users/nix/HCI/PIs/Thummel/GenesSeqs/wtEnrichedTop38-129.txt");
		Pattern pat = Pattern.compile("(CG\\d+)-.+");
		int start = 1500;
		int end = 200;
		
		//parse selectList
		LinkedHashSet selectGenes = IO.loadFileIntoLinkedHashSet(selectList);
		
		//parse gene table and split by chromosome, note using interbase numbering
		UCSCGeneModelTableReader genes = new UCSCGeneModelTableReader(geneFile, 0);
		
		//filter?
		if (pat !=null){
			UCSCGeneLine[] all = genes.getGeneLines();
			HashMap map = new HashMap();
			for (int i=0; i< all.length; i++){
				//parse text
				Matcher mat = pat.matcher(all[i].getName());
				if (mat.matches() == false) Misc.printExit("\nProblem matching gene text. "+all[i].getName());
				String trunkName = mat.group(1);
				all[i].setName(trunkName);
				//does it exist?
				if (map.containsKey(trunkName) == false) map.put(trunkName, all[i]);
				else {
					//check which tss is further 5' and save that one
					UCSCGeneLine existing = (UCSCGeneLine) map.get(trunkName);
					if (all[i].getStrand().equals("+")){
						//plus stranded thus want smallest tss
						int newTSS = all[i].getTxStart();
						int existingTSS = existing.getTxStart();
						if (newTSS < existingTSS) map.put(trunkName, all[i]);
					}
					else {
						//negative strand want largest tss
						int newTSS = all[i].getTxEnd();
						int existingTSS = existing.getTxEnd();
						if (newTSS > existingTSS) map.put(trunkName, all[i]);
					}
				}
			}
			System.out.println("Filtered out "+(all.length-map.size()));
			//convert hash to array
			all = new UCSCGeneLine[map.size()];
			int counter =0;
			Iterator it = map.keySet().iterator();
			while (it.hasNext()) all[counter++] = (UCSCGeneLine)map.get(it.next());
			genes.setGeneLines(all);
		}
		genes.splitByChromosome();
		HashMap chromSpecGenes = genes.getChromSpecificGeneLines();
		
		//for each chromosome
		Iterator it = chromSpecGenes.keySet().iterator();
		MultiFastaParser fastaParser = new MultiFastaParser();
		StringBuffer inList = new StringBuffer();
		StringBuffer notInList = new StringBuffer();
		
		while (it.hasNext()){
			//get gene lines
			String chrom = (String) it.next();
			//System.out.println("Processing "+chrom);
			UCSCGeneLine[] sub = (UCSCGeneLine[]) chromSpecGenes.get(chrom);
			
			//parse fasta and get sequence
			fastaParser.parseIt(new File (fastaDir, chrom+".fasta"));
			String seq = fastaParser.getSeqs()[0];
			
			//for each line, get sub sequence
			for (int i=0; i< sub.length; i++){
				int seqStart;
				int seqStop;
				//sense strand?
				if (sub[i].getStrand().equals("+")){
					int tss = sub[i].getTxStart();
					seqStart = tss - start;
					seqStop = tss + end;
				} 
				//nope, antisense
				else {
					int tss = sub[i].getTxEnd();
					seqStart = tss - end;
					seqStop = tss + start;
				}
				//make sure don't exceed bounds
				if (seqStart < 0) seqStart = 0;
				if (seqStop > seq.length()) seqStop = seq.length();
				
				//print unique promoters
				String subSeq = new String (seq.substring(seqStart, seqStop));
				String res = sub[i]+"\t"+seqStart+"\t"+seqStop+"\t"+subSeq+"\n";
				if (selectGenes.contains(sub[i].getName())) inList.append(res);
				else notInList.append(res);
			}
		}
		System.out.println("Genes not in select list:\n"+notInList);
		System.out.println("\n\nGenes in select list:\n"+inList);
	}
}
