package util.bio.parsers;
import java.io.*;
import java.util.*;

import edu.utah.seq.useq.data.Region;



public class XMerReadBiasParser {
	//fields
	int readSize = 26;
	int buffer = 200;
	double multiplier = 0;
	
	public XMerReadBiasParser( String[] args){
		
		File geneFile = new File (args[0]);
		File seqBedFile = new File (args[1]);
		String strand = args[2];
		multiplier = Double.parseDouble(args[3]);
		
		HashMap<String, UCSCGeneLine[]> genes = new UCSCGeneModelTableReader(geneFile,0).getChromSpecificGeneLines();
		HashMap<String, Region[]> seqs = Region.parseStartStops(seqBedFile, 0, 0, 0);
		
		Iterator<String> it = genes.keySet().iterator();
		
		System.out.println("Chrom\tStartGene\tStopGene\tMiddleEndoCluster\tNumXmersLessMiddle\tNumXmersMoreMiddle\tName\tStrand="+strand+"; Multiplier="+multiplier);
		//for each chromosome
		int totalNumLess = 0;
		int totalNumMore = 0;
		Random random = new Random(0);
		
		while (it.hasNext()){
			String chrom = it.next();
			UCSCGeneLine[] g = genes.get(chrom);
			Region[] s = seqs.get(chrom);
			if (s == null) continue;
			ArrayList<Region> al = new ArrayList<Region>();
			LinkedHashSet<String> parsedReads = new LinkedHashSet<String>();
			//for each gene
			for (int i=0; i< g.length; i++){
				
				//check strand 
				if (g[i].getStrand().equals(strand) == false) continue;
				int start = g[i].getTxStart() - buffer;
				int stop = g[i].getTxEnd() + buffer;
				al.clear();
				
				LinkedHashSet<String> parsedReadsForGene = new LinkedHashSet<String>();
				//for each sequence
				for (int j =0; j< s.length; j++){
					if (s[j].intersects(start, stop)) {
						al.add(s[j]);
						parsedReadsForGene.add(s[j].toString());
					}
				}
				//enough reads?
				if (al.size()<3) continue;
				Region[] reads = new Region[al.size()];
				al.toArray(reads);
				
				//calc middle
				double blockStart = reads[0].getStart();
				double blockEnd = reads[reads.length-1].getStop();
				double tenPercent = (blockEnd - blockStart) /10;
				
				//50
				int divider = (int)Math.round(((blockEnd-blockStart) / 2) +tenPercent* multiplier) + reads[0].getStart();
				
				//score Xmers
				int numLess =0;
				int numMore =0;
				for (int j=0; j<reads.length; j++){
					//if (reads[j].getLength() >19 && reads[j].getLength() <23 ) continue;
					//if (reads[j].getLength() != readSize ) continue;
					int readMiddle = reads[j].getMiddle();
					if (readMiddle < divider) numLess++;
					else if (readMiddle > divider) numMore++;
					else {
						double r = random.nextDouble();
						if (r <= 0.5) numLess++;
						else numMore++;
					}
				}
				
				//check if this is a double sampling of the same genes
				boolean duplicate = false;
				Iterator<String> itGene = parsedReadsForGene.iterator();
				while (itGene.hasNext()){
					String x = itGene.next();
					if (parsedReads.contains(x)) {
						duplicate = true;
						break;
					}
				}
				parsedReads.addAll(parsedReadsForGene);
				
				if (duplicate) continue;
				
				String dups = "";
				if (duplicate) dups = "\tduplicate";
				if (numLess !=0 || numMore !=0) {
				System.out.println(chrom+"\t"+start+"\t"+stop+"\t"+divider+"\t"+numLess+"\t"+numMore+"\t"+g[i].getName()+ dups);
					totalNumLess += numLess;
					totalNumMore += numMore;
				}
				
				
			}
		}
		
		System.out.println("Totals\t"+totalNumLess+"\t"+totalNumMore);
		
	}

	public static void main(String[] args) {
		new XMerReadBiasParser(args);

	}

}
