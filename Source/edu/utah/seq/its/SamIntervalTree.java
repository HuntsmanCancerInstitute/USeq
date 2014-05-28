package edu.utah.seq.its;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import util.bio.annotation.ExonIntron;
import util.bio.parsers.UCSCGeneLine;
import util.bio.parsers.UCSCGeneModelTableReader;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class SamIntervalTree {
	
	//fields
	File alignmentFile;
	SamReader samReader;
	
	public SamIntervalTree(File sortedSamBamFile) throws Exception {
		this.alignmentFile = sortedSamBamFile;
		
		SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		samReader = factory.open(alignmentFile);
		
		//has an index?
		if (samReader.hasIndex() == false){
			throw new Exception("Error: the following sam/bam file does not appear to have an associated index file, aborting -> "+alignmentFile);
		}
		
	}
	
	public void countHits(HashMap<String, UCSCGeneLine[]> regionsToScore) {
		
		//for each chromosome of genes
		for (String chrom : regionsToScore.keySet()){
			System.out.println("Counting "+chrom);
			UCSCGeneLine[] genes = regionsToScore.get(chrom);
			//find start and stop
			int[] minMax = UCSCGeneLine.findMinMax(genes);
			
			//make interval tree for this chromosome
			IntervalTree<String> tree = makeIntervalTree(chrom, minMax[0], minMax[1]);
			if (tree == null) {
				System.out.println("\tNo alignments to chromosome! Skipping.");
				continue;
			}
			
			//search tree for hits to regions
			search (tree, genes);
		}
		
		
	}

	private IntervalTree<String> makeIntervalTree(String chrom, int start, int stop) {
		//check to see if in header
		if (samReader.getFileHeader().getSequence(chrom) == null) {
			System.out.println("No alignments for "+chrom);
			return null;
		}
		
		//get interator
		SAMRecordIterator iterator = samReader.queryOverlapping(chrom, start, stop);

		//make container for intervals
		List<Interval<String>> intervals = new ArrayList<Interval<String>>();

		//for each sam record
		while (iterator.hasNext()){
			SAMRecord sam = iterator.next();
			String fragmentName = new String(sam.getReadName());
			
			//for each alignment block
			List<AlignmentBlock> blocks = sam.getAlignmentBlocks();
			for (AlignmentBlock b: blocks){
				int s = b.getReferenceStart() -1;
				int e = s + b.getLength();
				intervals.add(new Interval<String>(s, e, fragmentName));
			}
		}
		iterator.close();
		
		if (intervals.size() == 0) return null;
		return new IntervalTree<String>(intervals, false);
	}

	private void search(IntervalTree<String> tree, UCSCGeneLine[] genes) {
		HashSet<String> allNames = new HashSet<String>();
		HashSet<String> exonNames = new HashSet<String>();
		ArrayList<String> alignNames;
		//for each gene
		for (UCSCGeneLine gene: genes){
			allNames.clear();
			//for each exon
			ExonIntron[] ei = gene.getExons();
			int[] exonCounts = new int[ei.length];
			
			for (int i=0; i< ei.length; i++){
				exonNames.clear();
				alignNames = tree.search(ei[i].getStart(), ei[i].getEnd());
				if (alignNames.size() !=0) {
					exonNames.addAll(alignNames);
					allNames.addAll(exonNames);
				}
				exonCounts[i] = exonNames.size();
			}
			
		}
	}

	public static void main (String[] args){
		File bam = new File("/Users/u0028003/Desktop/dupFreeJessieSams_STP13MQ1N120A.bam");
		File ucscFile = new File("/Users/u0028003/HCI/Annotation/Hg19/Dec30_2013Hg19Ens/hg19EnsDoOver_Merged.ucsc.zip");
		try {
			UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader(ucscFile, 0);
			HashMap<String, UCSCGeneLine[]> chromGenes = reader.getChromSpecificGeneLines();
			
			long start = System.currentTimeMillis();
			
			SamIntervalTree sit = new SamIntervalTree(bam);
			sit.countHits(chromGenes);
			
			long diff = System.currentTimeMillis() - start;
			
			System.out.println("OneRep "+diff);
			
			/*
			HashMap<String, ExonIntron[]> map = new HashMap<String, ExonIntron[]>();
			ExonIntron[] exons = new ExonIntron[3];
			exons[0] = new ExonIntron(31767409, 31767473);
			exons[1] = new ExonIntron(31767400,31767466);
			exons[2] = new ExonIntron(31769045,31769046);
			map.put("chr20", exons);
			
			ExonIntron[] exons2 = new ExonIntron[1];
			exons2[0] = new ExonIntron(25334176, 25334178);
			
			map.put("chr22", exons2);
			
			SamIntervalTree sit = new SamIntervalTree(bam);
			sit.countHits(map);
			*/
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
}
