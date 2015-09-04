package edu.utah.seq.barcodes;

import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Random;

public class TestBarcodeClustering {

	
	/**For testing.*/
	public static void main (String[] args){
		
		//define the thresholds to use
		
		//phred base quality, if it is trusted then set to 20
		int minBaseQuality = 20;
		
		//to process, a barcode must have this number of bases >= the minBaseQuality
		//also, when comparing two barcodes for frac similarity, the number of non N base comparisons must be >= this threshold
		//prob good to set this to 7, 4^7 = 16,384
		int minNumBases = 7; 	
		
		//careful setting this too stringent, 6/7 = 0.857; 5/7 = 0.714
		double minFractionIdentity = 0.75; 
		
		//make a BCE object for each thread and call cluster(sams) multiple times
		BarcodeClusterEngine bce = new BarcodeClusterEngine(minBaseQuality, minNumBases, minFractionIdentity);

		//make some records with increasingly permuted barcodes
		Random random = new Random();
		SAMRecord[] sams = new SAMRecord[5000];
		char[] barcode = "GATCGATCGATCGATC".toCharArray();
		char[] gatc = {'G','A','T','C'};
		for (int i=0; i< sams.length; i++) {
			sams[i] = new SAMRecord(null);
			int seqIndex = random.nextInt(barcode.length);
			int baseIndex = random.nextInt(4);
			barcode[seqIndex] = gatc[baseIndex];
			sams[i].setReadName("NS500690:31:H5YG7BGXX:1:11101:22361:1038:BMF:"+new String(barcode)+"@C4@C54446355///");
		}
		long startTime = System.currentTimeMillis();
		
		//make the call to cluster the records
		//SAMRecord[][] clusters = bce.cluster(sams);
		
		//fetch records that were skipped due to bad quality barcodes
		//these should be deduped using standard procedures (random pick those with same unclipped start? use mate info too?)
		ArrayList<SAMRecord> skippedAlignments = bce.getSkippedRecords();
		
		//print results
		System.out.println(bce);
		
		
		//another test
		System.out.println("\nTest 2");
		SAMRecord[] sams2 = new SAMRecord[7];
		for (int i=0; i< sams2.length; i++) sams2[i] = new SAMRecord(null);
		sams2[0].setReadName("NS500690:31:H5YG7BGXX:1:11101:22361:1038:BMF:GGGGGGGGGGGGGGGG@@@@@@@@@@@@@@@@");
		sams2[1].setReadName("NS500690:31:H5YG7BGXX:1:11101:22361:1038:BMF:CCCCCCCCCCCCCCCC@@@@@@@@@@@@@@@@");
		sams2[2].setReadName("NS500690:31:H5YG7BGXX:1:11101:22361:1038:BMF:TTTTTTTTTTTTTTTT@@@@@@@@@@@@@@@@");
		sams2[3].setReadName("NS500690:31:H5YG7BGXX:1:11101:22361:1038:BMF:AAAAAAAAAAAAAAAA@@@@@@@@@@@@@@@@");
		sams2[4].setReadName("NS500690:31:H5YG7BGXX:1:11101:22361:1038:BMF:TTAAAAAAAAAAAAAA@@@@@@@@@@@@@@@@");
		sams2[5].setReadName("NS500690:31:H5YG7BGXX:1:11101:22361:1038:BMF:GATCGATCGATCGATC@@@@@@@@@@@@@@@@");
		sams2[6].setReadName("NS500690:31:H5YG7BGXX:1:11101:22361:1038:BMF:AAAAAATCGATCGATC@@@@@44444444444");//kill the base quality
		//bce.cluster(sams2);
		System.out.println(bce);
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime));
		System.out.println("Done! "+diffTime+" millsec\n");
	}
	
	
}
