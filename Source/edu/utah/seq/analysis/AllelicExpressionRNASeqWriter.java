package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import htsjdk.samtools.*;
import htsjdk.tribble.readers.TabixReader;
import util.bio.annotation.Bed;
import util.bio.parsers.UCSCGeneLine;
import util.bio.parsers.UCSCGeneModelTableReader;
import util.gen.*;
import edu.utah.seq.analysis.multi.PairedCondition;
import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.data.sam.SamLayoutForMutation;
import edu.utah.seq.vcf.VCFLookUp;
import edu.utah.seq.vcf.VCFParser;
import edu.utah.seq.vcf.VCFRecord;

/** 
 * @author Nix
 * */
public class AllelicExpressionRNASeqWriter {

	//user defined fields
	private File vcfFile;
	private File bamFile;
	private File resultsFile;
	private File tabixBedFile;
	private int minimumReadCoverage = 10;
	private int minimumBaseQuality = 20;
	private int minimumAltCoverage = 2;
	private int minimumRefCoverage = 2;
	private double minimumAF = 0.05;
	private double maximumAF = 0.95;
	private int readLength = 50;
	private boolean printHeader = true;
	
	//internal
	private HashMap<String, VCFLookUp> chrVcf = null;
	private ArrayList<VCFCount> passingRecords = new ArrayList<VCFCount>();
	private HashMap<String,ArrayList<VCFCount>> geneVcf = new HashMap<String,ArrayList<VCFCount>>();
	private ArrayList<VCFCount> noGeneHit = new ArrayList<VCFCount>();
	
	
	//constructor
	/**Stand alone.*/
	public AllelicExpressionRNASeqWriter(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		//load snp info via bed file
		System.out.println("Loading variant records from vcf file... ");
		loadVcfRecords();

		//walk through variant data looking for snvs with sufficient coverage
		System.out.print("\nReading vcf data and extracting alignments");
		walkVariantData();
		
		//intersect results with gene table
		System.out.println("\nAdding gene annotation...");
		annotateWithGenes();

		//write out summary
		System.out.println("\nWriting GeneiASE table...");
		writeOutSummaryTable();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Num.formatNumber(diffTime, 2)+" minutes\n");
	}


	private void annotateWithGenes() {
		try {
			TabixReader tr = new TabixReader( tabixBedFile.toString() );


			//for each good record, fetch bed records where name is the gene name
			for (VCFCount vc: passingRecords){
				int pos = vc.vcf.getPosition() - readLength;
				int end = vc.vcf.getPosition() + readLength;
				String q = vc.vcf.getChromosome()+":"+pos+"-"+end;
				TabixReader.Iterator it = tr.query(q);
				String bedLine;
				ArrayList<VCFCount> al = null; 
				//might be more than one
				while ((bedLine = it.next()) != null) {
					//chr begin end name score strand
					String[] tokens = Misc.TAB.split(bedLine);
					String geneName = tokens[3];
					al = geneVcf.get(geneName);
					if (al == null) {
						al = new ArrayList<VCFCount>();
						geneVcf.put(geneName, al);
					}
					al.add(vc);
				}
				//no hits?
				if (al == null) noGeneHit.add(vc);
			}
			tr.close();

		} catch (IOException e) {
			e.printStackTrace();
		} 
	}



	private void writeOutSummaryTable() {
		try {
			//make gzipper
			Gzipper out = new Gzipper(resultsFile);
			//header
			if (printHeader) out.println("gene\tsnp.id\talt.dp\tref.dp");
			//use to keep the snp.id column unique even though some are shared between genes
			int snpId = 0;
			//for each gene
			for (String geneName: geneVcf.keySet()){
				ArrayList<VCFCount> hits = geneVcf.get(geneName);
				//for each hit
				for (VCFCount vc: hits){
					out.print(geneName); out.print("\t");
					out.print(vc.vcf.getChrPosRefAlt()); out.print("_"); out.print(snpId++);out.print("\t");
					out.print(vc.alt); out.print("\t");
					out.print(vc.ref); out.println("\t");
				}
			}
			//for those with no associated gene
			for (VCFCount vc: noGeneHit){
				String noName = vc.vcf.getChrPosRefAlt();
				out.print(noName); out.print("\t");
				out.print(noName); out.print("_"); out.print(snpId++); out.print("\t");
				out.print(vc.alt); out.print("\t");
				out.print(vc.ref); out.println("\t");
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	
	



	private void walkVariantData() {
		try {
			//open sam reader
			SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			SamReader samReader = factory.open(bamFile);
			
			HashMap<String, Character> counts = new HashMap<String, Character>();

			for (String chr: chrVcf.keySet()) {
				System.out.print(".");
				VCFLookUp lookUp = chrVcf.get(chr);
				VCFRecord[] rec = lookUp.getVcfRecord();
				
				//for each record
				for (VCFRecord vcf: rec){
					//System.out.println(vcf.getOriginalRecord());
					//is it a snv?
					if (vcf.isSNP() == false) {
						//System.out.println("Skip: Not snp");
						continue;
					}
					
					//any nearby variants?
					VCFRecord[] n = lookUp.fetchVCFRecords(vcf.getPosition()- readLength, vcf.getPosition()+ readLength);
					if (n.length > 1) {
						//System.out.println("Skip: other vars");
						continue;
					}
					
					//fetch alignments
					counts.clear();
					SAMRecordIterator it = samReader.queryOverlapping(vcf.getChromosome(), vcf.getPosition(), vcf.getPosition()+1);
					while (it.hasNext()){
						SAMRecord sam = it.next();
						SamAlignment sa = new SamAlignment(sam.getSAMString().trim(), true);
						SamLayoutForMutation layout = new SamLayoutForMutation(sa);
						int index = layout.findBaseIndex(vcf.getPosition());
						if (index == -1) continue;
						//good score?
						if (layout.getQual()[index] < minimumBaseQuality) {
							//System.out.println("Bad base!");
							continue;
						}
						//a match?
						if (layout.getCall()[index] != 'M'){
							//System.out.println("Not a M!");
							continue;
						}
						//add to hash
						char base = layout.getSeq()[index];

						counts.put(sa.getName(), new Character(base));
					}
					it.close();

					if (counts.size() == 0) {
						//System.out.println("Skip: No counts");
						continue;
					}

					//check count
					int[] gatcn = countHash(counts);
					int[] refAltCounts = fetchRefAltCounts(vcf, gatcn);
					
					// min counts?
					if ( (refAltCounts[0] + refAltCounts[1]) < minimumReadCoverage){
						//System.out.println("Skip: Fail min cov");
						continue;
					}
					if (refAltCounts[0] < minimumRefCoverage) {
						//System.out.println("Skip: Fail min ref cov");
						continue;
					}
					if (refAltCounts[1] < minimumAltCoverage) {
						//System.out.println("Skip: Fail min alt cov");
						continue;
					}
					
					// af?
					double af = (double) refAltCounts[1]/ (double)(refAltCounts[0] + refAltCounts[1]);
					if (af < minimumAF || af > maximumAF) {
						//System.out.println("Skip: Fail AF "+af);
						continue;
					}
					
					//save it!
					passingRecords.add(new VCFCount(refAltCounts[0], refAltCounts[1], vcf));
					
				}
			}
			System.out.println();
			samReader.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private class VCFCount {
		int ref;
		int alt;
		VCFRecord vcf;
		
		public VCFCount(int ref, int alt, VCFRecord vcf){
			this.ref = ref;
			this.alt = alt;
			this.vcf = vcf;
		}
	}

	private int[] fetchRefAltCounts(VCFRecord vcf, int[] gatcn) {
		int refC = fetchCounts (vcf.getReference().toUpperCase(), gatcn);
		int altC = fetchCounts (vcf.getAlternate()[0].toUpperCase(), gatcn);
		return new int[]{refC, altC};
	}

	private int fetchCounts(String ref, int[] gatcn) {
		if (ref.equals("G")) return gatcn[0];
		if (ref.equals("A")) return gatcn[1];
		if (ref.equals("T")) return gatcn[2];
		if (ref.equals("C")) return gatcn[3];
		return gatcn[4];
	}

	private int[] countHash(HashMap<String, Character> counts) {
		int[] gatcn = new int[5];
		for (Character c : counts.values()){
			if (c == 'G') gatcn[0]++;
			else if (c == 'A') gatcn[1]++;
			else if (c == 'T') gatcn[2]++;
			else if (c == 'C') gatcn[3]++;
			else if (c == 'N') gatcn[4]++;
			else Misc.printErrAndExit("\nError: unrecognized base '"+c+"' from "+counts);
		}
		return gatcn;
	}



	private void loadVcfRecords() {
		VCFParser p = new VCFParser (vcfFile, true, false, false);
		chrVcf = p.getChromosomeVCFRecords();
	}

	
	public static void main(String[] args) {

		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AllelicExpressionRNASeqWriter(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bamFile = new File(args[++i]); break;
					case 'v': vcfFile = new File(args[++i]); break;
					case 't': tabixBedFile = new File(args[++i]); break;
					case 'o': resultsFile = new File(args[++i]); break;
					case 'l': readLength = Integer.parseInt(args[++i]); break;
					case 'c': minimumReadCoverage = Integer.parseInt(args[++i]); break;
					case 'a': minimumAltCoverage = Integer.parseInt(args[++i]); break;
					case 'r': minimumRefCoverage = Integer.parseInt(args[++i]); break;
					case 'q': minimumBaseQuality = Integer.parseInt(args[++i]); break;
					case 'm': minimumAF = Double.parseDouble(args[++i]); break;
					case 'x': maximumAF = Double.parseDouble(args[++i]); break;
					case 'p': printHeader = false; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		if (bamFile == null || bamFile.canRead() == false) Misc.printErrAndExit("\nError: can't find your bam alignment file?\n");
		if (vcfFile == null || vcfFile.canRead() == false) Misc.printErrAndExit("\nError: can't find your vcf variant file?\n");
		if (tabixBedFile == null || tabixBedFile.canRead() == false) Misc.printErrAndExit("\nError: can't find your tabixed bed file?\n");
		if (resultsFile == null) Misc.printErrAndExit("\nError: please provide a results output file?\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                      Allelic Expression RNASeq Writer :  Sept 2016               **\n" +
				"**************************************************************************************\n" +
				"Application for parsing count data for downstream Allele Specific Gene Expression\n"+
				"detection, e.g. GeneiASE. Avoids snvs with vars within the read length, skips INDELs.\n\n"+

				"Required Arguments:\n"+
				"-b Bam file with associated index from an RNASeq experiment after filtering for\n"+
				"      allelic alignment bias.\n"+
				"-v Vcf file containing snvs to use in extracting alignment counts from the bam. These\n"+
				"      will be filtered using the params below before saving.\n"+
				"-t Tabix gz indexed bed file of exons where the name column is the gene name, see\n"+
				"      ExportExons and https://github.com/samtools/htslib\n"+
				"-o Output file.\n"+
				
				"\nDefault Arguments:\n"+
				"-l Read length, defaults to 50\n"+
				"-c Minimum alignment depth for quality bases, defaults to 10\n"+
				"-q Minimum base quality, defaults to 20\n"+
				"-a Minimum Alt alignment depth, defaults to 2\n"+
				"-r Minimum Ref alignment depth, defaults to 2\n"+
				"-m Minimum allele frequency, defaults to 0.05\n"+
				"-x Maximum allele frequency, defaults to 0.95\n"+
				"-p Don't print header\n"+

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/AllelicExpressionRNASeqWriter -b proc.bam\n"+
				"-v lofreq.vcf.gz -t ~/Anno/b37EnsGenes7Sept2016_Exons.bed.gz -o forGeneiASE.txt.gz\n\n" +

				"**************************************************************************************\n");

	}
}
