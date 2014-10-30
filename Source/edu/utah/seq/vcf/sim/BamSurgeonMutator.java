package edu.utah.seq.vcf.sim;

import java.io.File;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.bio.annotation.Bed;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Creates snvs and indels for the BAMSurgeon.*/
public class BamSurgeonMutator {
	//fields
	private File bedFile;
	private char[] gatc = {'G','A','T','C'};
	private Gzipper mutSNVs;
	private Gzipper mutINDELs;
	private Random random = new Random(9);
	private int numberSNVs = 0;
	private int numberINDELs = 0;
	private int maxIndelSize = 6;
	private double allelicRatio = 0.5;
	private double minAllelicRatio = 0.1;

	public BamSurgeonMutator (String[] args){
		long startTime = System.currentTimeMillis();

		processArgs(args);

		try {
			//load bed and randomize
			Bed[] regions = Bed.parseFile(bedFile, 0,0);
			int totalToMutate = numberSNVs+numberINDELs;
			if (regions.length <= totalToMutate) Misc.printErrAndExit("\nError: More mutations requested than number of regions? Aborting!\n");
			Misc.randomize(regions, 9);
			
			//make gzippers for mutant info, fileName_allelicRatio_INDELS_number.txt.gz
			String name = Misc.removeExtension(bedFile.getName());
			name = name+"_"+(Num.formatNumber(allelicRatio, 2));
			if (numberSNVs !=0) mutSNVs = new Gzipper(new File( bedFile.getParentFile(), name+"SNVs"+numberSNVs+".txt.gz"));
			if (numberINDELs !=0) mutINDELs = new Gzipper(new File( bedFile.getParentFile(), name+"INDELs"+numberINDELs+".txt.gz"));
			
			//for each region, mutate
			System.out.println("Processing:");
			for (int i=0; i< numberSNVs; i++) mutateSNV(regions[i]);
			for (int i=numberSNVs; i< totalToMutate; i++) mutateINDEL(regions[i]);
			
			//close gzipper
			if (numberSNVs !=0) mutSNVs.close();
			if (numberINDELs !=0) mutINDELs.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		} 

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	private void mutateSNV(Bed region) throws Exception{
		int size = region.getLength();
		int pos = random.nextInt(size);
		int start = region.getStart()+pos;
		region.setStart(start);
		if (allelicRatio == 0) region.setScore(getRandomAllelicRatio());
		else region.setScore(allelicRatio);
		region.setStop(start+1);
		mutSNVs.println(region.toStringNoStrandNoName());
		
	}
	
	private void mutateINDEL(Bed region) throws Exception{
		int size = region.getLength();
		int pos = random.nextInt(size);
		int start = region.getStart()+pos;
		region.setStart(start);
		if (allelicRatio == 0) region.setScore(getRandomAllelicRatio());
		else region.setScore(allelicRatio);
		int lengthOfIndel = random.nextInt(maxIndelSize) + 1;

		//insertion 
		if (random.nextBoolean()){
			String seq = this.fetchRandomSeq(lengthOfIndel);
			region.setStop(start+1);
			region.setName("INS\t"+seq);
			mutINDELs.println(region.toStringScoreName());
		}
		//deletion
		else {
			int stop = start + lengthOfIndel;
			region.setStop(stop);
			region.setName("DEL");
			mutINDELs.println(region.toStringScoreName());
		}
	}

	
	private double getRandomAllelicRatio(){
		while (true) {
			double n = random.nextDouble();
			if (n >= minAllelicRatio) return n;
		}
		
	}
	
	private String fetchRandomSeq(int length){
		char[] s = new char[length];
		for (int i=0; i< length; i++) s[i] = fetchRandomBase();
		return new String(s);
	}
	
	private char fetchRandomBase(){
		int index = random.nextInt(4);
		return gatc[index];
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BamSurgeonMutator(args);
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
					case 'r': bedFile = new File(args[++i]); break;
					case 's': numberSNVs = Integer.parseInt(args[++i]); break;
					case 'i': numberINDELs = Integer.parseInt(args[++i]); break;
					case 'm': maxIndelSize = Integer.parseInt(args[++i]); break;
					case 'a': allelicRatio = Double.parseDouble(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		int numMuts = numberSNVs + numberINDELs;
		if (bedFile == null || bedFile.canRead()== false) Misc.printErrAndExit("\nPlease enter a bed file of regions to target for mutation.\n");
		if (numMuts == 0) Misc.printErrAndExit("\nPlease enter the number of SNPs and INDELs you would like to create.\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             BamSurgeon Mutator  : Oct 2014                       **\n" +
				"**************************************************************************************\n" +
				"Generates bed files of random SNVs and INDELs for BamSurgeon. Only one variant is made\n"+
				"per input region.\n\n" +

				"Required:\n"+
				"-r Full path bed file containing regions in which to create mutations.\n"+
				"-s Number of SNVs to make, set to zero for just INDELS.\n\n"+
				"-i Number of INDELs to make. This plus -s shouldn't exceed the number of regions.\n\n"+

				"Optional:\n" +
				"-m Max indel size, defaults to 6\n"+
				"-a Allelic ratio, defaults to 0.5. Set to zero for random from 0.1 - 1\n"+
				"\n"+

				"Example: java -Xmx10G -jar pathTo/USeq/Apps/RandomMutationGenerator \n\n" +

				"**************************************************************************************\n");

	}
}
