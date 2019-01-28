package edu.utah.seq.parsers.rnafold;

import java.io.File;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.bio.parsers.MultiFastaParser;
import util.gen.IO;
import util.gen.Misc;

public class RNACoFold {

	//fields
	private String rnaCoFoldApp;
	private String[] randomUtrs;
	private String[] mRNAs;
	private String[] hsur2;
	private File testTarget;

	private int numRandomTrials = 250;
	private Random random = new Random(0);
	private static Pattern mfePattern = Pattern.compile(".+\\(\\s*(.+)\\)");


	public RNACoFold(String[] args) {
		File minMRNA = new File(args[0]);
		File minHSUR2 = new File(args[1]);
		File randomUTRSeqs = new File(args[2]);
		rnaCoFoldApp = args[3];
		
		testTarget = new File(minMRNA.getParentFile(), "testTargetSeqDelme.fasta");
		testTarget.deleteOnExit();

		MultiFastaParser targetsMFP = new MultiFastaParser(minMRNA);
		mRNAs = targetsMFP.getSeqs();
		
		MultiFastaParser hsur2MFP = new MultiFastaParser(minHSUR2);
		hsur2 = hsur2MFP.getSeqs();
		
		MultiFastaParser utrsMFP = new MultiFastaParser(randomUTRSeqs);
		randomUtrs = utrsMFP.getSeqs();
		
		IO.pl("mRNATarget\tHSUR2\tMFE\tDeltaG\t#RandomMFEHits\tRandomMeanMFE\tpValMFE\t#RandomDeltaGHits\tRandomMeanDeltaG\tpValDeltaG");
		scoreTargets();
	}

	private void scoreTargets() {
		for (int x=0; x< mRNAs.length; x++){
			String mRNA = mRNAs[x];
			String hsur = hsur2[x];
			double[] realMfeDeltaG = score (mRNA, hsur);

			double mfeHits = 0;
			double deltaGHits = 0;
			double mfeTotal = 0;
			double deltaGTotal = 0;

			for (int i=0; i< numRandomTrials; i++){
				String randomMRNA = fetchRandom(mRNA.length());
				double[] randomMfeDeltaG = score (randomMRNA, hsur);
//IO.pl(randomMRNA+" & "+hsur+"\t"+ randomMfeDeltaG[0]+ "\t"+ randomMfeDeltaG[1]);
				if (randomMfeDeltaG[0]<= realMfeDeltaG[0]) mfeHits++;
				if (randomMfeDeltaG[1]<= realMfeDeltaG[1]) deltaGHits++;
				mfeTotal+= randomMfeDeltaG[0];
				deltaGTotal+= randomMfeDeltaG[1];
			}
			double aveMfe = mfeTotal/(double)numRandomTrials;
			double pvalMfe = mfeHits/(double)numRandomTrials;
			double aveDeltaG = deltaGTotal/(double)numRandomTrials;
			double pvalDeltaG = deltaGHits/(double)numRandomTrials;
			
			IO.pl(mRNA+"\t"+hsur+"\t"+ realMfeDeltaG[0]+ "\t"+ realMfeDeltaG[1]+ "\t"+ mfeHits + "\t"+ aveMfe+"\t"+ pvalMfe +"\t"+ deltaGHits+"\t"+aveDeltaG+"\t"+pvalDeltaG);
		
			
		}
		
	}

	private String fetchRandom(int len) {
		//find random utr of sufficient length
		int randomIndex = -1;
		while (true){
			randomIndex = random.nextInt(randomUtrs.length);
			if (randomUtrs[randomIndex].length()> len) break;
		}

		//select random subSeq
		int maxStart = randomUtrs[randomIndex].length() - len;
		int randStart = random.nextInt(maxStart);
		
		return randomUtrs[randomIndex].substring(randStart, len+randStart).toUpperCase();
	}

	private double[] score(String mRNA, String hsur2) {
		//wish I didn't have to write out the file, piping in java?
		IO.writeString(mRNA+ "&"+ hsur2, testTarget);
		String[] cmd = {rnaCoFoldApp, "-a", "-d2", "--noLP", testTarget.toString()};
		String[] res = IO.executeViaProcessBuilder(cmd, false);
		return parseMfeDeltaG (res);
	}

	private double[] parseMfeDeltaG(String[] res) {
		try {

			/*
		 WARNING: vrna_co_pf_probs: numeric instability detected, probability below zero! (sometimes!!!)
		 CUGAGAUGUCGUGUGUGU&ACACUACAUAUUUAUUGUUUAUUUAUACCUGAUAAUGCUGCUUUAAUCACACACGAGGUAUCGCUUCUACUA
		.(((.((.(((((((((.&............(((((((...........)))))))..........))))))))).)).)))......... (-19.70)
		.{((.((.(((((((((.&............{{,{((,...,.......)))))},..........))))))))).)).))}......... [-21.19]
		 frequency of mfe structure in ensemble 0.0896563; delta G binding=-14.88
		Free Energies:
		AB		AA		BB		A		B
		-21.186426	-7.009605	-18.565104	-0.428334	-5.877427
			 */
			int index = 0;
			for (int i=0; i< res.length; i++){
				if (res[i].startsWith("WARNING")) index++;
				else break;
			}
			
			double mfe;
			Matcher mat = mfePattern.matcher(res[index+1]);
			if (mat.matches()) mfe = Double.parseDouble(mat.group(1));
			else throw new Exception();

			String[] s = Misc.PATTERN_EQUALS.split(res[index+3]);
			double deltaG = Double.parseDouble(s[1]);

			return new double[]{mfe, deltaG};
		} catch (Exception e){
			System.out.println("\nERROR: Failed to parse RNAcofold result:");
			Misc.printArray(res);
			//e.printStackTrace();
			System.exit(1);
		}
		return null;
		
	}

	public static void main( String[] args){
		new RNACoFold(args);
	}

}