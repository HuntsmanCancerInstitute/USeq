package edu.utah.seq.parsers.rnafold;

import java.io.File;
import java.io.IOException;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.bio.parsers.MultiFastaParser;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;


/*
 * Used to run this app:
 * /Users/u0028003/HCI/Labs/Cazalla_Demian/BootStrappingMFold/targets.fasta 
/Users/u0028003/HCI/Labs/Cazalla_Demian/BootStrappingMFold/HSUR2.fasta 
/Users/u0028003/HCI/Labs/Cazalla_Demian/MotifSearching/Reference/calJac3_UTRsDeDuped_Rand10K.fasta.gz 
/usr/local/bin/RNAcofold
 * 
 u0028003$ head -n 6 targets.fasta
>QueryID1
GAGAGAGGTCGATATCGAAGATAATCCGATTTTGGAAGCAATGATTGTTTAAAGATACAAAAAATC
>QueryID2
TGGCAGAAGACACCAGAGCAGATGCAG
>QueryID3
GTAGAAGCATTCCTTCTTTGATAATGTTAAATTTGTAAGTTTCAGGTGACATGTGAAA

 * u0028003$ head -n 6 HSUR2.fasta 
>QueryID1
UAUUUAUUGUUUAUUUAUACCUGAUAAUGCUGCUUUAAUCACACACGAGGUAUCGCUUCUACUA
>QueryID2
CUGAUAAUGCUGCUUUAAUCACACACGAGGUAUCGCUUCUACUA
>QueryID3
GUUUAUUUAUACCUGAUAAUGCUGCUUUAAUCACACACGAGGUAUCGCUUCUAC

calJac3_UTRsDeDuped_Rand10K.fasta.gz
? random 10K UTRs?
 */


public class RNACoFold {

	//fields
	private String rnaCoFoldApp;
	private String[] randomUtrs;
	private String[] mRNAs;
	private String[] hsur2;
	private File testTarget;

	private int numRandomTrials = 1000;
	private Random random = new Random(0);
	private static Pattern mfePattern = Pattern.compile(".+ensemble ([\\d\\.]+);.+");


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

		IO.pl("mRNATarget\tHSUR\tMFE\tDeltaG\t#RandomMFEHits\tRandomMeanMFE\tpValMFE\t#RandomDeltaGHits\tRandomMeanDeltaG\tpValDeltaG");
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
		
		# 2.7.2
		[u0028003@redwood11:Delme]$ cat results.txt 
		>first
		CACCCAAUCUGCCCUUUUGCCCCCCCAGCCCUG&UAGGUACUGGGUGUAAAUAUGAUGACCGGUA
		.(((..(((......(((((...(((((.((((&))))..))))).)))))...)))....))). (-16.00)
		.(((..,((......(((((...(((((.((({&})))..))))).)))))...,,,},..))). [-17.29]
		frequency of mfe structure in ensemble 0.122576; delta G binding=-10.35
		Free Energies:
		AB              AA              BB              A               B
		-17.293675      -4.252051       -17.540123      -0.499205       -6.440962
		
			 */
			//updated for 2.7.2
			int index = -1;
			for (int i=0; i< res.length; i++){
				res[i] = res[i].trim();
				if (res[i].startsWith("frequency of mfe structure")) {
					index = i;
					break;
				}
			}
			if (index ==-1) throw new IOException("ERROR: failed to find the frequency line in "+Misc.stringArrayToString(res, "\n"));
			
			//IO.pl("\n"+Misc.stringArrayToString(res, "\n"));
			//IO.pl("'"+res[index-2]+"'");
			
			//MFE is the number on the end of the first structure, -16.00, sometimes there is a space before the - in 16.00
			String[] fs = Misc.WHITESPACE.split(res[index-2]);
			
			String mfeString = null;
			if (fs.length==2) mfeString = fs[1].substring(1, fs[1].length()-1);
			else if (fs.length == 3) mfeString = fs[2].substring(0, fs[2].length()-1);
			else throw new IOException("ERROR: failed to parse the mfe from "+res[index-2]);
			double mfe = Double.parseDouble(mfeString);

			//deltaG is the last number in the freq line, -10.35
			String[] s = Misc.PATTERN_EQUALS.split(res[index]);
			double deltaG = Double.parseDouble(s[1]);

			return new double[]{mfe, deltaG};
		} catch (Exception e){
			System.out.println("\nERROR: Failed to parse RNAcofold result:");
			Misc.printArray(res);
			e.printStackTrace();
			System.exit(1);
		}
		return null;
		
	}

	public static void main( String[] args){
		new RNACoFold(args);
	}

}