package util.bio.calc;
import java.text.*;

import util.bio.seq.*;
import util.gen.Num;
/**
 *Nucleic Acid parameter calculator and converter static methods.
*/
public class NucleicAcid {

	/**coverts ng/ul to uM given a DNA sequence*/
	public static double convertToMicroMolarDNA(String seq, double ngPerUl){
		int[] gatcn = Seq.countBases(seq);		
		double mw = calculateMolecularWtDNA(gatcn);
		return ngPerUl/ (mw/1000);
	}
	
	/**coverts uM to ng/ul given a DNA sequence*/
	public static double convertToNgPerUlDNA(String seq, double uM){
		int[] gatcn = Seq.countBases(seq);	
		double mw = calculateMolecularWtDNA(gatcn);		
		return uM*(mw/1000);
	}
	
	/**coverts ng/ul to uM given an RNA sequence*/
	public static double convertToMicroMolarRNA(String seq, double ngPerUl){
		String rna = seq.toLowerCase().replaceAll("u","t");
		int[] gatcn = Seq.countBases(rna);		
		double mw = calculateMolecularWtRNA(gatcn);
		return ngPerUl/ (mw/1000);
	}
	
	/**coverts uM to ng/ul given an RNA sequence*/
	public static double convertToNgPerUlRNA(String seq, double uM){
		String rna = seq.toLowerCase().replaceAll("u","t");
		int[] gatcn = Seq.countBases(rna);	
		double mw = calculateMolecularWtRNA(gatcn);		
		return uM*(mw/1000);
	}
	
	
	/**Calculates molecular weight of a DNA sequence given an int[] of gatcn's, no 5' phosphate*/
	public static double calculateMolecularWtDNA(int[] gatcn){
		return ((double)gatcn[0])*329.21 +((double)gatcn[1])*313.21+ 
		((double)gatcn[2])*304.2+ ((double)gatcn[3])*289.18+ ((double)gatcn[4])*308.95;
	}
	
	/**Calculates molecular weight of a RNA sequence given an int[] of gaucn's, no 5' phosphate*/
	public static double calculateMolecularWtRNA(int[] gaucn){
		return ((double)gaucn[0])*345.21 +((double)gaucn[1])*329.21+ 
		((double)gaucn[2])*306.17+ ((double)gaucn[3])*305.18+ ((double)gaucn[4])*321.44;
	}
	
	/**Calculates a GC concentration, case insensitive.*/
	public static double calculateFractionGC(String DNA){
		int[] gatcn = Seq.countBases(DNA);
		return ((double)(gatcn[0]+gatcn[3]))/((double)(DNA.length()));	
	}
	
	/**Takes a double in percent returns a String with one fraction 0.85842 -> 85.8% */
	public static String formatPercentOneFraction(double num){
		NumberFormat f = NumberFormat.getPercentInstance();
		f.setMaximumFractionDigits(1);
		return f.format(num);
	}
	
	/**Takes a double returns a String with one fractional digit  9.54443 -> 9.5*/
	public static String formatNumberOneFraction(double num){
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(1);
		return f.format(num);
		}
	
	/**Calculate basic Oligo Tm*/
	public static double calculateBasicOligoTm(String seq){
		int[] gatcn = Seq.countBases(seq);
		double length = seq.length();
		double numGC = gatcn[0]+gatcn[3];
		double numAT = gatcn[1]+gatcn[2];
		if (length<14) return (numAT*2) + (numGC*4);
		return 64.9+41*(((numGC - 16.4)/length));
	}
	
	/**Calculates a nearest neighbor Oligo Tm at 50nM oligo, 50mM Salt, also corrects case
	 * for real method.*/
	public static double calcDefaultNearestNeighborTm(String seq){
		return nearestNeighborOligoTm(seq.toLowerCase(), 0.00000005, 0.05); //for big oligos
	}	
	
	/**Calculate salt adjusted Oligo Tm*/
	public static double calculateSaltAdjustedOligoTm(String seq){
		double mono = 0.05; //monvalent salts in Moles, 50mM is standard for PCR
		int[] gatcn = Seq.countBases(seq);
		double length = seq.length();
		double numGC = gatcn[0]+gatcn[3];
		double numAT = gatcn[1]+gatcn[2];
		double gc = numGC/length;
		double monoCalc = 16.6 * (Math.log(mono)/Math.log(10));
		if (length<14) return (numAT*2) + (numGC*4); //assuming 0.05M Na+ 
		if (length<51) return 100.5 + (monoCalc) + (0.41 * gc * 100) - (820/ length);
		return 81.5 + (monoCalc) + (0.41 * gc * 100) - (500/ length); //for big oligos
	}
	
	/**Calculates nearest neighbor tm given a DNA oligonucleotide.
	 * Concentrations are in moles, suggest 0.00000005 for oligo (50nM) and 0.05 salt (50mM).
	 * Lower case is faster, only gatcGATC, nothing else.
	 * Tm = {(DH-3.4)/(DS+1.9872ln(1/[primer])} + 16.6*log10([Na+]) Ð273.15*/
	public static double nearestNeighborOligoTm(String seqLowerCase, double oligoConcentration, double saltConcentration){
		//fetch total deltas
		double totalDeltaH = getTotalDelta(dH, seqLowerCase);
		double totalDeltaS = getTotalDelta(dS, seqLowerCase);
		//calculate top
		double top= (totalDeltaH - 3.4)*1000;
		//calculate bottom
		double bottom = totalDeltaS + 1.9872 * Math.log(1.0/oligoConcentration) ;
		return (top/bottom) + 16.6* (Num.log10(saltConcentration)) - 273.15;
	}
	/**Returns the total delta given a lower case DNA sequence. 
	 * For tmNN().*/
	public static double getTotalDelta(double[][] delta, String seq){
		double total = 0;
		int numPairs = seq.length()-1;
		for (int i=0; i<numPairs; i++){
			int firstBase = baseNumber(seq.charAt(i));
			int secondBase = baseNumber(seq.charAt(i+1));
			total+= delta[firstBase][secondBase];
		}
		return total;
	}
	/**Returns 0 for a, 1 for c, 2 for g, 3 for t, -1 for everything else.
	 * For tmNN().*/
	public static int baseNumber(char b){
		if (b=='a') return 0;
		if (b=='c') return 1;
		if (b=='g') return 2;
		if (b=='t') return 3;
		if (b=='A') return 0;
		if (b=='C') return 1;
		if (b=='G') return 2;
		if (b=='T') return 3;
		System.out.println("Warning, bad char found in baseNumber() in the NucleicAcid.java class ->"+b);
		return -1;
	}
	/**delta H values for dinucleotide pairs.
	 * For tmNN().*/
	public static final double[][] dH = {
			//second base
			//A,  C,   G,   T,	    //First Base
			{8.0, 9.4, 6.6, 5.6},		//A
			{8.2, 10.9, 11.8, 6.6},	//C
			{8.8, 10.5, 10.9, 9.4},	//G
			{6.6, 8.8, 8.2, 8.0}		//T
	};
	/**delta S values for dinucleotide pairs.
	 * For tmNN().*/
	public static final double[][] dS = {
			//second base
			//A,    C,     G,   T,	//First Base
			{21.9, 25.5, 16.4, 15.2},	//A
			{21.0, 28.4, 29.0, 16.4},	//C
			{23.5, 26.4, 28.4, 25.5},	//G
			{18.4, 23.5, 21.0, 21.9}	//T
	};
	
	
}
