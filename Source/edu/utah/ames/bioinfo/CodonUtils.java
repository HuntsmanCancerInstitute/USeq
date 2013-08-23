package edu.utah.ames.bioinfo;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

//import dnaModels.DNASequence;

/**
 * A class with some static methods to help with basic codon operations. Also implements some functions related to a 
 * universal genetic code.  
 * @author brendan
 *
 */
public class CodonUtils {

	public static enum Codon {
		AAA, AAG, AAC, AAT, 
		AGA, AGG, AGC, AGT,
		ACA, ACG, ACC, ACT,
		ATA, ATG, ATC, ATT,
		
		GAA, GAG, GAC, GAT, 
		GGA, GGG, GGC, GGT,
		GCA, GCG, GCC, GCT,
		GTA, GTG, GTC, GTT,
		
		CAA, CAG, CAC, CAT, 
		CGA, CGG, CGC, CGT,
		CCA, CCG, CCC, CCT,
		CTA, CTG, CTC, CTT,
		
		TAA, TAG, TAC, TAT, 
		TGA, TGG, TGC, TGT,
		TCA, TCG, TCC, TCT,
		TTA, TTG, TTC, TTT,
	}
	
	public enum AAChar {
		F, S, Y, C, L, W, P, H, R, Q, 
		I, T, N, K, V, A, D, G, E, M, X
	}
	
	public enum AminoAcid {
		Phe, Ser, Tyr, Cys, Leu, Trp, 
		Pro, His, Arg, Ile, Thr, Asn, 
		Glu, Lys, Val, Ala, Asp, Gly,
		Met, Gln, Stop
	}
	
	Map<Codon, AAChar> codonCharMap;
	
	Map<AAChar, AminoAcid> codonAAMap;
	
//	Map<AminoAcid, String> AANameMap;
	
	
	public CodonUtils() {
		createGeneticCode();
		createAAMap();
		//createAANameMap();
	}

	/**
	 * Reverse-lookup the AAChar that corresponds to the given amino acid
	 * @param aa
	 * @return
	 */
	public AAChar lookupChar(AminoAcid aa) {
		for(AAChar c : codonAAMap.keySet()) {
			AminoAcid ca = codonAAMap.get(c);
			if (ca == aa)
				return c;
		}
		return null;
			
	}
	
	public char toChar(AminoAcid aa) {
		AAChar aac = lookupChar(aa);
		switch (aac) {
		case F : return 'F';
		case S : return 'S';
		case Y : return 'Y';
		case C : return 'C';
		case L : return 'L';
		case W : return 'W';
		case P : return 'P';
		case H : return 'H';
		case R : return 'R';
		case Q : return 'Q';
		case I : return 'I';
		case T : return 'T';
		case N : return 'N';
		case K : return 'K';
		case V : return 'V';
		case A : return 'A';
		case D : return 'D';
		case G : return 'G';
		case E : return 'E';
		case M : return 'M';
		case X : return 'X';
		}
		return '?';
	}
	
	public static AminoAcid translate(char first, char second, char third) {

		if (first == 'A') {
			if (second=='A') {
				if (third=='A') return AminoAcid.Lys;
				if (third=='G') return AminoAcid.Lys;
				return AminoAcid.Asp;
			}
			if (second=='C') {
				return AminoAcid.Thr;
			}
			if (second=='T') {
				if (third=='G') return AminoAcid.Met;
				return AminoAcid.Ile;
			}
			//Must be AGX
			if (third=='T') return AminoAcid.Ser;
			if (third=='C') return AminoAcid.Ser;
			return AminoAcid.Arg;
		}

		if (first=='T') {
			if (second=='A') { //TAX
				if (third=='T') return AminoAcid.Tyr;
				if (third=='C') return AminoAcid.Tyr;
				return AminoAcid.Stop;
			}
			if (second=='C') {  //TCX
				return AminoAcid.Ser;
			}
			if (second=='T') { //TTX
				if (third=='T') return AminoAcid.Phe;
				if (third=='C') return AminoAcid.Phe;
				return AminoAcid.Leu;
			}
			//Must be TGX
			if (third=='T') return AminoAcid.Cys;
			if (third=='C') return AminoAcid.Cys;
			if (third=='A') return AminoAcid.Stop;
			if (third=='G') return AminoAcid.Trp; 
		}
		
		if (first=='C') {
			if (second=='A') {
				if (third=='T') return AminoAcid.His;
				if (third=='C') return AminoAcid.His;
				return AminoAcid.Gln;
			}
			if (second=='T') {
				return AminoAcid.Leu;
			}
			if (second=='C') {
				return AminoAcid.Pro;
			}
			//Must be CGX
			return AminoAcid.Arg;
			
		}
		
		if (first=='G') {
			if (second=='A') {
				if (third=='A') return AminoAcid.Glu;
				if (third=='G') return AminoAcid.Glu;
				return AminoAcid.Asp;
			}
			if (second=='T') {
				return AminoAcid.Val;
			}
			if (second=='C') {
				return AminoAcid.Ala;
			}
			return AminoAcid.Gly;
		}
		
		throw new IllegalArgumentException("Could not identify AA for codon : " + first + "" + second + "" + third);
		
	}
	/**
	public static AminoAcid translate(DNASequence seq, int startSite) {
		char first = seq.getBaseChar(startSite);
		char second = seq.getBaseChar(startSite+1);
		char third = seq.getBaseChar(startSite+2);
		return translate(first, second, third);
	}
	*/
	
	public static AminoAcid translate(String codon) {
		if (codon.startsWith("A")) {
			if (codon.charAt(1)=='A') {
				if (codon.charAt(2)=='A') return AminoAcid.Lys;
				if (codon.charAt(2)=='G') return AminoAcid.Lys;
				return AminoAcid.Asp;
			}
			if (codon.charAt(1)=='C') {
				return AminoAcid.Thr;
			}
			if (codon.charAt(1)=='T') {
				if (codon.charAt(2)=='G') return AminoAcid.Met;
				return AminoAcid.Ile;
			}
			//Must be AGX
			if (codon.charAt(2)=='T') return AminoAcid.Ser;
			if (codon.charAt(2)=='C') return AminoAcid.Ser;
			return AminoAcid.Arg;
		}
		
		if (codon.startsWith("T")) {
			if (codon.charAt(1)=='A') { //TAX
				if (codon.charAt(2)=='T') return AminoAcid.Tyr;
				if (codon.charAt(2)=='C') return AminoAcid.Tyr;
				return AminoAcid.Stop;
			}
			if (codon.charAt(1)=='C') {  //TCX
				return AminoAcid.Ser;
			}
			if (codon.charAt(1)=='T') { //TTX
				if (codon.charAt(2)=='T') return AminoAcid.Phe;
				if (codon.charAt(2)=='C') return AminoAcid.Phe;
				return AminoAcid.Leu;
			}
			//Must be TGX
			if (codon.charAt(2)=='T') return AminoAcid.Cys;
			if (codon.charAt(2)=='C') return AminoAcid.Cys;
			if (codon.charAt(2)=='A') return AminoAcid.Stop;
			if (codon.charAt(2)=='G') return AminoAcid.Trp; 
		}
		
		if (codon.startsWith("C")) {
			if (codon.charAt(1)=='A') {
				if (codon.charAt(2)=='T') return AminoAcid.His;
				if (codon.charAt(2)=='C') return AminoAcid.His;
				return AminoAcid.Gln;
			}
			if (codon.charAt(1)=='T') {
				return AminoAcid.Leu;
			}
			if (codon.charAt(1)=='C') {
				return AminoAcid.Pro;
			}
			//Must be CGX
			return AminoAcid.Arg;
			
		}
		
		if (codon.startsWith("G")) {
			if (codon.charAt(1)=='A') {
				if (codon.charAt(2)=='A') return AminoAcid.Glu;
				if (codon.charAt(2)=='G') return AminoAcid.Glu;
				return AminoAcid.Asp;
			}
			if (codon.charAt(1)=='T') {
				return AminoAcid.Val;
			}
			if (codon.charAt(1)=='C') {
				return AminoAcid.Ala;
			}
			return AminoAcid.Gly;
		}
		
		throw new IllegalArgumentException("Could not identify AA for : " + codon);
	}
	
	private Map<AminoAcid, String> createAANameMap() {
		Map<AminoAcid, String> names = new HashMap<AminoAcid, String>();
		
		names.put(AminoAcid.Ala, "Alanaine");
		names.put(AminoAcid.Arg, "Arginine");
		names.put(AminoAcid.Asn, "Asparagine");
		names.put(AminoAcid.Asp, "Aspartic acid");
		names.put(AminoAcid.Cys, "Cysteine");
		names.put(AminoAcid.Glu, "Glutamic acid");
		names.put(AminoAcid.Gln, "Glutamine");
		names.put(AminoAcid.Gly, "Glycine");
		names.put(AminoAcid.His, "Histidine");
		names.put(AminoAcid.Ile, "Isoleucine");
		names.put(AminoAcid.Leu, "Leucine");
		names.put(AminoAcid.Lys, "Lysine");
		names.put(AminoAcid.Met, "Methionine");
		names.put(AminoAcid.Phe, "Phenylalanine");
		names.put(AminoAcid.Pro, "Proline");
		names.put(AminoAcid.Ser, "Serine");
		names.put(AminoAcid.Thr, "Threonine");
		names.put(AminoAcid.Trp, "Tryptophan");
		names.put(AminoAcid.Tyr, "Tyrosine");
		names.put(AminoAcid.Val, "Valine");
		return names;
	}

	private void createAAMap() {
		codonAAMap = new HashMap<AAChar, AminoAcid>();
		
		codonAAMap.put(AAChar.F, AminoAcid.Phe);
		codonAAMap.put(AAChar.S, AminoAcid.Ser);
		codonAAMap.put(AAChar.Y, AminoAcid.Tyr);
		codonAAMap.put(AAChar.C, AminoAcid.Cys);
		codonAAMap.put(AAChar.L, AminoAcid.Leu);
		codonAAMap.put(AAChar.W, AminoAcid.Trp);
		codonAAMap.put(AAChar.P, AminoAcid.Pro);
		codonAAMap.put(AAChar.H, AminoAcid.His);
		codonAAMap.put(AAChar.R, AminoAcid.Arg);
		codonAAMap.put(AAChar.Q, AminoAcid.Gln);
		codonAAMap.put(AAChar.I, AminoAcid.Ile);
		codonAAMap.put(AAChar.T, AminoAcid.Thr);
		codonAAMap.put(AAChar.N, AminoAcid.Asn);
		codonAAMap.put(AAChar.M, AminoAcid.Met);
		codonAAMap.put(AAChar.K, AminoAcid.Lys);
		codonAAMap.put(AAChar.V, AminoAcid.Val);
		codonAAMap.put(AAChar.A, AminoAcid.Ala);
		codonAAMap.put(AAChar.D, AminoAcid.Asp);
		codonAAMap.put(AAChar.G, AminoAcid.Gly);
		codonAAMap.put(AAChar.E, AminoAcid.Glu);
		
		codonAAMap.put(AAChar.X, AminoAcid.Stop);
	}
	
	private void createGeneticCode() {
		codonCharMap = new HashMap<Codon, AAChar>();
			
		codonCharMap.put(Codon.TTT, AAChar.F);
		codonCharMap.put(Codon.TTC, AAChar.F);

		codonCharMap.put(Codon.TAT, AAChar.Y);
		codonCharMap.put(Codon.TAC, AAChar.Y);
		codonCharMap.put(Codon.TGT, AAChar.C);
		codonCharMap.put(Codon.TGC, AAChar.C);

		codonCharMap.put(Codon.TCT, AAChar.S);
		codonCharMap.put(Codon.TCC, AAChar.S);
		codonCharMap.put(Codon.TCA, AAChar.S);
		codonCharMap.put(Codon.TCG, AAChar.S);
		
		codonCharMap.put(Codon.TGG, AAChar.W);

		codonCharMap.put(Codon.TTA, AAChar.L);
		codonCharMap.put(Codon.TTG, AAChar.L);
		codonCharMap.put(Codon.CTT, AAChar.L);
		codonCharMap.put(Codon.CTA, AAChar.L);
		codonCharMap.put(Codon.CTC, AAChar.L);
		codonCharMap.put(Codon.CTG, AAChar.L);
		
		codonCharMap.put(Codon.CCT, AAChar.P);
		codonCharMap.put(Codon.CCC, AAChar.P);
		codonCharMap.put(Codon.CCG, AAChar.P);
		codonCharMap.put(Codon.CCA, AAChar.P);
		
		codonCharMap.put(Codon.CAT, AAChar.H);
		codonCharMap.put(Codon.CAC, AAChar.H);
		
		codonCharMap.put(Codon.CGT, AAChar.R);
		codonCharMap.put(Codon.CGC, AAChar.R);
		codonCharMap.put(Codon.CGA, AAChar.R);
		codonCharMap.put(Codon.CGG, AAChar.R);
		codonCharMap.put(Codon.AGA, AAChar.R);
		codonCharMap.put(Codon.AGG, AAChar.R);

		codonCharMap.put(Codon.ATT, AAChar.I);
		codonCharMap.put(Codon.ATC, AAChar.I);
		codonCharMap.put(Codon.ATA, AAChar.I);
		
		codonCharMap.put(Codon.ATG, AAChar.M);
		
		codonCharMap.put(Codon.ACT, AAChar.T);
		codonCharMap.put(Codon.ACG, AAChar.T);
		codonCharMap.put(Codon.ACA, AAChar.T);
		codonCharMap.put(Codon.ACC, AAChar.T);
		
		codonCharMap.put(Codon.CAA, AAChar.Q);
		codonCharMap.put(Codon.CAG, AAChar.Q);
		
		codonCharMap.put(Codon.AAT, AAChar.N);
		codonCharMap.put(Codon.AAC, AAChar.N);
		
		codonCharMap.put(Codon.AAA, AAChar.K);
		codonCharMap.put(Codon.AAG, AAChar.K);
		
		codonCharMap.put(Codon.GTT, AAChar.V);
		codonCharMap.put(Codon.GTA, AAChar.V);
		codonCharMap.put(Codon.GTG, AAChar.V);
		codonCharMap.put(Codon.GTC, AAChar.V);
		
		codonCharMap.put(Codon.GCT, AAChar.A);
		codonCharMap.put(Codon.GCC, AAChar.A);
		codonCharMap.put(Codon.GCA, AAChar.A);
		codonCharMap.put(Codon.GCG, AAChar.A);
		
		codonCharMap.put(Codon.GAT, AAChar.D);
		codonCharMap.put(Codon.GAC, AAChar.D);
		codonCharMap.put(Codon.GAA, AAChar.E);
		codonCharMap.put(Codon.GAG, AAChar.E);
		
		codonCharMap.put(Codon.GGT, AAChar.G);
		codonCharMap.put(Codon.GGC, AAChar.G);
		codonCharMap.put(Codon.GGA, AAChar.G);
		codonCharMap.put(Codon.GGG, AAChar.G);
		
		codonCharMap.put(Codon.TAG, AAChar.X);
		codonCharMap.put(Codon.TAA, AAChar.X);
		codonCharMap.put(Codon.TGA, AAChar.X);
	}
	
}