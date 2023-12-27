package edu.utah.seq.run.invitae;

import java.io.IOException;
import htsjdk.samtools.reference.ReferenceSequence;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;

public class InvitaeVariant implements Comparable<InvitaeVariant> {
	
	private String varChromosome = null;
	private int varPosHg19 = 0;
	private String varRef = null;
	private String varAlt = null;
	private String varGene = null;
	private String varTranscript = null;
	private String varCDot = null;
	private String varPDot = null;
	private String varZygosity = null;
	private String varInterpretation = null;
	private String varMethod = null;
	private String geneEndHg38 = null;
	
	private boolean passed = false;


	public InvitaeVariant(String[] info, InvitaeBulkCsvParser idw) throws IOException {
		// Interp: Uncertain Significance; Pathogenic; Likely Pathogenic; Likely Benign; Pathogenic (low penetrance); Benign; Increased Risk Allele; Benign (reportable variant)
		String x = info[idw.varInterpretation].trim();
		if (x.length()!=0) varInterpretation = Misc.WHITESPACE.matcher(x).replaceAll("_");
		
		// Gene: ABCA4; ACADM; ACD; ACSF3; ACTA2; ADA; AIP; ALDH18A1; ALK; ALMS1; ANKRD26; AP3B1; AP3D1; APC; APOA4; ATM; ATP7B; AXIN2; B4GALT7; BAP1; BARD1; BLM; BMPR1A; BRCA1; BRCA2; BRIP1; BUB1B; C17ORF62; C2; C9; CASP10; CASP8; CASR; CCM2; CDC73; CDH1; CDKN1B; CDKN1C; CDKN2A; CDKN2A (P14ARF); CDKN2A (P16INK4A); CEBPA; CFTR; CHEK2; CLPB; CSF3R; CTC1; CTLA4; CTNNA1; CTRC; CTSC; CYBA; CYP21A2; DDX41; DIAPH1; DICER1; DIS3L2; DKC1; DNMT3B; DOCK2; DOCK8; DUOX2; EGFR; ELANE; ELP1; EPCAM; ERCC4; ERCC6L2; ETV6; EXT1; EYS; F2; F5; FANCA; FANCB; FANCC; FANCD2; FANCF; FANCI; FAS; FH; FKRP; FLCN; FLNC; GAA; GATA2; GBE1; GFI1; GJB2; GPC3; GREM1; HFE; HJV; HOXB13; HPS6; HRAS; IFIH1; IKZF1; IL10RA; IL17RA; IL2RB; ITGAM; JAK3; KCNA5; KCNJ2; KIT; KLHDC8B; LDLR; LOXHD1; LRRK2; LYST; LZTR1; MALT1; MAP2K1; MAX; MBD4; MC1R; MECOM; MEFV; MEN1; MET; MITF; MLH1; MLH3; MRE11A; MSH2; MSH3; MSH6; MUTYH; MYBPC3; MYH7; MYH9; NAF1; NBN; NF1; NF2; NFKB2; NOD2; NOP10; NOTCH1; NTHL1; PAH; PALB2; PALLD; PARN; PDGFRA; PHOX2B; PIK3R1; PMS2; PMS2 OR PMS2CL; POLD1; POLE; POT1; PRKAR1A; PRKDC; PROC; PROS1; PRSS1; PSMB4; PTCH1; PTEN; PTPN11; RAD50; RAD51C; RAD51D; RASGRP1; RASGRP2; RB1; RBM20; RECQL4; REST; RET; RFX5; RMRP; RPL11; RPL5; RPS19; RPS20; RTEL1; RUNX1; SAMD9; SAMD9L; SDHA; SDHB; SDHC; SDHD; SGPL1; SH2D1A; SLX4; SMAD4; SMARCA4; SMARCAL1; SMARCB1; SMARCE1; SPG11; SPINK1; STAT4; STIM1; STK11; STN1; STX11; STXBP2; SUFU; TAOK2; TCIRG1; TERC; TERT; TET2; TICAM1; TIMM50; TMC6; TMEM127; TNNT2; TP53; TSC1; TSC2; TUBB1; USB1; VHL; VPS13A; VPS13B; WAS; WRN; WT1; XPA; ZBTB24; 
		// Odd:  CDKN2A (P16INK4A); PMS2 OR PMS2CL
		x = info[idw.varGene].trim();
		if (x.length()!=0) varGene = Misc.WHITESPACE.matcher(x).replaceAll("_");
		else throw new IOException ("ERROR: no gene name in "+Misc.stringArrayToString(info, " "));
		//check them
		
	
		// Method
		x = info[idw.varMethod].trim();
		if (x.length()!=0) varMethod = Misc.WHITESPACE.matcher(x).replaceAll("_");
		
		// Transcript
		x = info[idw.varTranscript].trim();
		if (x.length()!=0) varTranscript = x;
		
		// cDot
		x = info[idw.varCDot].trim();
		if (x.length()!=0) varCDot = Misc.WHITESPACE.matcher(x).replaceAll("_");
		
		// pDot
		x = info[idw.varPDot].trim();
		if (x.length()!=0) varPDot = Misc.WHITESPACE.matcher(x).replaceAll("_");
		
		// zygosity
		x = info[idw.varZygosity].trim();
		if (x.length()!=0) varZygosity = Misc.WHITESPACE.matcher(x).replaceAll("_");
		
		// chromosome, 
		x = info[idw.varChromosome].trim();
		if (x.length()!=0) varChromosome = "chr"+x;
		
		// vcf pos
		x = info[idw.varPos].trim();
		if (x.length()!=0) varPosHg19 = Integer.parseInt(x);
		
		// vcf ref, can contain .
		x = info[idw.varRef].trim();
		if (x.length()!=0) varRef = x;
		
		// vcf alt, can contain <DEL> and <DUP>
		x = info[idw.varAlt].trim();
		if (x.length()!=0) varAlt = x;
		
		// missing pos ref or alt?
		if (varPosHg19 == 0 || varRef == null || varAlt == null) {
			IO.el("\t\tWARNING: failed to find a vcf_pos, vcf_ref, or vcf_alt, skipping -> "+Misc.stringArrayToString(info, " "));
			return;
		}
		varRef = varRef.toUpperCase();
		
		//missing chrom? or pos? <DEL> or <DUP>
		if (varChromosome == null || varAlt.contains("<")) {
			// InvitaeGeneName,Chrom,Start,Stop,RefseqName,Strand
			//       0           1     2     3      4        5
			Bed iGeneCoor = idw.getiGenCoor().get(varGene);
			if (iGeneCoor == null) throw new IOException ("ERROR: failed to find gene info for: "+Misc.stringArrayToString(info, " ")+"\nPlease add coordinate info for this in the InvitaeDataWrangler.jar ");
			if (varChromosome == null ) varChromosome = iGeneCoor.getChromosome();
			//set hg38 end position, not hg19 since crossmap doesn't convert these
			if (varAlt.contains("<")) geneEndHg38 = new Integer(iGeneCoor.getStop()).toString();
		}
		
		// check the sequence, this also checks for hg19
		ReferenceSequence p = idw.getFasta().getSubsequenceAt(varChromosome, varPosHg19, varPosHg19+varRef.length()-1);
		String ref = new String(p.getBases()).toUpperCase();
		//missing ref base?
		if (varRef.equals(".")) varRef = ref;
		else if (varRef.equals(ref) == false) throw new IOException ("ERROR: vcf_ref '"+varRef+"' in '"+Misc.stringArrayToString(info," ")+"' does not match hg19 ref base '"+ref+ "' for position "+varPosHg19);
		
		//made it this far so it passes
		passed = true;

	}
	
	public int compareTo(InvitaeVariant other) {
		//sort by chromosome
		int x = varChromosome.compareTo(other.varChromosome);
		if (x !=0) return x;
		//sort by position
		if (varPosHg19 < other.varPosHg19) return -1;
		if (varPosHg19 > other.varPosHg19) return 1;
		return 0;
	}
	
	/** Returns CHROM POS ID REF ALT QUAL FILTER INFO, no sample info though */
	public String toVcf(int id){
		StringBuilder sb = new StringBuilder();
		//CHROM
		sb.append(varChromosome); sb.append("\t");
		//POS
		sb.append(varPosHg19); sb.append("\tInvitae_");
		//ID
		sb.append(id);
		sb.append("\t");		
		//REF
		sb.append(varRef); sb.append("\t");
		//ALT
		sb.append(varAlt); sb.append("\t");
		//QUAL FILTER
		sb.append(".\t.\t");
		//INFO
		sb.append("iGene="); sb.append(varGene);
		sb.append(";iInterp="); sb.append(varInterpretation);
		sb.append(";iTrans="); sb.append(varTranscript);
		sb.append(";iCDot="); sb.append(varCDot);
		sb.append(";iPDot="); sb.append(varPDot);
		sb.append(";iZyg="); sb.append(varZygosity);
		sb.append(";iMeth="); sb.append(varMethod);
		
		if (geneEndHg38!=null) {
			sb.append(";Hg38End="); sb.append(geneEndHg38);
		}
		
		return sb.toString();
		
	}

	public boolean isPassed() {
		return passed;
	} 

	

}
