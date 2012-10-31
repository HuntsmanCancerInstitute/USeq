package edu.expr;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import edu.utah.seq.useq.data.Region;

import util.bio.annotation.*;
import util.gen.*;
import util.bio.parsers.*;
import util.bio.seq.Seq;

public class Allele extends Region {

	//fields
	/**Sequence with which to replace start-stop region.*/
	private String allele;
	/**Many be null*/
	private String notes;
	
	//for NovoalignIndelParser
	private StringBuilder posteriorProbabilites;
	private StringBuilder baseQualities;
	private StringBuilder readStartPositions;
	public static final Pattern comma = Pattern.compile(",");
	
	//for Alleler
	public static final HashMap<String,String> codon2AA = Seq.fetchCodonAAHashMap();
	public static final Pattern gatc = Pattern.compile("[GATC]");
	public static final HashMap<String,String> dnaAmbigBase2BasesHashMap = Seq.fetchDNAAmbigBase2BasesHashMap();
	

	//constructor
	public Allele (int start, int stop, String allele, String notes){
		super (start, stop);
		this.allele = allele;
		this.notes = notes;
	}

	public String toString(){
		return this.start+"\t"+this.stop+"\t"+allele+"\t"+notes;
	}
	
	/**Returns start stop allele posteriorProbabilities baseScores totalNumberReads numberUniqueReads
	 * provided the minimumNumberUniqueReads is met. Otherwise null.*/
	public String alleleTableLine(int minimumNumberUniqueReads){
		//see if meets minimum
		String[] strandStart = comma.split(readStartPositions);
		int numUniqueReads = Misc.loadHashSet(strandStart).size();
		if (numUniqueReads < minimumNumberUniqueReads) return null;
		//trim off final comma
		String pp = posteriorProbabilites.toString().substring(0,posteriorProbabilites.length()-1);
		String bs = baseQualities.toString().substring(0,baseQualities.length()-1);
		return start +"\t"+ stop +"\t"+ allele +"\t"+ pp +"\t"+ bs +"\t"+ strandStart.length +"\t"+ numUniqueReads; 
	}
	
	/**Start-1\tStop+1\tbasesInserted\t#UniqueOverlappingAlignments\tBlank*/
	public String fetchInsertionBedLine(){
		if (start != stop) return null;
		String[] strandStart = comma.split(readStartPositions);
		int numUniqueReads = Misc.loadHashSet(strandStart).size();
		return (start-1) + "\t" + (stop+1)+"\t"+allele+"\t"+numUniqueReads+"\t.";
	}

	/**Start\tStop\tBlank\t#UniqueOverlappingAlignments\tBlank*/
	public String fetchDeletionBedLine(){
		if (start == stop) return null;
		String[] strandStart = comma.split(readStartPositions);
		int numUniqueReads = Misc.loadHashSet(strandStart).size();
		return start + "\t" + stop +"\t.\t"+"\t"+numUniqueReads+"\t.";
	}
	
	public String affect(UCSCGeneLine geneModel, String chromSeq){
		String effects;
		int txStart = geneModel.getTxStart();
		int txEnd = geneModel.getTxEnd();
		boolean sense = geneModel.getStrand().equals("+");
		//Entirely 5' of gene?
		int diff = txStart-stop;
		if (diff>=0) {
			if (sense) effects = "\t\t"+diff +"bp 5' of gene\n";
			else effects = "\t\t"+diff +"bp 3' of gene\n";
			return effects;
		}
		//Entirely 3' of gene?
		diff = start - txEnd;
		if (diff>=0) {
			if (sense) effects = "\t\t"+diff +"bp 3' of gene\n";
			else effects = "\t\t"+diff +"bp 5' of gene\n";
			return effects;
		}
		StringBuilder sb = new StringBuilder();
		//Affects 5'
		int st = txStart -10000;
		if (st <0) st = 0;
		if (this.intersects(st, txStart)){
			if (sense) sb.append("\t\t5' of gene\n");
			else sb.append("\t\t3' of gene\n");
		}
		//Affects 3'
		if (this.intersects(txEnd, txEnd + 10000)){
			if (sense) sb.append("\t\t3' of gene\n");
			else sb.append("\t\t5' of gene\n");
		}

		int cdsStart = geneModel.getCdsStart();
		int cdsEnd = geneModel.getCdsEnd();
		boolean cds = true;

		//5' UTR
		diff = cdsStart -stop;
		if (diff>=0) {
			if (sense) sb.append("\t\t5' UTR\n");
			else sb.append("\t\t3' UTR\n");
			cds = false;
		}
		//3' UTR
		diff = start - cdsEnd;		
		if (diff>=0) {
			if (sense) sb.append("\t\t3' UTR\n");
			else sb.append("\t\t5' UTR\n");
			cds = false;
		}
		//exonic
		boolean aaHit = false;
		ExonIntron[] exons = geneModel.getExons();
		for (int i=0; i< exons.length; i++){
			if (intersects(exons[i].getStart(), exons[i].getEnd())){
				//hit coding?
				if (cds) {
					sb.append("\t\tCoding exon ");
					aaHit = true;
				}
				else sb.append("\t\tNon-coding exon ");
				if (sense) sb.append((i+1)+"");
				else sb.append( (exons.length-i) +"");
				//hit splice junction? 1st base or last 3
				if (sense){
					//1st base or last 3
					if ((intersects(exons[i].getStart(), exons[i].getStart()+1) && i!=0) || (intersects(exons[i].getEnd()-3, exons[i].getEnd()) && i!=(exons.length-1))) sb.append(" splice site");
				}
				else {
					//1st 3 bases or last 1
					if ((intersects(exons[i].getStart(), exons[i].getStart()+3) && i!=0) || (intersects(exons[i].getEnd()-1, exons[i].getEnd()) && i!=(exons.length-1))) sb.append(" splice site");
				}
				sb.append("\n");
			}
		}
		//calc aa change
		if (aaHit){
			//Find real seq under allele
			String base = chromSeq.substring(start, stop);
			String altAllele = allele.toUpperCase();
			//attempt to change ambigous base symbol to the opposite of the reference base
			if (allele.length() == 1){
				Matcher mat = gatc.matcher(altAllele);
				if (mat.matches() == false) {
					String x = Seq.findOppositeBase(base, altAllele, dnaAmbigBase2BasesHashMap);
					if (x == null) Misc.printExit("\nError: cannot find bases associated with ambiguous base -> "+allele);
					if (x.length()==1) altAllele = x;
				}
			}
			//collect the transcript sequence for exons
			StringBuilder mutSeq = new StringBuilder();
			StringBuilder refSeq = new StringBuilder();
			for (int i=0; i< exons.length; i++){
				//does it contain the cdsStart?
				if (exons[i].contains(cdsStart)){
					int endT;
					if (exons[i].contains(cdsEnd)) endT = cdsEnd;
					else endT = exons[i].getEnd();
					String seq = new String( chromSeq.substring(cdsStart, endT));
					refSeq.append(seq);
					//does allele effect this seq
					if (this.intersects(cdsStart, endT)){
						int startAllele = start - cdsStart;
						if (startAllele < 0) startAllele = 0;
						int endAllele = stop - cdsStart;
						if (endAllele> seq.length()) endAllele = seq.length();
						//make splice product
						String mutL = seq.substring(0,startAllele);
						String mutR = seq.substring(endAllele);
						mutSeq.append(mutL);
						mutSeq.append(altAllele);
						mutSeq.append(mutR);
						//System.out.println("CdsStart:\n\tmutL "+mutL+"\n\talle "+altAllele+"\n\tmutR "+mutR);
					}
					else mutSeq.append(seq);
					
				}
				//does it contain the cdsEnd?
				else if (exons[i].contains(cdsEnd)){
					String seq = new String( chromSeq.substring(exons[i].getStart(), cdsEnd));
					refSeq.append(seq);
					//does allele effect this seq
					if (this.intersects(exons[i].getStart(), cdsEnd)){
						int startAllele = start - exons[i].getStart();
						if (startAllele < 0) startAllele = 0;
						int endAllele = stop - exons[i].getStart();
						if (endAllele> seq.length()) endAllele = seq.length();
						//make splice product
						String mutL = seq.substring(0,startAllele);
						String mutR = seq.substring(endAllele);
						mutSeq.append(mutL);
						mutSeq.append(altAllele);
						mutSeq.append(mutR);
						//System.out.println("CdsEnd:\n\tmutL "+mutL+"\n\talle "+altAllele+"\n\tmutR "+mutR);
					}
					else mutSeq.append(seq);
				}
				//is it internal to cds's?
				else if (exons[i].getEnd()<= cdsEnd && exons[i].getStart()>= cdsStart) {
					String seq = exons[i].getSequence(chromSeq);
					refSeq.append(seq);
					//does allele effect this seq
					if (this.intersects(exons[i].getStart(), exons[i].getEnd())){
						int startAllele = start - exons[i].getStart();
						if (startAllele < 0) startAllele = 0;
						int endAllele = stop - exons[i].getStart();
						if (endAllele> seq.length()) endAllele = seq.length();
						//make splice product
						String mutL = seq.substring(0,startAllele);
						String mutR = seq.substring(endAllele);
						mutSeq.append(mutL);
						mutSeq.append(altAllele);
						mutSeq.append(mutR);
						//System.out.println("Internal:\n\tmutL "+mutL+"\n\talle "+altAllele+"\n\tmutR "+mutR);
					}
					else mutSeq.append(seq);
				}
			}
			//make sequences
			String refSeqString;
			String mutSeqString;
			if (sense == false){
				refSeqString = Seq.reverseComplementDNA(refSeq.toString());
				mutSeqString = Seq.reverseComplementDNA(mutSeq.toString());
			}
			else {
				refSeqString = refSeq.toString();
				mutSeqString = mutSeq.toString();
			}
			//compare seqs
			String aaChanges = compareSequences(refSeqString, mutSeqString);
			
//System.out.println("\nRef: "+ refSeqString);
//System.out.println("\nMut: "+ mutSeqString);
			
			sb.append("\t\t\t"+aaChanges+"\n");
		}
		
		//intronic
		ExonIntron[] introns = geneModel.getIntrons();
		if (introns != null){
			for (int i=0; i< introns.length; i++){			
				if (intersects(introns[i].getStart(), introns[i].getEnd())){
					sb.append("\t\tIntron ");
					if (sense) sb.append((i+1));
					else sb.append( (introns.length-i));
					//splice-site?
					if (intersects(introns[i].getStart(), introns[i].getStart()+2) || intersects(introns[i].getEnd()-2, introns[i].getEnd())){
						sb.append(" splice site");
					}
					sb.append("\n");
				}
				//look for insertion at junctions
				else if (start == stop){
					if (start == (introns[i].getStart()) || stop == introns[i].getEnd()){
						sb.append("\t\tIntron ");
						if (sense) sb.append((i+1));
						else sb.append( (introns.length-i));
						sb.append(" splice site\n");
					}
				}
			}
		}
		return sb.toString();
	}
	
	/**Takes two DNA sequences, translates them and finds the difference.*/
	public String compareSequences(String refSeqDNA, String mutSeqDNA){
		//translate
		String refAA = Seq.translate(refSeqDNA.toUpperCase(), codon2AA);
		String mutAA = Seq.translate(mutSeqDNA.toUpperCase(), codon2AA);
		
		//System.out.println("RefAA: "+refAA);
		//System.out.println("MutAA: "+mutAA);
		
		String[] nonCommon = Misc.trimCommon(new String[]{refAA, mutAA});
		if (nonCommon[0].equals(nonCommon[1])) return "Synonymous";
		return "Non-synonymous "+nonCommon[0]+" > "+nonCommon[1];
	}


	/**Parses a tab delimited file (chr, start, stop, allele, notes...), zip/ gz OK.
	 * @param alleleFile, skips empty lines and those starting with '#'
	 * @param subStart and subEnd are the number to subtract from the ends of each region
	 * @return a HashMap<Chr,sorted Region[]> or null in none are found
	 * */
	public static HashMap<String,Allele[]> parseAlleles (File alleleFile){
		HashMap<String,ArrayList<Allele>> ss = new HashMap<String,ArrayList<Allele>>();
		try{
			BufferedReader in = IO.fetchBufferedReader(alleleFile);
			String line;
			String[] tokens;
			ArrayList<Allele> al = null;
			Pattern tab = Pattern.compile("\\t");
			//chrom, start, stop
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() ==0 || line.startsWith("#")) continue;
				tokens = tab.split(line);
				//does chrom already exist?
				if (ss.containsKey(tokens[0])) al = ss.get(tokens[0]);
				else {
					al = new ArrayList<Allele>();
					ss.put(tokens[0], al);
				}
				int start = Integer.parseInt(tokens[1]);
				int stop = Integer.parseInt(tokens[2]);
				if (start > stop) throw new Exception("\nFound a start that is greater than stop!  Cannot parse file "+alleleFile+", bad line-> "+line);
				String allele = "";
				if (tokens.length > 3) allele = tokens[3];
				String notes = null;
				if (tokens.length > 4){
					StringBuilder sb = new StringBuilder(tokens[4]);
					for (int i=5; i< tokens.length; i++){
						sb.append("\t");
						sb.append(tokens[i]);
					}
					notes = sb.toString();
				}
				al.add(new Allele(start, stop, allele, notes));
			}
		}catch (Exception e){
			e.printStackTrace();
		}
		if (ss.size() == 0) return null;
		//make hashmap
		HashMap<String,Allele[]> ssReal = new HashMap<String,Allele[]>();
		Iterator<String> it = ss.keySet().iterator();
		while (it.hasNext()){
			String chrom = it.next();
			ArrayList<Allele> al = ss.get(chrom);
			Allele[] array = new Allele[al.size()];
			al.toArray(array);
			Arrays.sort(array);
			ssReal.put(chrom, array);
		}
		return ssReal;
	}

	public StringBuilder getPosteriorProbabilites() {
		return posteriorProbabilites;
	}

	public void setPosteriorProbabilites(StringBuilder posteriorProbabilites) {
		this.posteriorProbabilites = posteriorProbabilites;
	}

	public StringBuilder getBaseQualities() {
		return baseQualities;
	}

	public void setBaseQualities(StringBuilder baseQualities) {
		this.baseQualities = baseQualities;
	}

	public StringBuilder getReadStartPositions() {
		return readStartPositions;
	}

	public void setReadStartPositions(StringBuilder readStartPositions) {
		this.readStartPositions = readStartPositions;
	}

	public String getNotes() {
		return notes;
	}

	public void setNotes(String notes) {
		this.notes = notes;
	}
}
