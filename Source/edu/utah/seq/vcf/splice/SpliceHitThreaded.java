package edu.utah.seq.vcf.splice;

import java.util.ArrayList;

import util.bio.parsers.UCSCGeneLine;
import util.gen.Num;
import edu.utah.seq.vcf.VCFRecord;
import edu.utah.seq.vcf.xml.foundation.SimpleVcf;

public class SpliceHitThreaded {
	//fields
	private SimpleVcf vcf;
	private int vcfAltIndex;
	private UCSCGeneLine transcript;
	private ArrayList<SpliceJunction> affectedSpliceJunctions = null;
	
	//constructors
	public SpliceHitThreaded (SimpleVcf vcf, int vcfAltIndex, UCSCGeneLine transcript){
		this.vcf = vcf;
		this.vcfAltIndex = vcfAltIndex;
		this.transcript = transcript;
	}
	
	public void saveSpliceJunction(SpliceJunction sj){
		if (affectedSpliceJunctions == null) affectedSpliceJunctions = new ArrayList<SpliceJunction>();
		affectedSpliceJunctions.add(sj);
	}
	
	/**Returns SJType,SJPositon,GeneName,VCFAltSeq,SJ-10Log10(pval),SJRefScore,SJAltScore for every affected splice junctions that meet the categories
	 * 0 all, 1 damaging, 2 damaging and novel exonic and splice.*/
	public ArrayList<String> getVcfEntries(int category){
		String geneName = transcript.getDisplayNameThenName();
		ArrayList<String> anno = new ArrayList<String>();
		for (int i=0; i< affectedSpliceJunctions.size(); i++){			
			SpliceJunction sj = affectedSpliceJunctions.get(i);
			String sjType = sj.getType();
			//save it? 0 all, 1 just splice damage, 2 no intronic
			if (category == 1){
				if (sjType.startsWith("D") == false) continue;
			}
			else if (category == 2) {
				if (sjType.contains("I")) continue;
			}
			
			//TranscriptName,Type,Position,RefSeq,AltSeq,RefScore,AltScore,ScoreDelta
			int sjPos = sj.getPosition() +1; //convert back to 1 base
			String sjRefScore = Num.formatNumber(sj.getReferenceScore(), 2);
			String sjAltScore = Num.formatNumber(sj.getAlternateScore(), 2);
			String refSeq = sj.getReferenceSequence();
			String altSeq = sj.getAlternateSequence();
			String scoreDelta = Num.formatNumber(sj.getScoreDelta(), 2);
			anno.add(geneName+":"+sjType+","+sjPos+","+refSeq+","+altSeq+","+sjRefScore+","+sjAltScore+","+scoreDelta);
		}
		return anno;
	}
	
	//getters and setters
	public SimpleVcf getVcf() {
		return vcf;
	}
	public void setVcf(SimpleVcf vcf) {
		this.vcf = vcf;
	}
	public int getVcfAltIndex() {
		return vcfAltIndex;
	}
	public void setVcfAltIndex(int vcfAltIndex) {
		this.vcfAltIndex = vcfAltIndex;
	}
	public UCSCGeneLine getTranscript() {
		return transcript;
	}
	public void setTranscript(UCSCGeneLine transcript) {
		this.transcript = transcript;
	}

	public ArrayList<SpliceJunction> getAffectedSpliceJunctions() {
		return affectedSpliceJunctions;
	}
	
}
