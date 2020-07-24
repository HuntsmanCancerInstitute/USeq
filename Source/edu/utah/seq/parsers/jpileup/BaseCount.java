package edu.utah.seq.parsers.jpileup;

import java.io.IOException;
import java.util.HashSet;

import util.gen.Misc;

public class BaseCount {
	int bpPosition;
	int g = 0;
	int a = 0;
	int t = 0;
	int c = 0;
	int n = 0;
	int ins = 0;
	int del = 0;
	int failQual = 0;
	char ref;
	HashSet<String> readNames = new HashSet<String>();

	public BaseCount(int bpPosition, char ref){
		this.bpPosition = bpPosition;
		this.ref = ref;
	}
	
	public BaseCount(int bpPosition, char ref, String counts){
		this.bpPosition = bpPosition;
		this.ref = ref;
		String[] f = Misc.COMMA.split(counts);
		a= Integer.parseInt(f[0]);
		c= Integer.parseInt(f[1]);
		g= Integer.parseInt(f[2]);
		t= Integer.parseInt(f[3]);
		n= Integer.parseInt(f[4]);
		del= Integer.parseInt(f[5]);
		ins= Integer.parseInt(f[6]);
		failQual= Integer.parseInt(f[7]);
	}

	public double getPassingReadCoverage() {
		return g+a+t+c+ins+del;
	}
	public double getPassingReadCoverageSnv() {
		return g+a+t+c;
	}
	public double getTotalReadCoverage() {
		return g+a+t+c+ins+del+n+failQual;
	}
	public double getIndelCount() {
		return ins+del;
	}
	public void increment(char x) throws Exception{
		if (x == 'G') g++;
		else if (x == 'A') a++;
		else if (x == 'T') t++;
		else if (x == 'C') c++;
		else throw new Exception("Unrecognized base "+x);
	}
	/**A,C,G,T,N,Del,Ins,FailBQ*/
	public String debug(){
		return (bpPosition+1)+"\t"+ref+"\t"+a+","+c+","+g+","+t+","+n+","+del+","+ins+","+failQual;
	}
	/**A,C,G,T,N,Del,Ins,FailBQ*/
	public String toString(){
		return a+","+c+","+g+","+t+","+n+","+del+","+ins+","+failQual;
	}
	/**A,C,G,T,N,Del,Ins,FailBQ
	 * Returns whether any counts were found.*/
	public boolean loadStringBuilderWithCounts(StringBuilder sb) {
		if (getTotalReadCoverage()==0) {
			sb.append("0,0,0,0,0,0,0,0");
			return false;
		}
		sb.append(a); sb.append(",");
		sb.append(c); sb.append(",");
		sb.append(g); sb.append(",");
		sb.append(t); sb.append(",");
		sb.append(n); sb.append(",");
		sb.append(del); sb.append(",");
		sb.append(ins); sb.append(",");
		sb.append(failQual);
		return true;
	}
	/**Returns the passing read coverage minus the max count base or indel */
	public double fetchNonGermlinePassingReadCount() {
		int max = a;
		if (c> max) max= c;
		if (g> max) max= g;
		if (t> max) max= t;
		int indel = ins+del;
		if (indel> max) max = indel;
		return getPassingReadCoverage() - max;	
	}
	
	public double getAlleleFreqPassingBackground() {
		return fetchNonGermlinePassingReadCount() / getPassingReadCoverage();
	}
	
	public double getAlleleFreqIndel() {
		return (double)(ins+del) / getPassingReadCoverage();
	}
	
	/**Tests whether GATCID is max allele count, e.g. germline*/
	public boolean isAlleleMaxCount(char allele) {
		char germlineSnv = getMaxReadCountBase();
		if (germlineSnv == allele) return true;
		return false;
	}
	
	/**Returns the base (GATC) or Indel (ID) with the max passing read count */
	public char getMaxReadCountBase() {
		int max = a;
		char base = 'A';
		if (c> max) {
			max= c;
			base = 'C';
		}
		if (g> max) {
			max= g;
			base = 'G';
		}
		if (t> max) {
			base = 'T';
			max = t;
		}
		if (ins > max) {
			base = 'I';
			max = ins;
		}
		if (del > max) {
			base = 'D';
			max = del;
		}
		return base;	
	}

	public int getSnvCount(char allele) throws IOException {
		if (allele == 'A') return a;
		if (allele == 'C') return c;
		if (allele == 'G') return g;
		if (allele == 'T') return t;
		throw new IOException("Unrecognized allele? "+allele);
	}
	
	public int getIndelCount(char allele) throws IOException {
		if (allele == 'D') return del;
		if (allele == 'I') return ins;
		throw new IOException("Unrecognized allele? "+allele);
	}
	
	public double getSnvAlleleFreq(char allele) throws IOException {
		double aTotal = getPassingReadCoverageSnv();
		if (aTotal == 0) return 0;
		double aCount = getSnvCount(allele);
		return aCount/aTotal;
	}
	public double getIndelAlleleFreq(char allele) throws IOException {
		double aTotal = getPassingReadCoverage();
		if (aTotal == 0) return 0;
		double aCount = getIndelCount(allele);
		return aCount/aTotal;
	}
}
