package edu.utah.seq.vcf.loh;

import util.gen.Num;

public class PassingHet {
	
	private String chr;
	private int position;
	private double tumorCount;
	private double normalCount;
	private String tabixCoor;
	
	public PassingHet(String chr, int position, int tumorCount, int normalCount) {
		this.chr = chr;
		this.position = position;
		this.tumorCount = tumorCount;
		this.normalCount = normalCount;
		tabixCoor = chr+ ":"+ position+ "-"+ position;
	}
	
	public double fetchNormalizedLog2TNRatio (double tumorScalar, double normalScalar) {
		double tumCount = tumorScalar * tumorCount;
		double normCount = normalScalar * normalCount;
		return Num.log2(tumCount/normCount);
	}
	
	public double fetchNormalizedTNRatio (double tumorScalar, double normalScalar) {
		double tumCount = tumorScalar * tumorCount;
		double normCount = normalScalar * normalCount;
		return tumCount/normCount;
	}
	
	public String fetchBedLine(double tumorScalar, double normalScalar) {
		double rto = fetchNormalizedLog2TNRatio(tumorScalar, normalScalar);
		return chr+ "\t"+ (position-1)+ "\t"+ position+ "\t" + tumorCount+"x"+
				Num.formatNumber(tumorScalar, 3)+"/"+normalCount+"x"+ Num.formatNumber(normalScalar, 3)+ "\t"+ Num.formatNumber(rto, 3)+"\t.";
	}
	
	public double fetchNormal1KScalar() {
		return 1000.0 / normalCount;
	}
	
	public double fetchTumor1KScalar() {
		return 1000.0 / tumorCount;
	}

	public String getTabixCoor() {
		return tabixCoor;
	}

	public String getChr() {
		return chr;
	}

	public int getPosition() {
		return position;
	}

}
