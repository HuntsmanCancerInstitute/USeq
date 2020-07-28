package edu.utah.seq.vcf.fdr;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

class VCFFdrDataset {

	private static final String fdrInfo = "##INFO=<ID=dFDR,Number=1,Type=Float,Description=\"-10Log10(FDR) based on mock T/N contrasts; 20=0.01, 13=0.05, 10=0.1; see the USeq VCFFdrEstimator application.\">";
	private QualSortedVcf[] qualSortedVcfs = null;
	private float[] qualScores = null;
	private ArrayList<String> header = new ArrayList<String>();
	private File vcfFile;
	private double minFDR;
	private int numPass = 0;
	private int numFail = 0;
	private int numInMock = 0;
	private HashSet<String> mockVars;
	private boolean excludeMockMatches;

	public VCFFdrDataset (File vcfFile, double minFDR, HashSet<String> mockVars, boolean excludeMockMatches) {
		this.vcfFile = vcfFile;
		this.minFDR = minFDR;
		this.mockVars = mockVars;
		this.excludeMockMatches = excludeMockMatches;
		loadVcfRecords();
	}

	private void loadVcfRecords() {
		BufferedReader in = null;
		String line = null;
		try {
			in = IO.fetchBufferedReader(vcfFile);
			ArrayList<QualSortedVcf> quals = new ArrayList<QualSortedVcf>();
			boolean addFDRInfo = true;
			while ((line = in.readLine())!= null) {
				if (line.startsWith("#") == false) {
					//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
					String[] fields = Misc.TAB.split(line);
					Float qual = null;
					if (fields[5].equals(".") || fields[5].length() == 0) qual = new Float(0.0f);
					else qual = Float.parseFloat(fields[5]);
					quals.add( new QualSortedVcf(qual, fields));
				}
				else {
					//add header?
					if (addFDRInfo && line.startsWith("##INFO=")) {
						header.add(fdrInfo);
						addFDRInfo = false;
					}
					header.add(line);
				}
			}

			qualSortedVcfs = new QualSortedVcf[quals.size()];
			quals.toArray(qualSortedVcfs);
			Arrays.sort(qualSortedVcfs);
			qualScores = new float[qualSortedVcfs.length];
			for (int i=0; i< qualScores.length; i++) qualScores[i] = qualSortedVcfs[i].qual;

		} catch (Exception e) {
			if (in != null) IO.closeNoException(in);
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing QUAL score from "+vcfFile+" for VCF record "+line);
		} finally {
			if (in != null) IO.closeNoException(in);
		}
	}

	private class QualSortedVcf implements Comparable<QualSortedVcf> {
		float qual;
		String[] vcfFields;
		boolean inMock = false;

		QualSortedVcf (float qual, String[] vcfFields) {
			this.qual = qual;
			this.vcfFields = vcfFields;
			String coor = vcfFields[0]+"_"+vcfFields[1]+"_"+vcfFields[3]+"_"+vcfFields[4];
			if (mockVars.contains(coor)) {
				if (excludeMockMatches) inMock = true;
				numInMock++;
			}
		}

		public int compareTo(QualSortedVcf other) {
			//sort by qual
			if (qual < other.qual) return -1;
			if (qual > other.qual) return 1;
			return 0;
		}
	}
	
	public void addHeader(Gzipper outPass, Gzipper outFail) throws IOException {
		int num = header.size();
		for (int i=0; i< num; i++) {
			String h = header.get(i);
			outPass.println(h);
			outFail.println(h);
		}
	}

	public void estimateFdrs(Gzipper outPass, Gzipper outFail, VCFFdrEstimator vcfFdrEstimator) throws IOException {
		double priorPhred = -1;
		float priorQual = -1;
		double numRealCounts = 0;
		double numMockCounts = 0;
		String priorPhredString = null;
		for (int i=0; i< qualScores.length; i++) {
			if (priorQual != qualScores[i]) {
				priorQual = qualScores[i];
				numRealCounts = VCFFdrEstimator.countGreaterThanEqualTo(qualScores[i], qualScores);
				numMockCounts = vcfFdrEstimator.estimateMockCounts(qualScores[i]);
				//watch out for zero counts
				if (numRealCounts != 0 && numMockCounts != 0) {
					double fdr = numMockCounts/ numRealCounts;
					double phred = Num.minus10log10(fdr);
					//only assign to prior if it's better
					if (phred > priorPhred) {
						priorPhred = phred;	
						priorPhredString = "dFDR=" + Num.formatNumber(priorPhred, 2)+";";
					}
				}
			}

			qualSortedVcfs[i].vcfFields[7] = priorPhredString + qualSortedVcfs[i].vcfFields[7];
			//IO.pl("Qual "+priorQual+"  #Real "+numRealCounts+"   #NumMock "+numMockCounts+"   Phred "+priorPhredString +"\t"+(Misc.stringArrayToString(qualSortedVcfs[i].vcfFields, "\t")));
			if (qualSortedVcfs[i].inMock == false && priorPhred >= minFDR) {
				outPass.println(Misc.stringArrayToString(qualSortedVcfs[i].vcfFields, "\t"));
				numPass++;
			}
			else {
				outFail.println(Misc.stringArrayToString(qualSortedVcfs[i].vcfFields, "\t"));
				numFail++;
			}
		}
		IO.pl("\t"+numPass+"\t"+numFail+"\t"+numInMock);
	}

	public int getNumPass() {
		return numPass;
	}

	public int getNumFail() {
		return numFail;
	}
}

