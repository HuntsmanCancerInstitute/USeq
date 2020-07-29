package edu.utah.seq.vcf.fdr;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

class VCFFdrDataset {

	private static final String fdrInfo = "##INFO=<ID=dFDR,Number=1,Type=Float,Description=\"Estimated false discovery rate (-10Log10(Ave#Mock/#Real)) based on mock T/N contrasts, see the USeq VCFFdrEstimator application. 10=10%, 13=5%, 20=1%, 30=0.1%\">";
	private QualSortedVcf[] qualSortedVcfs = null;
	private float[] qualScores = null;
	private ArrayList<String> header = new ArrayList<String>();
	private File vcfFile;
	private double minFDR;

	public VCFFdrDataset (File vcfFile, double minFDR) {
		this.vcfFile = vcfFile;
		this.minFDR = minFDR;
		loadVcfRecords();
	}

	private void loadVcfRecords() {
		BufferedReader in = null;
		String line = null;
		try {
			in = IO.fetchBufferedReader(vcfFile);
			ArrayList<QualSortedVcf> quals = new ArrayList<QualSortedVcf>();
			boolean addFDRInfo = true;
			float order = 0;
			while ((line = in.readLine())!= null) {
				if (line.startsWith("#") == false) {
					//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
					String[] fields = Misc.TAB.split(line);
					Float qual = null;
					if (fields[5].equals(".") || fields[5].length() == 0) qual = new Float(0.0f);
					else qual = Float.parseFloat(fields[5]);
					quals.add( new QualSortedVcf(qual, fields, order++));
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
		float order;

		QualSortedVcf (float qual, String[] vcfFields, float order) {
			this.qual = qual;
			this.vcfFields = vcfFields;
			this.order = order;
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
		double priorFdr = 1000;
		float priorQual = -1;
		double numRealCounts = 0;
		double numMockCounts = 0;
		String priorFdrString = "dFDR=100;";
		ArrayList<QualSortedVcf> pass = new ArrayList<QualSortedVcf>();
		ArrayList<QualSortedVcf> fail = new ArrayList<QualSortedVcf>();
		for (int i=0; i< qualScores.length; i++) {
			if (priorQual != qualScores[i]) {
				priorQual = qualScores[i];
				numRealCounts = VCFFdrEstimator.countGreaterThanEqualTo(qualScores[i], qualScores);
				numMockCounts = vcfFdrEstimator.estimateMockCounts(qualScores[i]);
				//watch out for zero counts
				if (numRealCounts != 0 && numMockCounts != 0) {
					double fdr = numMockCounts/ numRealCounts;
					//only assign to prior if it's smaller, e.g. decreasingFDR dFDR
					if (fdr < priorFdr) {
						priorFdr = fdr;	
						priorFdrString = "dFDR=" + Num.formatNumberJustMax(Num.minus10log10(fdr), 2)+";";
					}
				}
			}

			qualSortedVcfs[i].vcfFields[7] = priorFdrString + qualSortedVcfs[i].vcfFields[7];
			//IO.pl("Qual "+priorQual+"  #Real "+numRealCounts+"   #NumMock "+numMockCounts+"   Phred "+priorPhredString +"\t"+(Misc.stringArrayToString(qualSortedVcfs[i].vcfFields, "\t")));
			if (priorFdr <= minFDR) pass.add(qualSortedVcfs[i]);
			else fail.add(qualSortedVcfs[i]);
		}
		sortSave(pass, outPass);
		sortSave(fail, outFail);
		IO.pl("\t"+pass.size()+"\t"+fail.size());
	}

	private void sortSave(ArrayList<QualSortedVcf> qsv, Gzipper out) throws IOException {
		if (qsv.size() == 0) {
			out.close();
			IO.deleteFile(out.getGzipFile());
		}
		QualSortedVcf[] toSort = new QualSortedVcf[qsv.size()];
		qsv.toArray(toSort);
		for (int i=0; i< toSort.length; i++) toSort[i].qual= toSort[i].order;
		Arrays.sort(toSort);
		for (int i=0; i< toSort.length; i++) out.println(Misc.stringArrayToString(toSort[i].vcfFields, "\t"));	
		out.close();
	}
}

