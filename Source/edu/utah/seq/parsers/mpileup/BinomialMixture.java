package edu.utah.seq.parsers.mpileup;

import java.io.File;
import java.util.ArrayList;

import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;
import util.gen.Passwords;

/**Wrapper for an algorithm by Alun Thomas, University of Utah to call confidence on non reference base observations in a tumor/test sample given a panel of normals.
 * Uses an EM algorithm to fit a binomial mixture of normal observations that may contain outlier error prone normals.*/
public class BinomialMixture {
	
	//params
	private File r;
	private File tempDir;
	private boolean deleteTempFiles = true;
	
	public BinomialMixture (File rExe, File tempDir, boolean deleteTempFiles){
		this.r = rExe;
		this.tempDir = tempDir;
		this.deleteTempFiles = deleteTempFiles;
	}
	
	/**Call (repeatedly) to fetch -10Log10(AdjPVals) using Alun Thomas' R alogorithm to estimate the probability that a tumor Allele Frequency observations nonRef/(ref+nonRef) comes 
	 * from a mixed binomial distribution of normals where some may be extra noisy.
	 * 
	 * The data file should contain two rows per test, the first is the non reference base observations, the second the total # observations.
	 * The first column of data should come from the tumor sample, all subsequent columns, the normals, 
	 * e.g. 1 tumor sample (35/2865), 14 normals (0/6046, 4/1573, etc.)
	 * 35 0 4 1 0 0 0 0 0 2 1 0 1 1 5
	 * 2865 6046 1573 3869  868  970  952  962  915  976  825  821 5119 3914 6272
	 * 
	 * The results contain Phred transformed, Benjamini-Hochberg adjusted pvals, e.g. -10Log10(adjPvals)
	 * For INF values, the max score is assigned.
	 * Returns null if an error was thrown.*/
	public double[] estimatePValues(File dataFile){
		try {
			String rndWord = Passwords.createRandowWord(6);
			
			//build cmd script
			File rResults = new File (tempDir, rndWord+"_RResults.txt");
			String cmd = fetchRCmd (dataFile, rResults);
			File rScriptFile = new File(tempDir, rndWord+"_RScript.R");
			IO.writeString(cmd, rScriptFile);

			//make command
			File rOutFile = new File (tempDir, rndWord+"_RScript.out");

			String[] command = new String[] {
					r.getCanonicalPath(),
					"CMD",
					"BATCH",
					"--no-save",
					"--no-restore",
					rScriptFile.getCanonicalPath(),
					rOutFile.getCanonicalPath()};			
			//execute
			String[] output = IO.executeViaProcessBuilder(command, true);
			if (output == null || rResults.exists() == false) return null;

			//load results
			double[] values = Num.loadDoubles(rResults);
			if (values == null) return null;

			//cleanup if no errors
			if (deleteTempFiles){
				rResults.deleteOnExit();
				rScriptFile.deleteOnExit();
				rOutFile.deleteOnExit();
			}

			return values;
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}
	}
	

	private String fetchRCmd (File dataFile, File resultsFile){
		ArrayList<String> s = new ArrayList<String>();
		
		s.add("fitbinmix = function(x,n,steps = 100) {");
		s.add("	p = sum(x)/sum(n)");
		s.add("	loglike0 = sum(log(dbinom(x,n,p)))");
		s.add("	q = 2*p");
		s.add("	y = dbinom(x,n,p) / ( dbinom(x,n,p) + dbinom(x,n,q) )");
		s.add("	for (i in 1:steps) {");
		s.add("		p = sum(y*x) / sum(y*n)");
		s.add("		q = sum((1-y)*x) / sum((1-y)*n)");
		s.add("		y = dbinom(x,n,p) / ( dbinom(x,n,p) + dbinom(x,n,q) )");
		s.add("	}");
		s.add("	r = mean(y)");
		s.add("	loglike1 = sum ( log ( r*dbinom(x,n,p) + (1-r)*dbinom(x,n,q) ) )");
		s.add("	llratio = 2 * (loglike1 - loglike0)");
		s.add("	pvalue = pchisq(llratio,2,lower.tail=F)");
		s.add("	list(pvalue=pvalue, llratio=llratio, p=p, q=q, r=r, y=y)");
		s.add("}");
		s.add("findnonrefs = function(x,n,xnorm,nnorm,steps=100){");
		s.add("	mix = fitbinmix(xnorm,nnorm,steps)");
		s.add("	if (!is.na(mix$pvalue) && mix$pvalue < 0.05){");
		s.add("		p = mix$p");
		s.add("	}");
		s.add("	else {");
		s.add("		p = sum(xnorm)/sum(nnorm)");
		s.add("	}");
		s.add("    pvalues = pbinom(x-1,n,p,lower.tail=F)");
		s.add("	pvalues");
		s.add("}");

		//working code to load the data table, NAs are added to ragged rows
		s.add("nr = read.table('"+ dataFile + "', header=F, fill=T)");
		s.add("numberRows = nrow(nr)");
		s.add("p = vector('numeric', length = numberRows/2)");
		s.add("counter=1");

		//for each variant
		s.add("for (i in seq(1,numberRows,2)){");
		s.add("	nonRef = nr[i,]");
		s.add("	nonRef = nonRef[!is.na(nonRef)]");
		s.add("	total = nr[i+1,]");
		s.add("	total = total[!is.na(total)]");
		s.add("	len = length(total)");
		s.add("	p[counter] = findnonrefs(nonRef[1], total[1], nonRef[2:len], total[2:len])");
		s.add("	counter = counter+1");
		s.add("}");

		// Adjust pvals with Benjamini-Hochberg
		s.add("adjP = p.adjust(p, method='BH')");

		// Phred the adjPvals
		s.add("finalP = -10*log10(adjP)");
		s.add("write(finalP, file='"+ resultsFile + "', sep='\t', ncolumns = 1)");

		return Misc.stringArrayListToString(s, "\n");
	}
	
}
