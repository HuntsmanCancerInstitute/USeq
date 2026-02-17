package util.gen;

import java.io.*;

/**For generating pvalues for intersection of lists. This could be genes or samples.
 * @author Nix
From https://rdrr.io/bioc/GeneOverlap/man/GeneOverlap.html , hypergeometric calc is faster than fishers exact

n = total number of genes
a = number of genes in set A
b = number of genes in set B
t = number of genes in common to A and B
p_upper = sum(dhyper(t:b, a, n - a, b)) for over-representation

Using built in cumulative
p_upper = 1 - phyper(t - 1, a, n - a, b) for over-representation
p_lower = phyper(t, a, n - a, b) for under-representation

n = 200
a = 70
b = 30
t = 10
sum(dhyper(t:b  , a , n - a   , b))
sum(dhyper(10:30, 70, 200 - 70, 30))
0.6561562

n = 20000
a = 500
b = 15
t = 10
sum(dhyper(t:b  , a , n - a   , b))
sum(dhyper(10:15, 500, 20000 - 500, 15))
2.342473e-13

underrepresentation
n = 100
a = 50
b = 60
t = 10
1 - phyper(t - 1, a, n - a, b)
1
phyper(t, a, n - a, b)
7.472794e-19


 */
public class IntersectListsHypergeometric {

	//fields
	private File tempDirectory;
	private File fullPathToR;

	public IntersectListsHypergeometric(File tempDirectory, File fullPathToR){
		this.tempDirectory = tempDirectory;
		this.fullPathToR = fullPathToR;
	}

	public static void main(String[] args) throws IOException {
		File tempDir = new File("/Users/u0028003/Downloads/Test");
		File r = new File("/usr/local/bin/R");
		IntersectListsHypergeometric ih = new IntersectListsHypergeometric(tempDir, r);
		int[][] nabt = {
				{200,70,30,10},
				{20000,500,15,10},
				{100, 50, 60, 10}
		};

		//over rep
		double[] pvals = ih.calculateOverRepresentationPValues(nabt);
		for (double p: pvals) IO.pl(p);
		// using sum and dhyper
		// 0.656156221606744
		// 2.34247276269585e-13
		// 0.999999
		// using phyper
		// 0.656156221606744
		// 2.34257058195908E-13
		// 1.0
		
		//under rep
		IO.pl();
		double[] underPvals = ih.calculateUnderRepresentationPValues(nabt);
		for (double p: underPvals) IO.pl(p);
		// 0.505902152540219
		// 0.999999999999997
		// 7.47279441121441E-19
		
		//both
		IO.pl();
		double[] overUnderPvals = ih.calculateOverUnderRepresentationPValues(nabt);
		for (int i=0; i< overUnderPvals.length; i++) {
			IO.pl(overUnderPvals[i++]+"\t"+overUnderPvals[i]);
		}
	}

	/**Returns the p-value for a single intersection test. Don't use for testing many intersections! Slow.
		n = total number of genes
		a = number of genes in set A
		b = number of genes in set B
		t = number of genes in common to A and B, the intersection
	 * @throws IOException */
	public double calculateOverRepresentationPValue(int n, int a, int b, int t) throws IOException{
		int[][] v = new int[1][4];
		v[0] = new int[]{n, a, b, t};
		double[] p = calculateOverRepresentationPValues(v);
		return p[0];
	}

	/**Returns the upper tail, over representation p-values for many intersection tests. Faster.
	1 n = total number of genes
	2 a = number of genes in set A
	3 b = number of genes in set B
	4 t = number of genes in common to A and B, the intersection
	 * @throws IOException */
	public double[] calculateOverRepresentationPValues(int[][] nabt) throws IOException {

		String rndWrd = Misc.getRandomString(10);
		File counts = writeOutFourCounts(nabt, rndWrd);
		File pvals = new File (tempDirectory, rndWrd+"_pvals.txt");
	
		File rOut = new File(tempDirectory, rndWrd+"_Script.txt.Rout");
		File scriptFile = new File(tempDirectory, rndWrd+"_Script.txt");

		//build R script
		StringBuilder script = new StringBuilder();
		script.append("counts = read.table('"+counts.getCanonicalPath()+"')\n");
		
		script.append("numRows = nrow(counts)\n");
		script.append("r = runif(numRows)\n");

		script.append("for (i in 1:numRows){\n");
		//script.append("    r[i] = sum(dhyper(counts[i,4]:counts[i,3]  , counts[i,2] , counts[i,1] - counts[i,2]   , counts[i,3]))\n");
		script.append("    r[i] = 1 - phyper(counts[i,4] - 1, counts[i,2], counts[i,1] - counts[i,2], counts[i,3])\n");
		script.append("}\n");

		script.append("write.table(r, file='"+pvals.getCanonicalPath()+"',row.names = FALSE, col.names = FALSE)\n");

		//write script
		IO.writeString(script.toString(), scriptFile);
		//make command
		String[] command = new String[] {
				fullPathToR.getCanonicalPath(),
				"CMD",
				"BATCH",
				"--no-save",
				"--no-restore",
				scriptFile.getCanonicalPath(),
				rOut.getCanonicalPath()};			
		//execute
		IO.executeCommandLine(command);
		//read in results 
		double[] pvalCalc = new double[nabt.length];

		if (pvals.exists() != false){
			String line;
			BufferedReader in = new BufferedReader ( new FileReader(pvals));
			int counter =0;
			while ((line=in.readLine()) != null){
				pvalCalc[counter++] = Double.parseDouble(line.trim());
			}
			in.close();
			if (counter != pvalCalc.length) throw new IOException("\nIncorrect length of R results for pval estimation. See "+pvals);
		}
		else throw new IOException("\nR results file doesn't exist. Check tempFiles to debug.\n");

		//cleanup temp files
		counts.delete();
		pvals.delete();
		rOut.delete();
		scriptFile.delete();

		return pvalCalc;
	}
	
	/**Returns the lower tail, under representation p-values for many intersection tests. Faster.
	1 n = total number of genes
	2 a = number of genes in set A
	3 b = number of genes in set B
	4 t = number of genes in common to A and B, the intersection
	 * @throws IOException */
	public double[] calculateUnderRepresentationPValues(int[][] nabt) throws IOException {

		String rndWrd = Misc.getRandomString(10);
		File counts = writeOutFourCounts(nabt, rndWrd);
		File pvals = new File (tempDirectory, rndWrd+"_pvals.txt");
	
		File rOut = new File(tempDirectory, rndWrd+"_Script.txt.Rout");
		File scriptFile = new File(tempDirectory, rndWrd+"_Script.txt");

		//build R script
		StringBuilder script = new StringBuilder();
		script.append("counts = read.table('"+counts.getCanonicalPath()+"')\n");
		
		script.append("numRows = nrow(counts)\n");
		script.append("r = runif(numRows)\n");

		script.append("for (i in 1:numRows){\n");
		//                 4  2  1   2  3
		//p_lower = phyper(t, a, n - a, b) for under-representation
		script.append("    r[i] = phyper(counts[i,4], counts[i,2], counts[i,1] - counts[i,2], counts[i,3])\n");
		script.append("}\n");

		script.append("write.table(r, file='"+pvals.getCanonicalPath()+"',row.names = FALSE, col.names = FALSE)\n");

		//write script
		IO.writeString(script.toString(), scriptFile);
		//make command
		String[] command = new String[] {
				fullPathToR.getCanonicalPath(),
				"CMD",
				"BATCH",
				"--no-save",
				"--no-restore",
				scriptFile.getCanonicalPath(),
				rOut.getCanonicalPath()};			
		//execute
		IO.executeCommandLine(command);
		//read in results 
		double[] pvalCalc = new double[nabt.length];

		if (pvals.exists() != false){
			String line;
			BufferedReader in = new BufferedReader ( new FileReader(pvals));
			int counter =0;
			while ((line=in.readLine()) != null){
				pvalCalc[counter++] = Double.parseDouble(line.trim());
			}
			in.close();
			if (counter != pvalCalc.length) throw new IOException("\nIncorrect length of R results for pval estimation. See "+pvals);
		}
		else throw new IOException("\nR results file doesn't exist. Check tempFiles to debug.\n");

		//cleanup temp files
		counts.delete();
		pvals.delete();
		rOut.delete();
		scriptFile.delete();

		return pvalCalc;
	}
	
	/**Returns the over and under representation p-values for intersection.
	1 n = total number of genes
	2 a = number of genes in set A
	3 b = number of genes in set B
	4 t = number of genes in common to A and B, the intersection
	 * @throws IOException */
	public double[] calculateOverUnderRepresentationPValues(int[][] nabt) throws IOException {

		String rndWrd = Misc.getRandomString(10);
		File counts = writeOutFourCounts(nabt, rndWrd);
		File pvals = new File (tempDirectory, rndWrd+"_pvals.txt");
	
		File rOut = new File(tempDirectory, rndWrd+"_Script.txt.Rout");
		File scriptFile = new File(tempDirectory, rndWrd+"_Script.txt");

		//build R script
		StringBuilder script = new StringBuilder();
		script.append("counts = read.table('"+counts.getCanonicalPath()+"')\n");
		
		script.append("numRows = nrow(counts)\n");
		script.append("r = runif(numRows*2)\n");
		script.append("counter = 1\n");
		script.append("for (i in 1:numRows){\n");
		script.append("    r[counter] = 1 - phyper(counts[i,4] - 1, counts[i,2], counts[i,1] - counts[i,2], counts[i,3])\n");
		script.append("    counter = counter + 1\n");
		script.append("    r[counter] = phyper(counts[i,4], counts[i,2], counts[i,1] - counts[i,2], counts[i,3])\n");
		script.append("    counter = counter + 1\n");
		script.append("}\n");

		script.append("write.table(r, file='"+pvals.getCanonicalPath()+"',row.names = FALSE, col.names = FALSE)\n");

		//write script
		IO.writeString(script.toString(), scriptFile);
		//make command
		String[] command = new String[] {
				fullPathToR.getCanonicalPath(),
				"CMD",
				"BATCH",
				"--no-save",
				"--no-restore",
				scriptFile.getCanonicalPath(),
				rOut.getCanonicalPath()};			
		//execute
		IO.executeCommandLine(command);
		//read in results 
		double[] pvalCalc = new double[nabt.length * 2];

		if (pvals.exists() != false){
			String line;
			BufferedReader in = new BufferedReader ( new FileReader(pvals));
			int counter =0;
			while ((line=in.readLine()) != null){
				pvalCalc[counter++] = Double.parseDouble(line.trim());
			}
			in.close();
			if (counter != pvalCalc.length) throw new IOException("\nIncorrect length of R results for pval estimation. See "+pvals);
		}
		else throw new IOException("\nR results file doesn't exist. Check tempFiles to debug.\n");

		//cleanup temp files
		counts.delete();
		pvals.delete();
		rOut.delete();
		scriptFile.delete();

		return pvalCalc;
	}



	private File writeOutFourCounts(int[][] nabt, String randomWord) throws IOException {
		File tempData = new File (tempDirectory, randomWord+"_countData.txt");
		PrintWriter out = new PrintWriter( new FileWriter(tempData));
		for (int[] c: nabt) {
			out.print(c[0]);
			out.print("\t");
			out.print(c[1]);
			out.print("\t");
			out.print(c[2]);
			out.print("\t");
			out.println(c[3]);
		}
		out.close();
		return tempData;
	}

}
