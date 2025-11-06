package util.gen;

import java.io.*;

/**For generating pvalues for intersection of lists.
 * @author Nix
From https://rdrr.io/bioc/GeneOverlap/man/GeneOverlap.html , hypergeometric calc is faster than fishers exact

n = total number of genes
a = number of genes in set A
b = number of genes in set B
t = number of genes in common to A and B
pval = sum(dhyper(t:b, a, n - a, b))

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
				{20000,500,15,10}
		};

		double[] pvals = ih.calculatePValues(nabt);
		for (double p: pvals) IO.pl(p);

		//0.656156221606744
		//2.34247276269585e-13
	}

	/**Returns the p-value for a single intersection test. Don't use for testing many intersections! Slow.
		n = total number of genes
		a = number of genes in set A
		b = number of genes in set B
		t = number of genes in common to A and B, the intersection
	 * @throws IOException */
	public double calculatePValue(int n, int a, int b, int t) throws IOException{
		int[][] v = new int[1][4];
		v[0] = new int[]{n, a, b, t};
		double[] p = calculatePValues(v);
		return p[0];
	}

	/**Returns the p-values for a many intersection tests. Faster.
	n = total number of genes
	a = number of genes in set A
	b = number of genes in set B
	t = number of genes in common to A and B, the intersection
	 * @throws IOException */
	public double[] calculatePValues(int[][] nabt) throws IOException {

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
		script.append("    r[i] = sum(dhyper(counts[i,4]:counts[i,3]  , counts[i,2] , counts[i,1] - counts[i,2]   , counts[i,3]))\n");
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
