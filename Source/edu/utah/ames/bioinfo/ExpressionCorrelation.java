package edu.utah.ames.bioinfo;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Misc;
import util.gen.Num;
import util.gen.PearsonCorrelation;

/**
 * 
 * @author darren.ames@hci.utah.edu
 *
 */

public class ExpressionCorrelation {

	//fields
	private String expressionFile = "/Users/darren/Desktop/10318R_146KOvsWT_2K.txt";
	private ArrayList<Gene> gL = new ArrayList<Gene>();
	private ArrayList<Double> cL1 = new ArrayList<Double>();
	private ArrayList<Double> cL2 = new ArrayList<Double>();
	private double xTot=0;
	private double yTot=0;
	private double sqrXTot =0;
	private double sqrYTot =0;
	private double xYTot=0;
	private double totalN = 0;

	//constructor
	public ExpressionCorrelation(String[] args) {
		processArgs(args);
	}

	public static void main(String[] args) {
		//check for args
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		System.out.println("Starting...");
		ExpressionCorrelation ec = new ExpressionCorrelation(args);
		ec.parseExpressionFile();
		ec.calculateCorrelation();
		System.out.println("Finished!");
	}

	/**
	 * Use this method to sequentially add multiple pairs for correlation
	 * @param x
	 * @param y
	 */
	public void addPairsToCorrelate (double[] x, double[] y){
		double N = x.length;
		for (int i=0; i<N; i++){
			xTot += x[i];
			yTot += y[i];
			sqrXTot += Math.pow(x[i],2);
			sqrYTot += Math.pow(y[i],2);
			xYTot += (x[i] * y[i]);
		} 
		totalN += N;
	}

	/**
	 * Following loading of this class using the addPairsToCorrelate(),
	 * call this method to calculate to correlation coefficient (R) on the total.
	 * Returns null if an error is encountered.
	 * @return
	 */
	public Double calculateAdditivePairCorrelation(){
		double top = (totalN * xYTot) - (xTot * yTot);
		double botLeft = Math.sqrt( (totalN * sqrXTot) - Math.pow(xTot,2) );
		double botRight = Math.sqrt( (totalN * sqrYTot) - Math.pow(yTot,2) );
		double bot = botLeft*botRight;
		if (bot == 0) {
			System.err.println("\nERROR calculating correlation:");
			System.err.println("\tTop "+top);
			System.err.println("\tBotL "+botLeft);
			System.err.println("\tBotR "+botRight);
			System.err.println("\tBot "+bot);
			return null;
		}
		return new Double(top/bot);
	}

	/**
	 * Calculates the correlation coefficient on the 2 conditions and prints output
	 * correlation file. 
	 * @param gL
	 */
	public void calculateCorrelation() {
		try {
			//output file name
			String outName = null;
			String dirName = null;

			//fetch input path
			Pattern p2 = Pattern.compile("/.+/");
			Matcher m2 = p2.matcher(expressionFile.toString());
			if (m2.find()) {
				dirName = m2.group();
			}

			//strip path from input file
			int index = expressionFile.lastIndexOf(File.separatorChar);
			String name = expressionFile.substring(index+1);

			//strip extension from input file
			if (name.indexOf(".") > 0) {
				outName = name.substring(0, name.lastIndexOf("."));
			}

			//good output		
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(dirName, outName + ".corr.txt")));
			out.write("Condition 1\tCondition 2\n");	

			//iterate over the arrays
			for (int i = 0; i < gL.size(); i++) {
				for (int j = i+1; j < gL.size(); j++) {

					//check for invalid input
					if (!gL.get(i).isValid(0) || !gL.get(i).isValid(1) || !gL.get(j).isValid(0) || !gL.get(j).isValid(1)) {
						continue;
					}
					//compare all pairwise comparisons and calculate correlation coefficients
					//double cor1 = Num.correlationCoefficient(gL.get(i).getArray(0), gL.get(j).getArray(0));
					//double cor2 = Num.correlationCoefficient(gL.get(i).getArray(1), gL.get(j).getArray(1));

					double cor1 = PearsonCorrelation.correlationCoefficient(gL.get(i).getArray(0), gL.get(j).getArray(0));
					double cor2 = PearsonCorrelation.correlationCoefficient(gL.get(i).getArray(1), gL.get(j).getArray(1));

					//check for corr coef errors
					if (cor1 == -2 || cor2 == -2) {
						System.out.println(Arrays.toString(gL.get(i).getArray(0)));
						System.out.println(Arrays.toString(gL.get(i).getArray(1)));
						System.out.println(Arrays.toString(gL.get(j).getArray(0)));
						System.out.println(Arrays.toString(gL.get(j).getArray(1)));
						System.out.println(i + " " + j);
						System.exit(0);
					}

					//add squared correlation coefficients to the array lists
					//cL1.add(cor1);
					//cL2.add(cor2);

					//write output
					out.write((cor1*cor1) + "\t" + (cor2*cor2) + "\n");
				}
			}
			out.close();
		}
		//catch exceptions
		catch (FileNotFoundException ioe) {
			System.out.println("File not found.");
			System.exit(0);
		}
		catch (IOException e) {
			System.out.println("Error reading file.");
			System.exit(0);
		}
	}

	/**
	 * Method to parse the FPKM values from an output file from USeq's DRDS app
	 * @param s
	 * @return
	 */
	public ArrayList<ArrayList<Integer>> parseString(String s) {
		//split contents of tab-delimited file 
		String[] headerList = s.split("\t");

		//array list for conditions
		ArrayList<ArrayList<Integer>> conds = new ArrayList<ArrayList<Integer>>();
		//array list for counts within an individual condition
		ArrayList<Integer> c = new ArrayList<Integer>();

		//iterate over the header to find condition-specific FPKM columns
		for (int i = 0; i < headerList.length; i++) {
			//patterns to match FPKM fields in header
			String thing1 = "FPKM_(.+?)";
			String thing2 = "FPKM_(.+?)0";
			Pattern p1 = Pattern.compile(thing1);
			Pattern p2 = Pattern.compile(thing2);
			Matcher m1 = p1.matcher(headerList[i]);
			Matcher m2 = p2.matcher(headerList[i]);

			//add first FPKM condition
			if (m1.matches()) {
				if (!m2.matches()) {
					c.add(i);
				}
				else {
					//add new conditions
					if (c.size() > 0) {
						conds.add(c);
					}
					c = new ArrayList<Integer>();
					c.add(i);
				}
			}
		}
		conds.add(c);

		for (ArrayList<Integer> cond : conds) {
			for (Integer w : cond) {
			}
		}
		return conds;
	}

	/**
	 * Method to read in the output expression file from USeq's DRDS
	 */
	public void parseExpressionFile() {

		try {
			//create buffered reader to read in the file
			BufferedReader br = new BufferedReader(new FileReader(expressionFile));

			String line = br.readLine(); 

			ArrayList<ArrayList<Integer>> conds = this.parseString(line);

			while ((line = br.readLine()) != null) {
				addToArrayList(line, conds);
			}

			//close the reader
			br.close();
		}
		//catch exceptions
		catch (FileNotFoundException e) {
			System.out.println("File not found.");
			System.exit(0);
		}
		catch (IOException e) {
			System.out.println("Error reading file.");
			System.exit(0);
		}
	}

	/**
	 * 
	 * @param line
	 * @param conds
	 */
	public void addToArrayList(String line, ArrayList<ArrayList<Integer>> conds) {

		Gene g = new Gene();
		String[] items = line.split("\t");
		for (int i = 0; i < items.length; i++) {

			//get condition 1
			if (conds.get(0).contains(i)) {
				g.addIndex(0, Double.parseDouble(items[i]));
			}
			//get condition 2
			else if (conds.get(1).contains(i)) {
				g.addIndex(1, Double.parseDouble(items[i]));
			}
		}
		gL.add(g);
	}

	/**
	 * This method will process each argument and assign new variables.
	 * @param args
	 */
	public void processArgs(String[] args) {

		Pattern pat = Pattern.compile("-[a-z]");
		String programArgs = Misc.stringArrayToString(args, ",");
		boolean verbose = false;
		if (verbose) System.out.println("\nArguments: " + programArgs + "\n");
		for (int i = 0; i < args.length; i++) {
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()) {
				char test = args[i].charAt(1);
				try {
					switch (test) {
					case 'f': expressionFile = new String(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem--unknown option used!" + mat.group());
					}
				}
				catch (Exception e) {
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -" + test + "\n");
				}
			}
		}
	}

	public static void printDocs() {
		System.out.println("\n" +
				"**********************************************************************************\n" +
				"**                       ExpressionCorrelation: July 2013                       **\n" +
				"**********************************************************************************\n" + 
				"" +
				"\nParameters: \n\n" +
				" -f input expression file\n" +

				"Usage:\n\n" +
				"java -jar pathTo/ExpressionCorrelation.jar -f inputFile.xls\n\n" +
				"Questions or comments? Contact: darren.ames@hci.utah.edu\n" +
				"**********************************************************************************\n");
	}
}
