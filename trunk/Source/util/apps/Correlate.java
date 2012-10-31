package util.apps;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;
import util.bio.cluster.*;

/**
 * Calculates a Pearson correlation coefficient between all pairs of serializae float[].
 */
public class Correlate {
	
	//fields
	private File directory;
	private boolean cluster = false;
	private boolean collapseArrays = false;
	private File[] files;
	 
	//constructor
	public Correlate (String[] args){
		processArgs (args);
		//get floats
		files = IO.extractFiles(directory);
		//convert?
		File tempDir = null;
		if (collapseArrays){
			//make tempDirectory
			String rnd = Passwords.createRandowWord(8);
			tempDir = new File (files[0].getParentFile(), "TempDir_"+rnd);
			tempDir.mkdir();
			File[] convertedFiles = new File[files.length];
			for (int i=0; i< files.length; i++){
				float[][] oneD = (float[][])IO.fetchObject(files[i]);
				float[] one = Num.collapseFloatArray(oneD);
				convertedFiles[i] = new File (tempDir, files[i].getName());
				IO.saveObject(convertedFiles[i], one);
			}
			files = convertedFiles;
		}
		if (files.length < 2) Misc.printExit("\nLess than two files were found, aborting.\n");
		correlate();
		//cluster?
		if (cluster) new HierarchicalClustering(files);
		
		if (collapseArrays) {
			//move ClusterPlot.png
			File toMove = new File (tempDir, "ClusterPlot.png");
			File moved = new File (directory, "ClusterPlot.png");
			toMove.renameTo(moved);
			IO.deleteDirectory(tempDir);;
		}
	}
	
	public Correlate(File[] files, File hcName){
		this.files = files;
		correlate();
		new HierarchicalClustering (files, hcName);
	}
	
	//methods
	public void correlate() {
		//make hash fileName:ArrayList for R values
		LinkedHashMap hash = new LinkedHashMap(files.length);
		for (int i=0; i<files.length; i++){
			hash.put(files[i].getName(), new ArrayList());
		}
		
		//calc corr coef
		System.out.println("\nPairwise correlations:");
		for (int i=0; i<files.length; i++){
			int j = i+1;
			float[] one = (float[])IO.fetchObject(files[i]);
			ArrayList oneAl = (ArrayList)hash.get(files[i].getName());
			if (j<files.length) System.out.println(files[i].getName()+" vs ");
			for (; j<files.length; j++){
				float[] two = (float[])IO.fetchObject((files[j]));
				//print results
				double pcc = PearsonCorrelation.correlationCoefficient(one, two);
				String r = Num.formatNumber(pcc,3);
				System.out.println("\t"+files[j].getName()+"\t"+r+ "\t"+Num.formatPercentOneFraction(pcc*pcc));
				//load ArrayLists
				oneAl.add(r);
				ArrayList twoAl = (ArrayList)hash.get(files[j].getName());
				twoAl.add(r);
			}
			System.out.println();
		}
		
		//print hash
		System.out.println("File Name\tAssociated R Values");
		Iterator it = hash.keySet().iterator();
		while (it.hasNext()){
			String name = (String)it.next();
			ArrayList rs = (ArrayList)hash.get(name);
			System.out.println(name+"\t"+rs);
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new Correlate(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'd': directory = new File (args[i+1]); i++; break;
					case 'c': cluster = true; break;
					case 'a': collapseArrays = true; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (directory == null || directory.isDirectory() == false) Misc.printExit("\nError: cannot find your directory! "+directory);
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Correlate:    Nov 2008                            **\n" +
				"**************************************************************************************\n" +
				"Calculates all pair-wise Pearson correlation coefficients (r) and if indicated will\n" +
				"perform a hierarchical clustering on the files.\n\n"+				
				
				"Parameters:\n" +
				"-d The full path directory text containing serialized java float[] files (xxx.celp\n"+
				"      see CelProcessor app).\n"+
				"-a Files provided are float[][] files (xxx.cela) and need to be collapsed to float[]\n"+
				"-c Cluster files.\n\n" +
				
				"Example: java -Xmx256M -jar pathTo/T2/Apps/Correlate -d /Mango/PCels/ -c -a\n\n" +

				
		"**************************************************************************************\n");		
	}	
}
