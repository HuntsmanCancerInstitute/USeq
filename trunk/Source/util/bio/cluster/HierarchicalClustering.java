package util.bio.cluster;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import util.gen.*;

/**
 * Calculates a Pearson correlation coefficient between all pairs of serialized float[]s, uses hierarchical clustering to
 * group like arrays. Displays on screen and saves as a png.
 */
public class HierarchicalClustering {

	//fields
	private File[] celpFiles;
	private ArrayList clusters = new ArrayList();
	private HashMap correlationCoefficients = new HashMap();
	private Cluster megaCluster;
	private double minimalCorrelationCoefficient = 0.75;
	private boolean savePNG = true;
	private boolean displayClusterPlot = true;
	private File pngResultFile;
	private String title  = "";
	private StringBuffer results = new StringBuffer();
	private int maxSizeToFlagCluster = 1;
	private HashSet badFileNumbers = new HashSet();

	/**For stand alone use.*/
	public HierarchicalClustering (String[] args) {
		processArgs(args);
		cluster();
		System.out.println(results);
	}

	/**For integration with CelFileQualityControl App.*/
	public HierarchicalClustering (double minRValue, int maxSizeToFlagCluster, StringBuffer results, File pngResultFile, boolean displayClusterPlot){
		//set params
		minimalCorrelationCoefficient = minRValue;
		this.results = results;
		this.pngResultFile = pngResultFile;
		this.maxSizeToFlagCluster = maxSizeToFlagCluster;
		this.displayClusterPlot = displayClusterPlot;
	}

	/**For integration with Correlate App.  Fires stand alone settings using the provided serialized float[] files.*/
	public HierarchicalClustering (File[] floatArrayFiles){
		this.celpFiles =floatArrayFiles;
		displayClusterPlot = false;
		cluster();
	}

	/**For integration with CelProcessor.*/
	public HierarchicalClustering (File[] floatArrayFiles, File pngResultFile){
		this.pngResultFile = pngResultFile;
		this.celpFiles =floatArrayFiles;
		displayClusterPlot = false;
		cluster();
	}

	public void cluster(){
		try {
			results.append("\nLaunching ...\n");
			//cluster files
			clusterNormalizedFloatArrayFiles();
			//break out clusters
			results.append("Ordering...\n");
			breakOutClusters();
			//display cluster on screen?
			if (savePNG && pngResultFile == null) pngResultFile = new File(celpFiles[0].getParentFile(),title+"ClusterPlot.png");
			if (displayClusterPlot) new ClusterDrawFrame(this);
			else {
				//kill the need to display
				System.setProperty("java.awt.headless","true");
				new ClusterDrawPanel(this);
			}
		} catch (Exception e){
			results.append("\nError: problem with making cluster HC plot.\n");
			System.out.println("\nError: problem with making cluster HC plot.\n");
			e.printStackTrace(); 
		}
		results.append("\nDone!\n");
	}

	/**Performs actual clustering for CelFileQualityControl app.*/
	public File[] clusterVirtualCelFiles(File[] virtualCelFiles, int[][] pmCoordinates){
		File[] baddies = clusterCelaFiles(virtualCelFiles, pmCoordinates);
		breakOutClusters();
		//display cluster on screen?
		if (savePNG && pngResultFile == null) pngResultFile = new File(virtualCelFiles[0].getParentFile(),title+"ClusterPlot.png");
		if (displayClusterPlot) new ClusterDrawFrame(this);
		else {
			//kill the need to display
			System.setProperty("java.awt.headless","true");
			new ClusterDrawPanel(this);
		}
		return baddies;
	}

	/**Clusters cela files, extracts PM intensities, median normalizes. Returns the names of virtual cel files with
	 * bad correlation coefficients.*/
	public File[] clusterCelaFiles(File[] virtualCelFiles, int[][] pmCoordinates){
		Arrays.sort(virtualCelFiles);
		//get trimmed file names
		String[] trimmedNames = IO.getTruncatedNames(virtualCelFiles);
		//make temp directory and copy files into it, also create ArrayList to hold Clusters
		File tempDir = null;
		results.append("\nMaking temp files...\n");
		tempDir = new File(virtualCelFiles[0].getParent(),"HRClusterTempDir");
		if (tempDir.exists()) IO.deleteDirectory(tempDir);
		tempDir.mkdir();
		for (int i=0; i<virtualCelFiles.length; i++){
			File x = new File (tempDir,(i+1)+"");
			//fetch virtual cel matrix
			float[][] virtualCel = (float[][])IO.fetchObject(virtualCelFiles[i]);
			//extract and median normalize pm intensities
			float[] pmValues = Num.fetchMatrixValues(virtualCel, pmCoordinates);
			pmValues = Num.medianNormalize(pmValues, 100);
			//save in temp directory
			IO.saveObject(x, pmValues);
			//make VCluster
			Cluster cluster = new Cluster(x,trimmedNames[i]);
			clusters.add(cluster);
			results.append("\t"+x.getName()+"\t"+virtualCelFiles[i].getName()+"\n");
		}

		//run thru clusters making all pairwise corr coefs, merging best two, saving badFiles
		results.append("\nClustering...\nR Squared * 100\tCluster 1 \t Cluster 2\t Merged Cluster\n");
		while (clusters.size()>1){
			double[] ccBestPair = findBestPair();
			//merge best
			mergeClusters(ccBestPair);
		}

		megaCluster = (Cluster)clusters.get(0);

		//delete temp dir
		IO.deleteDirectory(tempDir);

		//convert baddies into File[]
		if (badFileNumbers.size() == 0) return null;
		File[] badFiles = new File[badFileNumbers.size()];
		int index = 0;
		Iterator it = badFileNumbers.iterator();
		while (it.hasNext()){
			String number = (String)it.next();
			int fileIndex = Integer.parseInt(number) - 1;
			badFiles[index] = virtualCelFiles[fileIndex];
			index++;
		}
		return badFiles;
	}


	/**Clusters float[] files, does not median normalize*/
	public void clusterNormalizedFloatArrayFiles(){
		Arrays.sort(celpFiles);
		//get trimmed file names
		String[] trimmedNames = IO.getTruncatedNames(celpFiles);
		//make temp directory and copy files into it, also create ArrayList to hold Clusters
		File tempDir = null;
		results.append("\nMaking temp files...\n");
		String random = Passwords.createRandowWord(10);
		tempDir = new File(celpFiles[0].getParent(),"HRClusterTempDir_"+random);
		if (tempDir.exists()) IO.deleteDirectory(tempDir);
		tempDir.mkdir();
		for (int i=0; i<celpFiles.length; i++){
			File x = new File (tempDir,(i+1)+"");
			IO.copyViaFileChannel(celpFiles[i],x);
			clusters.add(new Cluster(x,trimmedNames[i]));
			results.append("\t"+x.getName()+"\t"+celpFiles[i].getName()+"\n");
		}

		//run thru clusters making all pairwise corr coefs, merging best two
		results.append("\nClustering...\n");
		while (clusters.size()>1){
			double[] ccBestPair = findBestPair ();
			//merge best
			mergeClusters(ccBestPair);
		}
		megaCluster = (Cluster)clusters.get(0);

		//delete temp dir
		IO.deleteDirectory(tempDir);
	}

	public void breakOutClusters() {
		//order clusters by scanning and spliting parents
		int counter = 10000;
		while (counter-- > 0){ 
			boolean noParents = true;
			for (int i=0; i<clusters.size(); i++){
				Cluster parentCluster = (Cluster)clusters.get(i);
				//parent found?
				if (parentCluster.getTotalClusters() != 1){
					noParents = false;
					//fetch top and bottom
					Cluster top = parentCluster.getTopCluster();
					Cluster bottom = parentCluster.getBottomCluster();
					//replace parent with children
					clusters.remove(i);
					clusters.add(i, bottom);
					clusters.add(i,top);
				}
			}
			if (noParents){
				break;
			}
		}
	}	

	/**Merges Clusters from best pair notes whether they have a bad cc.*/
	public void mergeClusters(double[] ccBestPair){
		double cc = ccBestPair[0];
		int indexFirst = (int)ccBestPair[1];
		int indexSecond = (int)ccBestPair[2];
		Cluster one = (Cluster)clusters.get(indexFirst);
		Cluster two = (Cluster)clusters.get(indexSecond);

		//check for bad cc and set
		//same number names of bad clusters
		if (cc == 1 || cc< minimalCorrelationCoefficient) {
			String[] ones = one.getFloatArrayFile().getName().split("_");
			String[] twos = two.getFloatArrayFile().getName().split("_");
			//save baddies
			//both with one file, add both
			if (twos.length == ones.length && ones.length ==1) {
				badFileNumbers.add(ones[0]);
				badFileNumbers.add(twos[0]);
			}
			else {
				if (twos.length < ones.length) ones = twos;
				if (ones.length < maxSizeToFlagCluster ){
					for (int x=0; x< ones.length; x++){
						badFileNumbers.add(ones[x]);
					}
				}
			}
		}
		//make new cluster
		Cluster combo = new Cluster (one, two, cc, results);
		//remove old 
		if (indexFirst > indexSecond) {
			clusters.remove(indexFirst);
			clusters.remove(indexSecond);
		}
		else {
			clusters.remove(indexSecond);
			clusters.remove(indexFirst);
		}
		//add new
		clusters.add(combo);
	}

	/**Takes an ArrayList of VCluster and finds the best pair based on a 
	 * PearsonCorrelationCoefficient.  Returns a double[3]{cc, indexOne, indexTwo}.
	 * The HashMap contains calculated ccs.  These will be used if the same comparison is made again.*/
	public double[] findBestPair (){
		int num = clusters.size();
		double maxCC = -1;
		int maxOne = -1;
		int maxTwo = -1; 
		for (int i=0; i<num; i++){
			int j = i+1;
			Cluster c1 = (Cluster)clusters.get(i);
			float[] one = c1.getValues();
			for (; j<num; j++){
				Cluster c2 = (Cluster)clusters.get(j);
				//look for previously calculated correlation coefficient
				Object objFor = correlationCoefficients.get(c1.getFloatArrayFile().getName()+":"+c2.getFloatArrayFile().getName());
				double pcc;
				//calculate new cc
				if (objFor == null){
					float[] two = c2.getValues();
					pcc = PearsonCorrelation.correlationCoefficient(one, two);
					//save cc in hash
					correlationCoefficients.put(c1.getFloatArrayFile().getName()+":"+c2.getFloatArrayFile().getName(), new Double(pcc));
				}
				//already exists
				else pcc = ((Double)objFor).doubleValue();

				//set max values?
				if (pcc>maxCC) {
					maxCC = pcc;
					maxOne = i;
					maxTwo = j;
				}
			}
		}
		return new double[] {maxCC, maxOne, maxTwo};
	}


	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		File dir = null;
		boolean convert = false;
		boolean antiLog = false;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'd': dir = new File(args[i+1]); i++; break;
					case 'c': convert = true; break;
					case 'a': antiLog = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group()+"\n");
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (dir == null || dir.canRead() == false) {
			Misc.printExit("\nCannot find or read your float[] directory!\n");
		}
		if (convert){
			File[] toConvert = IO.extractFiles(dir);
			for (int i=0; i< toConvert.length; i++){
				File convertedFile = new File (toConvert[i].getParentFile(), Misc.removeExtension(toConvert[i].getName())+".celp");
				float[] values = Num.loadFloats(toConvert[i],0);
				if (values == null) System.out.println("\nError: something is wrong with this float file "+toConvert[i].getName()+", skipping\n");
				else {
					if (antiLog) values = Num.antiLog(values, 2);
					IO.saveObject(convertedFile, values);
				}
			}
		}
		celpFiles = IO.extractFiles(dir, "celp");
		if (celpFiles.length<2) Misc.printExit("\nOnly one file found?! Too few for clustering.\n");
	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Hierarchical Clustering: June 2012                        **\n" +
				"**************************************************************************************\n" +
				"HC hierarchically clusters serialized float[] arrays using a Pearson correlation\n" +
				"coefficient (r) as a metric.  For each round, all pairwise r values between the arrays\n" +
				"are calculated, the pair of arrays with the highest r value is removed from the pool,\n" +
				"their intensities are averaged, and the averaged array is added back to the pool.\n" +
				"Rounds of clustering continue until only one cluster remains. R values are typically\n" +
				"squared and multiplied by 100 to give a similarity percentage. For ChIP chip\n" +
				"experiments, clusters with an r < 0.75 (50% similar) are processed separately.\n\n"+ 

				"-d Full path directory text containing serialized float[] arrays, 'xxx.celp' files.\n" +
				"-c These are text files containing one column of floats, convert to 'xxx.celp'\n" +
				"-a AntiLog base 2 the values.\n\n"+

				"Example: java -Xmx1500M -jar pathTo/T2/HierarchicalClustering -d /affy/CelpFiles\n\n" +

		"**************************************************************************************\n");		
	}

	public static void main(String[] args) {
		if (args.length==0){
			printDocs();
			System.exit(0);	  
		}
		new HierarchicalClustering(args);
	}

	public ArrayList getClusters() {
		return clusters;
	}

	public Cluster getMegaCluster() {
		return megaCluster;
	}

	public double getMinimalCorrelationCoefficient() {
		return minimalCorrelationCoefficient;
	}

	public boolean savePNG() {
		return savePNG;
	}

	public void setSavePNG(boolean savePNG) {
		this.savePNG = savePNG;
	}

	public File getPngResultFile() {
		return pngResultFile;
	}

	public void setPngResultFile(File pngResultFile) {
		this.pngResultFile = pngResultFile;
	}
	public String getTitle() {
		return title;
	}
	public void setTitle(String panelTitle) {
		this.title = panelTitle;
	}

	public boolean displayClusterPlot() {
		return displayClusterPlot;
	}

	public void setDisplayClusterPlot(boolean displayClusterPlot) {
		this.displayClusterPlot = displayClusterPlot;
	}




}
