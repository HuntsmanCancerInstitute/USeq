package edu.utah.seq.barcodes;


import htsjdk.samtools.SAMRecord;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class BarcodeClusterEngine {

	//fields
	private int minBaseQuality = 20;
	protected double minNumBases = 6;
	protected double minFractionIdentity = 0.75; 
	private HashMap<String, ArrayList<SAMRecord>> nBarSams = new HashMap<String, ArrayList<SAMRecord>>();
	private ArrayList<SAMRecord> mergedRecords = new ArrayList<SAMRecord>();
	private static final Pattern SEQ_QUAL = Pattern.compile(":BMF:([\\x21-\\x7E]+)");
	
	//results
	private ArrayList<SAMRecord> skippedRecords = new ArrayList<SAMRecord>();
	private SAMRecord[][] clusteredRecords = null;

	
	
	//constructors
	public BarcodeClusterEngine(int minBaseQuality, double minNumBases, double minFractionIdentity) {
		this.minBaseQuality = minBaseQuality;
		this.minNumBases = minNumBases;
		this.minFractionIdentity = minFractionIdentity;
	}
	
	//methods
	/**After instantiating a ClusterBarcodedSams object, call this method to cluster the records.  
	 * It'll return a SAMRecord[cluster number][SAMRecords in a given cluster].
	 * Records with barcodes failing the minNumBases with minBaseQuality are dumped into the skippedRecords obj and ignored.*/
	public SAMRecord[][] cluster(ArrayList<SAMRecord> records){
		
		//load the HashMap<String, ArrayList<SAMRecord> of N'ed barcodes, this will collapse based on exact identity
		loadHash(records);

		//all collapsed? some might have been skipped...
		if (nBarSams.size() == 1){
			clusteredRecords = new SAMRecord[1][];
			ArrayList<SAMRecord> al = nBarSams.values().iterator().next();
			clusteredRecords[0] = new SAMRecord[al.size()];
			al.toArray(clusteredRecords[0]);
		}
		//nope, then take unique N'ed barcodes and cluster the hash
		else clusterHash();

		return clusteredRecords;
	}

	/*Called when there are actually barcodes needing fuzzy clustering*/
	private void clusterHash() {
		
		//make a BarcodeCluster for each hash key
		ArrayList<BarcodeCluster> clusters = makeBarcodeClusters();
		
		//cluster the single entry clusters, this self collapses
		clusterClusters(clusters);
		
		//make final clusters
		makeFinalClusters(clusters);
	}

	/*Merges the fuzzy clustering results with those with identical barcodes.*/
	private void makeFinalClusters(ArrayList<BarcodeCluster> clusters) {
		//System.out.println("\nFinal clusters:");
		clusteredRecords = new SAMRecord[clusters.size()][];
		//for each merged cluster
		for (int i=0; i< clusteredRecords.length; i++){
			mergedRecords.clear();
			BarcodeCluster c = clusters.get(i);
			//System.out.println(i+"\t"+c);

			//for each barcode in merged cluster, fetch associated SAMRecords from the starting hash and combine
			Iterator<String> it = c.getFamilyMembers().keySet().iterator();
			while (it.hasNext()) mergedRecords.addAll(nBarSams.get(it.next()));
			
			//make array
			clusteredRecords[i] = new SAMRecord[mergedRecords.size()];
			mergedRecords.toArray(clusteredRecords[i]);
		}
	}

	/*Primary method for performing the fuzzy clustering on N'ed barcodes.*/
	private void clusterClusters(ArrayList<BarcodeCluster> clusters) {
		//for each cluster
		int num = clusters.size();
		boolean clustersJoined = true;

		while (clustersJoined){
			//this boolean will determine when all clusters have been made
			clustersJoined = false;
			//for each cluster
			for (int i=0; i< num; i++){
				BarcodeCluster first = clusters.get(i);
				if (first == null) continue;
				//System.out.println("\nNewFirst "+first);
				//look at all subsequent clusters
				for (int j=i+1; j< num; j++){
					BarcodeCluster second = clusters.get(j);
					if (second == null) continue;
					//System.out.println(first+" Comp "+second);
					//are they members?
					if (second.isMemberOf(first)){
						//System.out.println(first+" Joining "+second);
						//yes, add seconds members to first and null second
						first.getFamilyMembers().putAll(second.getFamilyMembers());
						clusters.set(j, null);
						clustersJoined = true;
					}
				}
			}
			//remove nulls? this can be expensive
			if (clustersJoined){
				for (int i=0; i< clusters.size(); i++) {
					if (clusters.get(i) == null) {
						clusters.remove(i);
						i--;
					}
				}
				num = clusters.size();
			}
		}
	}

	private ArrayList<BarcodeCluster> makeBarcodeClusters() {
		ArrayList<BarcodeCluster> clusters = new ArrayList<BarcodeCluster>();
		Iterator<String> it = nBarSams.keySet().iterator();
		while (it.hasNext()) clusters.add(new BarcodeCluster(it.next(), this));
		return clusters;
	}

	private void loadHash(ArrayList<SAMRecord> records) {
		//clear old data
		nBarSams.clear();
		skippedRecords.clear();
		//for each sam
		for (SAMRecord s : records){
			String nBarcode = extractAndNBarcode(s);
			if (nBarcode == null) skippedRecords.add(s);
			else {
				ArrayList<SAMRecord> al = nBarSams.get(nBarcode);
				if (al == null){
					al = new ArrayList<SAMRecord>();
					nBarSams.put(nBarcode, al);
				}
				al.add(s);
			}
		}
	}

	private String extractAndNBarcode(SAMRecord samRecord) {
		Matcher mat = SEQ_QUAL.matcher(samRecord.getReadName());
		if (mat.find()) {
			String seqQual = mat.group(1);
			int halfLen = seqQual.length()/2;
			char[] seq = seqQual.substring(0, halfLen).toCharArray();
			int[] qualScores = BarcodeClusterEngine.convertSangerQualityScores(seqQual.substring(halfLen));
			int numGoodBases = 0;
			for (int i=0; i< qualScores.length; i++){
				if (qualScores[i]< minBaseQuality) seq[i] = 'N';
				else numGoodBases++;
			}
			if (numGoodBases < minNumBases) return null;
			return new String(seq);
		}
		//should never hit this
		System.err.println("ERROR: could not extract a BMF barcode from the following alignment:\n"+samRecord.getSAMString().trim());
		System.exit(1);
		return null;
	}
	
	/**Summarizes the results.*/
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append("Clustered Records:\n");
		for (int i=0; i< clusteredRecords.length; i++){
			sb.append("Cluster "+i+"\n");
			for (int j=0; j< clusteredRecords[i].length; j++){
				sb.append("\t"+clusteredRecords[i][j].getReadName()+"\n");
			}
		}
		sb.append("\nRecords skipped for having bad barcodes: "+skippedRecords.size()+"\n");
		if (skippedRecords.size() !=0){
			for (SAMRecord s : skippedRecords) {
				sb.append(s.getSAMString()+" ");
			}
			sb.append("\n");
		}
		return sb.toString();
	}
	
	//copied over from the USeq Seq utility class
	/**Returns a map of asci text character to it's associated Sanger fastq base quality score.
	 * ! = 0*/
	public static HashMap<String,Integer> asci2FastQScore(){	
		String[] acii = new String[]{"!", "\"", "#", "$", "%", "&", "'", "(", ")", "*", "+", ",", "-", ".", "/", "0", "1", "2", "3", "4", "5", 
			"6", "7", "8", "9", ":", ";", "<", "=", ">", "?", "@", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", 
			"O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "[", "\\", "]", "^", "_", "`", "a", "b", "c", "d", "e", "f", "g", 
			"h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "{", "|", "}", "~"};	
		HashMap<String,Integer> map = new HashMap<String,Integer>();
		for (int i=0; i< acii.length; i++) map.put(acii[i], new Integer(i));
		return map;
	}
	public static final HashMap<String,Integer> asci2FastQScore = asci2FastQScore();
	/**Converts the ascii quality scores to numeric scores. Sanger fastq and illumina 1.8+
	 * See http://onetipperday.blogspot.com/2012/10/code-snip-to-decide-phred-encoding-of.html */
	public static int[] convertSangerQualityScores(String seqQual){
		int[] scores = new int[seqQual.length()];
		for (int i=0; i< seqQual.length(); i++){
			String sub = seqQual.substring(i, i+1);
			Integer val = asci2FastQScore.get(sub);
			if (val == null) {
				System.err.println("\nError converting seq quality character -> "+sub+" from "+seqQual);
				return null;
			}
			scores[i] = val.intValue();
		}
		return scores;
	}

	public ArrayList<SAMRecord> getSkippedRecords() {
		return skippedRecords;
	}

	public SAMRecord[][] getClusteredRecords() {
		return clusteredRecords;
	}
}


