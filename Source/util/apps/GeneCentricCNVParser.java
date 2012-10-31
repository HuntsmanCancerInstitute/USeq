package util.apps;
import java.io.*;
import util.gen.*;
import java.util.*;
import java.util.regex.*;

public class GeneCentricCNVParser {

	HashMap<String, Scores> scores = new HashMap<String, Scores>();
	static final Pattern tabPattern = Pattern.compile("\\t");
	static final Pattern whiteSpacePattern = Pattern.compile("\\s+");
	boolean printScores = true;
	
	
	public GeneCentricCNVParser(String[] args){
		
		//pull directories
		File norm = new File (args[0]);
		File oligo = new File (args[1]);
		File azoo = new File (args[2]);
		
		//pull parsed files
		File[] psNorm = IO.extractFiles(norm, ".txt");
		File[] psOligo = IO.extractFiles(oligo, ".txt");
		File[] psAzoo = IO.extractFiles(azoo, ".txt");
		
		//for each parsed file load
		loadScores (psNorm, "norm");
		loadScores (psOligo, "oligo");
		loadScores (psAzoo, "azoo");
		
		//set sorted scores
		setScores();
		
		//find max length
		//findMaxLength();
		
		//print scores
		printResults();
		
		
		
	}
	
	public void printResults(){
		Iterator<String> regions = scores.keySet().iterator();
		System.out.println("Name\tMean\t5th\t25th\t50th\t75th\t95th");
		String[] toFetch = {"chr19:1332426-1451407", "chr5:178373455-178714935", "chr15:25929782-26250890", "chr5:70305175-70434776", "chr6:34212627-34331986", "chr12:116035361-116293965", "chr19:5474163-5629489", "chr12:7655355-7771419", "chr2:48667416-48846384", "chr19:422324-544493", "chr1:151017421-151134147", "chr2:213757360-214993470", "chr11:9339088-9516647", "chr1:144024547-144138902", "chr11:76416957-76613934", "chr11:63263470-63443629", "chr8:72172221-72447021", "chr7:143626971-143748253", "chr5:1206286-1358162", "chr17:37877876-37993375", "chr1:37998955-38195316", "chr2:119216216-119332229", "chr11:340279-491387", "chr19:358496-480654"};
		HashSet names = Misc.loadHashSet(toFetch);
		while (regions.hasNext()){
			//text
			String name = regions.next();
			if (names.contains(name)== false) continue;
			System.out.println(name);
			Scores scr = scores.get(name);
			//norm scores
			float[] norm = scr.sortedNormScores;
			if (printScores) System.out.println("norm = c("+Misc.floatArrayToString(norm, ",")+")");
			//Azoo scores
			float[] azoo = scr.sortedAzooScores;
			if (printScores) System.out.println("azoo = c("+Misc.floatArrayToString(azoo, ",")+")");
			//norm scores
			float[] oligo = scr.sortedOligoScores;
			if (printScores) System.out.println("oligo = c("+Misc.floatArrayToString(oligo, ",")+")");
			
			System.out.println("Norm \t"+Num.mean(norm)+"\t"+Num.percentile(norm, 0.05)+"\t"+Num.doubleArrayToString(Num.quartiles(norm), 1, "\t")+"\t"+Num.percentile(norm, 0.95));
			System.out.println("Azoo \t"+Num.mean(azoo)+"\t"+Num.percentile(azoo, 0.05)+"\t"+Num.doubleArrayToString(Num.quartiles(azoo), 1, "\t")+"\t"+Num.percentile(azoo, 0.95));
			System.out.println("Oligo\t"+Num.mean(oligo)+"\t"+Num.percentile(oligo, 0.05)+"\t"+Num.doubleArrayToString(Num.quartiles(oligo), 1, "\t")+"\t"+Num.percentile(oligo, 0.95));
			System.out.println();
			//System.exit(0);
			
			
			
		}
	}
	

	
	public void setScores(){
		Iterator<String> region = scores.keySet().iterator();
		while (region.hasNext()){
			Scores scr = scores.get(region.next());
			//norm
			float[] norm = Num.arrayListOfFloatToArray(scr.norm);
			Arrays.sort(norm);
			scr.sortedNormScores = norm;
			//oligo
			float[] oligo = Num.arrayListOfFloatToArray(scr.oligo);
			Arrays.sort(oligo);
			scr.sortedOligoScores = oligo;
			//azoo
			float[] azoo = Num.arrayListOfFloatToArray(scr.azoo);
			Arrays.sort(azoo);
			scr.sortedAzooScores = azoo;
		}
	}
	
	public void loadScores(File[] files, String dataset){
		for (int i=0; i<files.length; i++){
			String[] lines = IO.loadFileIntoStringArray(files[i]);
			for (int j=0; j< lines.length; j++){
				String[] keyVal = tabPattern.split(lines[j]);
				Scores s;
				//is in Hash?
				if (scores.containsKey(keyVal[0])){
					s = scores.get(keyVal[0]);
				}
				else {
					s = new Scores();
					scores.put(keyVal[0], s);
				}
				//add scores
				String[] split = whiteSpacePattern.split(keyVal[1]);
				ArrayList<Float> al;
				if ("oligo".equals(dataset)) al = s.oligo;
				else if ("azoo".equals(dataset)) al = s.azoo;
				else al = s.norm;
				for (int k=0; k< split.length; k++) al.add(new Float(split[k]));
			}
		}
	}

	public static void main(String[] args) {
		new GeneCentricCNVParser (args);

	}
	
	private class Scores{
		ArrayList<Float> norm = new ArrayList<Float>();
		ArrayList<Float> oligo = new ArrayList<Float>();
		ArrayList<Float> azoo = new ArrayList<Float>();
		float[] sortedNormScores;
		float[] sortedOligoScores;
		float[] sortedAzooScores;
	}
}
