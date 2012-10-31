package trans.roc;
import java.io.*;
import java.util.*;

/**Parses out non overlapping {@link RocWindow}s.*/
public class ExtractNonOverlappingWindows {
	//constructor
	public ExtractNonOverlappingWindows(File txtRocFile){
		RocWindow[] rocWindows = ParsePatternedWindows.loadTxtRocFile(txtRocFile);
		Arrays.sort(rocWindows);
		System.out.println("\nNum Roc Windows: "+rocWindows.length);
		RocWindow[] nonOverlapping = extractNonOverlappingWindows(rocWindows);
		System.out.println("\nNum Non Overlapping Roc Windows: "+nonOverlapping.length);
		ParsePatternedWindows.writeTxtRocFile(nonOverlapping, new File (txtRocFile+"NonOL"));
	}
	
	/**Assumes the windows are sorted*/
	public static RocWindow[] extractNonOverlappingWindows(RocWindow[] r){
		ArrayList uni = new ArrayList();
		//walk through left to right
		RocWindow first = r[0];
		for (int i=1; i<r.length; i++){
			//do they intersect
			if (first.intersects(r[i])){
				//is the score better
				if (r[i].getScore()>= first.getScore()){
					//replace first
					first = r[i];
				}
				//add if last
				if ( (i+1) == r.length) uni.add(first);
			}
			//no they don't intersect
			else{
				uni.add(first);
				first = r[i];
			}
		}
		RocWindow[] newWin = new RocWindow[uni.size()];
		uni.toArray(newWin);
		return newWin;
	}
	
	public static void main(String[] args){
		if (args.length ==0) {
			System.out.println("\nEnter a full path file text for a txt roc file (tab delimited: chrom, start, stop, middle(not used), score)\n");
			System.exit(0);
		}
		new ExtractNonOverlappingWindows(new File(args[0]));
	}
	
}
