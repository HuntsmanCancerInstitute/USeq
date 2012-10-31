package edu.utah.seq.data;

import java.io.*;
import java.util.*;
import util.gen.*;

/**Merges Point Data and then divides using random sampling.*/
public class DividePointData {


	public static void main(String[] args) {
		if (args.length ==0) Misc.printExit("\nUsage: comma delimited list of PointData directories to merge and divide, the number of divisions, where to save them, numberBPToShift3Prime\n");
		//fetch and merge data
		File[] pointDirs = IO.extractFiles(args[0]);
		//only one directory look deeper
		if (pointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(pointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) pointDirs = otherDirs;
		}
		System.out.println("Merging...");
		HashMap<String, PointData> data = PointData.mergePointData(PointData.fetchPointData(pointDirs,null, false), false, true);
		System.out.println("Total points "+PointData.totalObservations(data));
		
		//split it
		HashMap<String, PointData>[] split = PointData.divide(data, Integer.parseInt(args[1]));
		
		//save em
		File saveDir = new File (args[2]);
		saveDir.mkdir();
		
		//splits
		for (int i=0; i< split.length; i++){
			File sf = new File(saveDir,"Split"+i);
			sf.mkdir();
			System.out.println("Saving split "+i+" "+PointData.totalObservations(split[i]));
			PointData.writePointData(split[i], sf);
		}
		

	}


}
