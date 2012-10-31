package edu.utah.seq.data;
import util.gen.*;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.parsers.Tag2Point;

public class RemovePointDataScores {

	//fields
	private File saveDirectory;
	private File[] pointDataFiles;

	
	/**Randomly samples out reads to match a desired number.*/
	public RemovePointDataScores(String[] args){
		pointDataFiles = IO.extractFiles(new File(args[0]), ".bar.zip");
		if (pointDataFiles == null || pointDataFiles.length ==0){
			pointDataFiles = IO.extractFiles(new File(args[0]), ".bar");
		}
		saveDirectory = new File(args[1]);
		saveDirectory.mkdir();
		
		//for each point directory
		for (int i=0; i< pointDataFiles.length; i++){
			System.out.println(pointDataFiles[i]);
			//fetch the data
			PointData pd = new PointData(pointDataFiles[i], true);
			float[] scores = pd.getScores();
			Arrays.fill(scores, 1);
			pd.writePointData(saveDirectory);
		}
		System.out.println("Done!");
	}

	public static void main(String[] args) {
		if (args.length ==0){
			System.out.println("\nEnter a PointData directory and a save directory.\n");
			System.exit(0);
		}
		new RemovePointDataScores(args);
	}		

	


}
