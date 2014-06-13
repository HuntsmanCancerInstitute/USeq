package edu.utah.ames.bioinfo;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import util.gen.Num;

public class MathStuff {

	private final String filename1 = "/Users/darren/Desktop/proteaseInhibitorOptimization/BC17.optimizer.txt";
	private final String filename2 = "/Users/darren/Desktop/proteaseInhibitorOptimization/BC18.optimizer.txt";
	
	public static void main(String[] args) {

		
	}
	
	public HashMap<String,String> loadMaps() throws IOException{

		HashMap<String,String> hm = new HashMap<String,String>();
		
		//make new hashmap
		try{
			//read in the file of index codes and names
			BufferedReader br1 = new BufferedReader(new FileReader(filename1));
			String line;
			
			//skip first two lines
			br1.readLine();
			br1.readLine();
			
			//assign values
			String[] keyValue;
			//split file on tabs
			while ((line = br1.readLine()) != null){
				keyValue = line.split("\t");
				//put the key/value pairs in the hashmap
				hm.put(keyValue[0].trim(), keyValue[1].trim());
			}
		} 
		catch(Exception e) {
			System.out.println("Problem loading the Hash");
			e.printStackTrace();
		}
		return hm;
	}
	
}
