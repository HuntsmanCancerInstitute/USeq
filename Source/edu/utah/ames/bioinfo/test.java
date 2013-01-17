package edu.utah.ames.bioinfo;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;

public class test {
	
	private static String x = "/Users/darren/Desktop/novoalignerTestDir/novoindexNames.txt";

	public static HashMap<String,String> loadNovoindex(){
		HashMap<String,String> indices = new HashMap<String,String>(10000);
		try{
			BufferedReader br = new BufferedReader(new FileReader(x));
			String line;
			String[] keyValue;
			while ((line = br.readLine())!=null){
				keyValue = line.split("\t");
				indices.put(keyValue[0].trim(), keyValue[1].trim());
			}
		}catch(Exception e){
			System.out.println("Problem loading the Hash");
			e.printStackTrace();
		}
		System.out.println(indices.entrySet());
		return indices;
	}
	
	public static void main(String[] args) {
		loadNovoindex();
	}
}
