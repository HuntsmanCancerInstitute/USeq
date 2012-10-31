package util.apps;

import java.io.*;
import java.util.*;

/**Counts the frequency of words in the master list, then assigns these frequencies to words in the test list.*/
public class GeneNameCounter {

	public static void main(String[] args) {
		if (args.length==0){
			System.out.println("\n Enter a full path file text for the master list and the test list.\n");
			System.exit(0);
		}
		//load words into an array
		String[] words = loadFile(new File (args[0]));
		
		//run thru list adding words to hash, incrementing counter.
		HashMap keyValue = new HashMap();
		for (int i=0; i< words.length; i++){
			
			if (keyValue.containsKey(words[i])){
				Integer integer = (Integer)keyValue.get(words[i]);
				int value = integer.intValue() + 1;
				integer = new Integer (value);
				keyValue.put(words[i], integer);
			}
			else keyValue.put(words[i], new Integer(1));
		}
		
		//print hash
		System.out.println("Counted Master List: "+keyValue); 
		
		//load test list
		String[] test = loadFile(new File (args[1]));
		System.out.println("\nCounts for Test List:");
		//look for test words in master list, print number of repeats
		for (int i=0; i<test.length; i++){
			if (keyValue.containsKey(test[i])){
				int value = ((Integer)keyValue.get(test[i])).intValue();
				System.out.println(test[i]+"\t"+ value);
			}
			else {
				System.out.println(test[i]+"\t"+ 0); 
			}
		}
		
		
	}
	/**Loads a file's lines into a String[], will save blank lines.*/
	public static String[] loadFile(File file){
		ArrayList a = new ArrayList();
		try{
			BufferedReader in = new BufferedReader(new FileReader(file));
			String line;
			while ((line = in.readLine())!=null){
				line = line.trim();
				a.add(line);
			}
		}catch(Exception e){
			e.printStackTrace();
		}
		String[] strings = new String[a.size()];
		a.toArray(strings);
		return strings;
	}
}
