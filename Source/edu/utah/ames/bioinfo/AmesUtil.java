package edu.utah.ames.bioinfo;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.zip.GZIPOutputStream;

public class AmesUtil {
	
	public PrintWriter getPrintWriter(File file){
		int fileNumber = 0;
		PrintWriter out =null;
		try {
			String name = file.getName();
			name = name.replaceAll(".gz", "");
			name = name.replaceAll(".zip", "");
			out = new PrintWriter(new FileWriter(new File(file.getParentFile(),fileNumber+"_"+name)));
			fileNumber++;
		}
		catch (IOException e){
			e.printStackTrace();
		}
		return out;
	}

	public GZIPOutputStream getGZipOutputStream(File file){
		int fileNumber = 0;
		GZIPOutputStream out =null;
		try {
			String name = file.getName();
			name = name.replaceAll(".gz", "");
			name = name.replaceAll(".zip", "");
			name = name+".gz";
			out = new GZIPOutputStream(new FileOutputStream(new File(file.getParentFile(),fileNumber+"_"+name)));
			fileNumber++;
		}
		catch (IOException e){
			e.printStackTrace();
		}
		return out;
	}
	
	/**
	 * TODO DEBUG: doesn't work like I thought it would!!
	 * returns filepath stripped of extension 
	 * @param filePath
	 * @return
	 */
	public static String removeExtension(String filePath) {
		File f = new File(filePath);
		
		//if it's a directory, don't remove the extension
		if (f.isDirectory()) return filePath;
		
		String name = f.getAbsolutePath();
		
		//we know it's a file
		final int lastPeriodPos = name.lastIndexOf('.');
		if (lastPeriodPos <= 0) {
			
			//no period after first character, so return name as it was passed in
			return filePath;
		}
		else {
			//remove the last period and everything after it
			File renamed = new File(name.substring(0, lastPeriodPos));
			return renamed.getPath();
		}
	}
	
	/**
	 * Sort by value of hashmap
	 * @param map
	 * @return
	 */
	static <K,V extends Comparable<? super V>> List<Entry<K,V>> entriesSortedByValues(Map<K,V> map) {

		List<Entry<K,V>> sortedEntries = new ArrayList<Entry<K,V>>(map.entrySet());

		Collections.sort(sortedEntries, new Comparator<Entry<K,V>>() {
			@Override
			public int compare(Entry<K,V> e1, Entry<K,V> e2) {
				return e2.getValue().compareTo(e1.getValue());
			}
		}
				);

		return sortedEntries;
	}
	
	/**
	 * HashMap containing the amino acid codon table
	 */
	@SuppressWarnings("serial")
	public static HashMap<String, String> aminoAcidTable = new HashMap<String, String>(){{
		//put in codon translations
		//alanine
		put("GCT","A");
		put("GCC","A");
		put("GCA","A");
		put("GCG","A");
		//arginine
		put("CGT","R");
		put("CGC","R");
		put("CGA","R");
		put("CGG","R");
		put("AGA","R");
		put("AGG","R");
		//asparagine
		put("AAT","N");
		put("AAC","N");
		//aspartic acid
		put("GAT","D");
		put("GAC","D");
		//cysteine
		put("TGT","C");
		put("TGC","C");
		//glutamine
		put("CAA","Q");
		put("CAG","Q");
		//glutamic acid
		put("GAA","E");
		put("GAG","E");
		//glycine
		put("GGT","G");
		put("GGC","G");
		put("GGA","G");
		put("GGG","G");
		//histidine
		put("CAT","H");
		put("CAC","H");
		//isoleucine
		put("ATT","I");
		put("ATC","I");
		put("ATA","I");
		//leucine
		put("TTA","L");
		put("TTG","L");
		put("CTT","L");
		put("CTC","L");
		put("CTA","L");
		put("CTG","L");
		//lysine
		put("AAA","K");
		put("AAG","K");
		//methionine
		put("ATG","M");
		//phenylalanine
		put("TTT","F");
		put("TTC","F");
		//proline
		put("CCT","P");
		put("CCC","P");
		put("CCA","P");
		put("CCG","P");
		//serine
		put("TCT","S");
		put("TCC","S");
		put("TCA","S");
		put("TCG","S");
		put("AGT","S");
		put("AGC","S");
		//threonine
		put("ACT","T");
		put("ACC","T");
		put("ACA","T");
		put("ACG","T");
		//tryptophan
		put("TGG","W");
		//tyrosine
		put("TAT","Y");
		put("TAC","Y");
		//valine
		put("GTT","V");
		put("GTC","V");
		put("GTA","V");
		put("GTG","V");
		//start
		put("ATG","(start)");
		//stop
		put("TAA","(stop-ochre)");
		put("TGA","(stop-opal)");
		put("TAG","(stop-amber)");
	}};
	
	/**
	 * rounds up a double to return an int
	 * @param d
	 * @return
	 */
	private static int roundUp(double d) {
		return (d > (int) d) ? (int) d + 1 : (int) d;
	}
	
	/**
	 * Generates a Map of results (Map of key to boolean). Works nicely regardless of
	 * different keys and key sort order
	 * @param map1
	 * @param map2
	 * @return
	 */
	public static <K extends Comparable<? super K>, V> Map<K, Boolean> 
	compareEntries(final Map<K, V> map1, final Map<K, V> map2) {
		
		final Collection<K> allKeys = new HashSet<K>();
		allKeys.addAll(map1.keySet());
		allKeys.addAll(map2.keySet());
		final Map<K, Boolean> result = new TreeMap<K, Boolean>();
		
		for (final K key : allKeys) {
			result.put(key,  map1.containsKey(key) == map2.containsKey(key) &&
					Boolean.valueOf(equal(map1.get(key), map2.get(key))));
		}
		return result;
	}
	
	/**
	 * Usage below in commented section
	 * @param obj1
	 * @param obj2
	 * @return
	 */
	private static boolean equal(final Object obj1, final Object obj2) {
		return obj1 == obj2 || (obj1 != null && obj1.equals(obj2));
	}
	
	/**
	final Map<Integer, Boolean> comparisonResult =
	        compareEntries(map1, map2);
	    for(final Entry<Integer, Boolean> entry : comparisonResult.entrySet()){
	        System.out.println("Entry:" + entry.getKey() + ", value: "
	            + entry.getValue());
	    } **/
	    
}
