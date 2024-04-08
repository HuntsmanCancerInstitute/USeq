package util.gen;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.apache.log4j.Logger;


/**
 * A variety of static methods.
 *
 */
public class Misc {
	
	public static final Pattern PATTERN_SEMICOLON = Pattern.compile(";");
	public static final Pattern PATTERN_EQUALS = Pattern.compile("=");
	public static final Pattern FORWARD_SLASH = Pattern.compile("/");
	public static final Pattern TAB = Pattern.compile("\t");
	public static final Pattern TAB_MULTIPLE = Pattern.compile("\t+");
	public static final Pattern RETURN = Pattern.compile("\n");
	public static final Pattern WHITESPACE = Pattern.compile("\\s+");
	public static final Pattern COMMA_WHITESPACE = Pattern.compile("[,\\s]+");
	public static final Pattern COMMA_WHITESPACE_FWDSLASH = Pattern.compile("[,\\s/]+");
	public static final Pattern UNDERSCORE = Pattern.compile("_");
	public static final Pattern DASH = Pattern.compile("-");
	public static final Pattern COLON = Pattern.compile(":");
	public static final Pattern COMMA = Pattern.compile(",");
	public static final Pattern PERIOD = Pattern.compile("\\.");
	public static final Pattern SINGLE_QUOTE = Pattern.compile("'");
	public static final Pattern DOUBLE_QUOTE = Pattern.compile("\"");
	public static final Pattern SEMI_COLON = Pattern.compile(";");
	public static final Pattern SEMI_COLON_SPACE = Pattern.compile("\\s*;\\s*");
	public static final Pattern PIPE = Pattern.compile("\\|");
	public static final Pattern FORWARD_PARENTHESIS = Pattern.compile("\\(");
	public static final Pattern GREATER_THAN = Pattern.compile(">");
	public static final Pattern AMP = Pattern.compile("&");
	public static final Pattern CSV = Pattern.compile("\",\"");
	
	/**For identifying MolecularDataPatientIds, see edu.utah.hci.bioinfo.smm.CoreId */
	public static final Pattern CORE_ID_Pattern = Pattern.compile("[A-Za-z&&[^OoIilL]]{2,3}[2-9][A-Za-z&&[^OoIilL]]{2}[2-9][A-Za-z&&[^OoIilL]]{2,3}");

	/**Attempts to parse a MolecularDataPatientID from the name separated by an underscore or from the path where the MDPID is used as a dir in the path
	 * Returns null if not found or throws an error if a problem is an issue.*/
	public static String fetchMDPID(File f) {
		try {
			//attempt to get it from the first part of the name, e.g. from hG9fu8Jx_TN22-135333_2022-05-20_10_47_deid_Neeraj_Agarwal
			String[] split = UNDERSCORE.split(f.getName());
			if (CORE_ID_Pattern.matcher(split[0]).matches()) return split[0];
			
			//attempt to get it from a sub dir in the path, e.g. /scratch/general/pe-nfs1/u0028003/Tempus/Runs/11Aug2022/TJobs/Rk7Ps4er/Tempus/TL-22-VUTSJR3Y/ClinicalReport/TL-22-VUTSJR3Y_XT.V4_2022-07-12_deid_Joe_Mendez_M.json
			String cp = f.getCanonicalPath();
			split = FORWARD_SLASH.split(cp);
			String id = null;
			for (String testId: split) {
				if (CORE_ID_Pattern.matcher(testId).matches()) {
					//not set?
					if (id == null) id = testId;
					//not same?
					else if (id.equals(testId) == false) throw new IOException("ERROR, two different coreIds found: "+id+" "+testId);
				}
			}
			return id;
		} catch (IOException e) {
			System.err.println("ERROR fetching coreId from "+f);
			e.printStackTrace();
			return null;
		}
	}
	
	/**Splits a comma and quotation split csv line where cells are delimited by "," with a " on either end, e.g. "xxxxxxx","yyyyyyyy"
	 * */
	public static String[] splitCsvLine(String line) {
		String[] split = CSV.split(line);
		split[0] = split[0].substring(1);
		split[split.length-1] = split[split.length-1].substring(0, split[split.length-1].length()-1);
		return split;
	}
	
	/**Calls garbage collection then returns total - free.*/
	public static long fetchUsedMemory(){
	    System.gc(); System.gc(); System.gc(); System.gc();
	    System.gc(); System.gc(); System.gc(); System.gc();
	    System.gc(); System.gc(); System.gc(); System.gc();
	    System.gc(); System.gc(); System.gc(); System.gc();
	    return Runtime.getRuntime().totalMemory() -
	      Runtime.getRuntime().freeMemory();
	}
	
	/**If present removes a leading or trailing /   /results/go/ -> results/go */
	public static String removeLeadingTrailingForwardSlashes(String toTrim){
		int beginningIndex = 0;
		int endingIndex = toTrim.length();
		if (toTrim.startsWith("/")) beginningIndex++;
		if (toTrim.endsWith("/")) endingIndex--;
		return toTrim.substring(beginningIndex, endingIndex);
	}
	
	/**Finds Strings common to all, any zero sized or null members will cause a return of a zero sized AL.*/
	public static ArrayList<String> commonToAll (ArrayList<String>[] als){
		//find smallest and copy it
		ArrayList<String> common = new ArrayList<String>();
		if (als[0] == null) return common;
		int smallestIndex = 0;
		int smallestNum = als[0].size();
		for (int i=1; i< als.length; i++){
			if (als[i] == null) return common;
			int size = als[i].size();
			if (size < smallestNum){
				smallestNum = size;
				smallestIndex = i;
			}
		}
		//any zeros?
		if (smallestNum == 0) return common;
		
		//copy
		common.addAll(als[smallestIndex]);
		
		for (int i=0; i< als.length; i++){
			if (i != smallestIndex)  common.retainAll(als[i]);
		}
		return common;
	}
	
	/**Flips key:value to value:key.*/
	public static HashMap invert(HashMap h){
		HashMap inverted = new HashMap (h.size());
		Iterator it = h.keySet().iterator();
		while (it.hasNext()){
			Object o = it.next();
			inverted.put(h.get(o), o);
		}
		return inverted;
	}
	
	/**Reverses a String using a StringBuilder.*/
	public static String reverse(String x){
		return new StringBuilder(x).reverse().toString();
	}
	
	public static void printMemoryUsed(){
		System.gc();
		System.gc();
		System.gc();
		double total = Runtime.getRuntime().totalMemory();
		double used  = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
		String percent = Num.formatPercentOneFraction(used/total);
		System.out.println(percent+" Memory Used");
	}
	
	
	public static void randomize (boolean[] array, long seed){
	    Random rng = new Random(seed);       
	    // n is the number of items left to shuffle
	    for (int n = array.length; n > 1; n--) {
	        // Pick a random element to move to the end
	        int k = rng.nextInt(n);  // 0 <= k <= n - 1.
	        // Simple swap of variables
	        boolean tmp = array[k];
	        array[k] = array[n - 1];
	        array[n - 1] = tmp;
	    }
	}
	
	public static void randomize (Object[] array, long seed){
	    Random rng = new Random(seed);       
	    // n is the number of items left to shuffle
	    for (int n = array.length; n > 1; n--) {
	        // Pick a random element to move to the end
	        int k = rng.nextInt(n);  // 0 <= k <= n - 1.
	        // Simple swap of variables
	        Object tmp = array[k];
	        array[k] = array[n - 1];
	        array[n - 1] = tmp;
	    }
	}
	
	/**For multiple iterations.*/
	public static void randomize (Object[] array, Random rng){     
	    // n is the number of items left to shuffle
	    for (int n = array.length; n > 1; n--) {
	        // Pick a random element to move to the end
	        int k = rng.nextInt(n);  // 0 <= k <= n - 1.
	        // Simple swap of variables
	        Object tmp = array[k];
	        array[k] = array[n - 1];
	        array[n - 1] = tmp;
	    }
	}

	/**For multiple iterations.*/
	public static void randomize (int[] array, Random rng){     
	    // n is the number of items left to shuffle
	    for (int n = array.length; n > 1; n--) {
	        // Pick a random element to move to the end
	        int k = rng.nextInt(n);  // 0 <= k <= n - 1.
	        // Simple swap of variables
	        int tmp = array[k];
	        array[k] = array[n - 1];
	        array[n - 1] = tmp;
	    }
	}
	
	/**Does not check if indexes are out of bounds!
	 * @return the number of trues between the start and stop
	 * @param startIndex is included
	 * @param stopIndex is excluded
	 * */
	public static int countTrues(boolean[] b, int startIndex, int stopIndex){
		int trueCount = 0;
		for (int i=startIndex; i< stopIndex; i++) if (b[i]) trueCount++;
		return trueCount;
	}
	
	/**Fetch random string consisting of ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890*/
	private static String SALTCHARS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890";
	public static String getRandomString(int length) {
        StringBuilder salt = new StringBuilder();
        Random rnd = new Random();
        while (salt.length() < length) { // length of the random string.
            int index = (int) (rnd.nextFloat() * SALTCHARS.length());
            salt.append(SALTCHARS.charAt(index));
        }
        String saltStr = salt.toString();
        return saltStr;
    }
	
	/**Prints key + seperator + value+ return to out.*/
	public static void printHashMap(HashMap hash, String seperator){
		Iterator it = hash.keySet().iterator();
		while (it.hasNext()){
			String key = (String) it.next();
			System.out.println(key+seperator+hash.get(key));
		}
	}
	
	/**Prints key + seperator + value+ return to out.*/
	public static void printTreeMap(TreeMap hash, String seperator, String prepender){
		Iterator it = hash.keySet().iterator();
		while (it.hasNext()){
			String key = (String) it.next();
			System.out.println(prepender+key+seperator+hash.get(key));
		}
	}
	
	public static void printCollection(Collection hash){
		Iterator it = hash.iterator();
		while (it.hasNext()){
			System.out.println(it.next().toString());
		}
	}
	
	/**Converts a hash to a String[].*/
	public static String[] hashSetToStringArray(HashSet<String> hash){
		Iterator<String> it = hash.iterator();
		String[] s = new String[hash.size()];
		int counter =0;
		while (it.hasNext()){
			s[counter++] = it.next();
		}
		return s;
	}
	
	/**Converts a linked hash to a String[]. Not sure if this is necessary?*/
	public static String[] linkedHashSetToStringArray(LinkedHashSet<String> hash){
		Iterator<String> it = hash.iterator();
		String[] s = new String[hash.size()];
		int counter =0;
		while (it.hasNext()){
			s[counter++] = it.next();
		}
		return s;
	}
	
	/**Converts a hash to a String using the delimiter to join.*/
	public static String hashSetToString(HashSet hash, String delimiter){
		if (hash.size() == 0) return "";
		Iterator it = hash.iterator();
		StringBuilder s = new StringBuilder();
		Object obj = (Object)it.next();
		s.append(obj.toString());
		while (it.hasNext()){
			obj = it.next();
			s.append(delimiter);
			s.append(obj.toString());
		}
		return s.toString();
	}
	
	/**Converts a hash to a String using the delimiter to join.*/
	public static String treeSetToString(TreeSet hash, String delimiter){
		if (hash.size() == 0) return "";
		Iterator it = hash.iterator();
		StringBuilder s = new StringBuilder();
		Object obj = (Object)it.next();
		s.append(obj.toString());
		while (it.hasNext()){
			obj = it.next();
			s.append(delimiter);
			s.append(obj.toString());
		}
		return s.toString();
	}
	
	/**Converts a hash to a String using the delimiter to join.*/
	public static String linkedSetToString(LinkedHashSet hash, String delimiter){
		if (hash.size() == 0) return "";
		Iterator it = hash.iterator();
		StringBuilder s = new StringBuilder();
		Object obj = (Object)it.next();
		s.append(obj.toString());
		while (it.hasNext()){
			obj = it.next();
			s.append(delimiter);
			s.append(obj.toString());
		}
		return s.toString();
	}
	
	
	/**Converts a hash to a String[].*/
	public static String[] setToStringArray(Set<String> hash){
		Iterator<String> it = hash.iterator();
		String[] s = new String[hash.size()];
		int counter =0;
		while (it.hasNext()){
			s[counter++] = it.next();
		}
		return s;
	}
	
	/**Converts a hashMap to a String key1=value1;key2=value2...*/
	public static String hashMapToString(HashMap<String, String> hash){
		StringBuilder sb = new StringBuilder();
		Iterator<String> it = hash.keySet().iterator();
		String key = it.next();
		String value = hash.get(key);
		sb.append(key);
		sb.append("=");
		sb.append(value);
		while (it.hasNext()){
			key = it.next();
			value = hash.get(key);
			sb.append(";");
			sb.append(key);
			sb.append("=");
			sb.append(value);
		}
		return sb.toString();
	}
	
	/**Sleeps for given number of seconds*/
	public static void sleep(int seconds){
		try{
			java.lang.Thread.sleep( seconds * 1000 );
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	/**Looks for key, if not found prints error message and exists.*/
	public static String fetchExit(String key, HashMap map){
		if (map.containsKey(key) == false) Misc.printExit("\nError: the key is not present in the HashMap -> "+ key);
		return (String) map.get(key);
	}
	
	/**Returns a new array with the newValue appended on the beginning.*/
	public static String[] prependString (String[] d, String newValue){
		String[] added = new String[d.length+1];
		System.arraycopy(d,0,added,1,d.length);
		added[0] = newValue;
		return added;
	}
	
	/**Concatinates the toRepeat together with a StringBuilder.*/
	public static String concatinateRepeats(String toRepeat, int numberOfRepeats){
		StringBuilder sb = new StringBuilder();
		for (int i=0; i< numberOfRepeats; i++) sb.append(toRepeat);
		return sb.toString();
	}
	
	/**Inserts the String into the char[].*/
	public static char[] insertString (char[] seq, String toInsert, int indexToInsert){
		int altSize = toInsert.length();
		int seqSize = seq.length;
		int nSize = seqSize + altSize;
		char[] nSeq = new char[nSize];
		//Arrays.fill(nSeq, 'x');
		//copy left side
		System.arraycopy(seq, 0, nSeq, 0, indexToInsert);
		//add middle
		int counter = indexToInsert;
		for (int i=0; i < altSize; i++) nSeq[counter++] = toInsert.charAt(i);
		//copy right side
		System.arraycopy(seq, indexToInsert, nSeq, indexToInsert+altSize, seqSize-indexToInsert);
		return nSeq;
	}
	
	/**Capitalizes the first letter in a String.*/
	public static String capitalizeFirstLetter(String s){
		char[] first = s.toCharArray();		
		if (Character.isLetter(first[0])){
			first[0] = Character.toUpperCase(first[0]);
			return new String(first);
		}
		return s;
	}
	
	/**Removes .gz, .zip. and then an extension if found xxx.txt.gz -> xxx
	 * If none found returns the original.
	 */
	public static String removeExtension(String txt) {
		txt = txt.replaceAll("\\.gz", "");
		txt = txt.replaceAll("\\.zip", "");
		int index = txt.lastIndexOf(".");
		if (index != -1)  return txt.substring(0,index);
		return txt;
	}
	
	/**Given a String of 'treat1=/file1,/file2*treat2=/file3,/file4*treat3=/file5,file6'
	 * Will return an ArrayList contining first a String[] of the keys {treat1,treat2,treat3}
	 * and second a String[][] of the parsed values matching the keys.
	 * Returns null if an error is found.*/
	public static ArrayList parseKeyValueGroups(String line){		
		String[] groups = line.split("\\*");
		String[] names = new String[groups.length];
		String[][] items = new String[groups.length][];		
		for (int x=0; x<groups.length; x++){
			String[] keyValue = groups[x].split("=");
			if (keyValue.length !=2 )return null;
			names[x] = keyValue[0];
			items[x] = keyValue[1].split(",");
		}
		ArrayList al = new ArrayList(2);
		al.add(names);
		al.add(items);
		return al;
	}
	
	/**Trims characters common to all from start of lines.*/
	public static String[] trimCommonStart(String[] lines){
		//trim front
		boolean go = true;
		int clip = 0;
		while (go){
			//check length
			if (lines[0].length() <= clip){
				clip--;
				break;
			}
			char test = lines[0].charAt(clip);
			//for each line
			for (int i=1; i< lines.length; i++){
				//check if long enough
				if (lines[i].length() <= clip) {
					clip--;
					go = false;
					break;
				}
				
				else if (test != lines[i].charAt(clip) ){
					go = false;
					break;
				}
			}
			if (go) clip++;
		}
		String[] clipped = new String[lines.length];
		for (int i=0; i< lines.length; i++){
			clipped[i] = lines[i].substring(clip);
		}
		return clipped;
	}
	
	/**Trims the comon chars from the front and back of each line, won't entirely delete line.*/
	public static String[] trimCommon(String[] lines){
		//this is sort of lazy Nix! :)
		//trim front
		String[] trimmed = trimCommonStart(lines);
		//reverse lines
		for (int i=0; i< lines.length; i++){
			trimmed[i] = new StringBuffer(trimmed[i]).reverse().toString();
		}
		//trim back
		trimmed = trimCommonStart(trimmed);
		//reverse lines
		for (int i=0; i< lines.length; i++){
			trimmed[i] = new StringBuffer(trimmed[i]).reverse().toString();
		}
		return trimmed;
	}
	
	/**Trims the comon chars from the back of each line, won't entirely delete line.*/
	public static String[] trimCommonEnd(String[] lines){
		//this is sort of lazy Nix! :)
		//reverse
		String[] trimmed = new String[lines.length];
		for (int i=0; i< lines.length; i++){
			trimmed[i] = new StringBuffer(lines[i]).reverse().toString();
		}
		//trim start of reversed lines
		trimmed = trimCommonStart(trimmed);
		//reverse lines
		for (int i=0; i< lines.length; i++){
			trimmed[i] = new StringBuffer(trimmed[i]).reverse().toString();
		}
		return trimmed;
	}
	

	
	/**Given a String of 'treat1=/file1,treat2=/file3, etc', no spaces,
	 * returns a LinkedHashMap of keys=values.
	 * Returns null if an error is found.*/
	public static LinkedHashMap parseKeyValues(String line){
		String[] keyValues = line.split(",");
		if (keyValues.length==0)return null;
		LinkedHashMap hash = new LinkedHashMap(keyValues.length);
		for (int i=0; i<keyValues.length; i++){
			String[] kv = keyValues[i].split("=");
			if (kv.length!=2)return null;			
			hash.put(kv[0],kv[1]);
		}
		return hash;
	}
	
	/**Breaks a key1=value1;key2=value2 String by semicolon and then by equals sign and loads the data into a hash map.
	 * @return null if something bad happened as well as prints error message to System.err*/
	public static HashMap<String,String> parseKeyValuesToHashMap(String keyValues){
		HashMap<String,String> map = new HashMap<String,String>();
		String[] kvs = PATTERN_SEMICOLON.split(keyValues);
		for (int i=0; i< kvs.length; i++){
			String[] kv = PATTERN_EQUALS.split(kvs[i]);
			if (kv.length != 2) {
				System.err.println("Problem parsing key values, could not split on '=' sign?! from "+keyValues);
				return null;
			}
			map.put(kv[0].trim(), kv[1].trim());
		}
		return map;
	}
	
	/**Attempts to match the stop of the word with the match, case insensitive, replaces it if found,
	 * otherwise returns null.*/
	public static String replaceEnd(String word, String match, String replacement){
		Pattern p = Pattern.compile(match +"$", Pattern.CASE_INSENSITIVE);
		Matcher m = p.matcher(word);
		if (m.find()){
			return m.replaceFirst(replacement);
		}
		return null;
	}
	
	/**Prints message to screen, then exits.*/
	public static void printErrAndExit (String message){
		System.err.println (message);
		System.exit(1);
	}
	/**Prints message to screen, then exits.*/
	public static void printExit (String message){
		System.out.println (message);
		System.exit(0);
	}
	/**Counts the total length of an Object[][]*/
	public static int totalLength(Object[][] array){
		int length = 0;
		int num = array.length;
		for (int i=0; i<num; i++){
			length += array[i].length;
		}
		return length;
	}
	
	/**Removes .txt from stop of String if present*/
	public static String trimTxt(String s){
		if (s.endsWith(".txt")) return s.substring(0,s.length()-4);
		return s; 
	}
	
	/**Prepends zeros to a text such that its total length = 10.  Use for converting a String number to a form that can be 
	 * correctly sorted.  ie would convert 39 to 0000000039*/
	public static String addZeros(String x){
		int numDigits = 10;
		String[] zeros = {"", "0", "00", "000", "0000", "00000", "000000", "0000000", "00000000", "000000000", "0000000000"};
		int length = numDigits - x.length();
		if (length !=0 ) x = zeros[length]+x;
		return x;
	}

	/**Prepends zeros to each String such that its total length = 10.  Use for converting a String number to a form that can be 
	 * correctly sorted.  ie would convert 39 to 0000000039*/
	public static String[] addZeros(String[] x){
		int num = x.length;
		int numDigits = 10;
		String[] zeros = {"", "0", "00", "000", "0000", "00000", "000000", "0000000", "00000000", "000000000", "0000000000"};
		for (int i=0; i<num; i++){
			int length = numDigits - x[i].length();
			if (length !=0 ) x[i] = zeros[length]+x[i];
		}
		return x;
	}
	
	/**Loads a HashSet with and Object array.*/
	public static HashSet loadHashSet(Object[] s){
		int num = s.length;
		HashSet hash = new HashSet(num);
		for (int i=0; i< num; i++) hash.add(s[i]);
		return hash;
	}
	
	/**Trims "" from the beginning and stop of String if both are found.  
	 * Excel when saving as a tab delimited file, adds quotes !$#%#$!!!*/
	public static String trimQuotes(String cell){
		if (cell.matches("^\".*\"$")) return cell.substring(1, cell.length()-1);
		return cell;
	}
	
	/**Splits a String on a regular expression, watching out for escaped characters
	 * ie Cat\=Mouse would not be split when = is given as the regex, a String[0], 
	 * {"Cat\=Mouse"} would be returned.  The method will not
	 * return empty Strings, Cat===Mouse would return {"Cat","Mouse}. Lastly,
	 * each String in the String[] is .trim()ed prior to saving, so no flanking 
	 * white space. Use this method when you need to allow for character escaping.
	 */
	public static String[] splitString (String concat, String regex){
		if (Misc.isEmpty(regex)) return new String[]{concat};
		
		Pattern p = Pattern.compile(regex);
		Matcher m = p.matcher(concat);
		ArrayList al = new ArrayList();
		int start = 0;
		int stop = 0;
		char slash = '\\';
		String token;
		
		//scan seq for match
		while (m.find()){
			//find a match
			stop = m.start();
			//System.out.println(stop);
			//is it escaped with a \
			if (stop==0) { start = 1;}
			else if (concat.charAt(stop-1) != slash ) {
				token = concat.substring(start, stop).trim();
				if (token.equals("") == false) al.add(token);
				start = stop +1;
			}
		}
		
		//add last
		token = concat.substring(start).trim();
		if (token.equals("") == false) al.add(token);

		//convert to String[]
		String [] cells = new String[al.size()];
		al.toArray(cells);
		return cells;
	}
	
	
	/**Loose matching of substring guess from larger String[]s.  
	 * Returns index of first matching choice, -1 if nothing found.
	 * Case insensitive. */
	public static int fetchMatchingStringArrayIndex(String guess, String[] choices){
		int num = choices.length;
		String lcGuess = guess.toLowerCase().trim();
		for (int i=0; i< num; i++){
			 if (choices[i].toLowerCase().indexOf(lcGuess)!=-1) return i;
		}
		return -1;
	}
	
	/**Lower cases a String[]*/
	public String[] toLowerCase(String[] strings){
		int num = strings.length;
		String[] ls = new String[num];
		for (int i=0; i<num; i++) ls[i] = strings[i].toLowerCase();
		return ls;
	}

	/**Prints a Object[] to System.out*/
	public static void printArray(Object[] array){
		if (array==null){
			System.out.println("null");
			return;
		}
		int len = array.length;
		for (int i=0; i<len; i++) System.out.println(array[i]);
	}
	
	/**Prints a long[] to System.out*/
	public static void printArray(long[] array){
		if (array==null){
			System.out.println("null");
			return;
		}
		int len = array.length;
		for (int i=0; i<len; i++) System.out.print(array[i]+"\t");
		System.out.println();
	}
	
	/**Prints a byte[] to System.out*/
	public static void printArray(byte[] array){
		if (array==null){
			System.out.println("null");
			return;
		}
		int len = array.length;
		for (int i=0; i<len; i++) System.out.print(array[i]+"\t");
		System.out.println();
	}
	
	/**Prints a short[] to System.out on one line*/
	public static void printArray(short[] array){
		if (array==null){
			System.out.println("null");
			return;
		}
		int len = array.length;
		System.out.print(array[0]);
		for (int i=1; i<len; i++) {
			System.out.print(" ");
			System.out.print(array[i]);
		}
		System.out.println();
	}
	
	/**Prints a boolean[] to System.out*/
	public static void printArray(boolean[] array){
		if (array==null){
			System.out.println("null");
			return;
		}
		int len = array.length;
		for (int i=0; i<len; i++) System.out.println(array[i]);
	}

	/**Prints a Object[][] to System.out*/
	public static void printArray(Object[][] array){
		if (array==null){
			System.out.println("null");
			return;
		}
		int len = array.length;
		for (int i=0; i<len; i++) {
			int len2= array[i].length;
			for (int j=0; j< len2; j++){
				System.out.println(array[i][j].toString());
			}
			System.out.println();
		}
	}
	
	/**Prints a long[][] to System.out*/
	public static void printArray(long[][] array){
		if (array==null){
			System.out.println("null");
			return;
		}
		int len = array.length;
		for (int i=0; i<len; i++) {
			int len2= array[i].length;
			for (int j=0; j< len2; j++){
				System.out.print(Long.toString(array[i][j]));
				System.out.print("\t");
			}
			System.out.println();
		}
	}

	/**Prints a String[] to System.out*/
	public static void printArray(String[] array){
		if (array==null){
			System.out.println("null");
			return;
		}
		int len = array.length;
		for (int i=0; i<len; i++) System.out.println(array[i]);
	}
	
	/**Prints a String[] to System.out*/
	public static void printArrayLine(String[] array){
		if (array==null){
			System.out.println("null");
			return;
		}
		int len = array.length;
		for (int i=0; i<len; i++) System.out.print(array[i]+" ");
		System.out.println();
	}

	/**Prints ArrayList to System.out*/
	public static void printArray(ArrayList array){
		if (array==null){
			System.out.println("null");
			return;
		}
		int len = array.size();
		for (int i=0; i<len; i++) System.out.println(array.get(i));
	}

	/**Prints a int[] to System.out*/
	public static void printArray(int[] array){
		if (array==null){
			System.out.println("null");
			return;
		}
		int len = array.length;
		for (int i=0; i<len; i++) System.out.print(array[i]+"\t");
		System.out.println();
	}
	/**Prints a char[] to System.out*/
	public static void printArray(char[] array){
		if (array==null){
			System.out.println("null");
			return;
		}
		int len = array.length;
		for (int i=0; i<len; i++) System.out.print(array[i]+" ");
		System.out.println();
	}
	/**Prints a char[] to System.out*/
	public static void printArrayWithTabs(char[] array){
		if (array==null){
			System.out.println("null");
			return;
		}
		int len = array.length;
		for (int i=0; i<len; i++) System.out.print(array[i]+"\t");
		System.out.println();
	}

	/**Prints a int[] to System.out*/
	public static void printArray(float[] array){
		if (array==null){
			System.out.println("null");
			return;
		}
		int len = array.length;
		for (int i=0; i<len; i++) System.out.print(array[i]+ " ");
		System.out.println();
	}	
	/**Prints a double[] to System.out on a single line.*/
	public static void printArray(double[] array){
		if (array==null){
			System.out.println("null");
			return;
		}
		int len = array.length;
		for (int i=0; i<len; i++) System.out.print(array[i]+"\t");
		System.out.println();
	}

	/**Splits an object[] into chunks containing the minNumEach. Any remainder is evenly distributed over the prior.
	 * Note this is by reference, the array is not copied. */
	public static Object[][] chunk (Object[] s, int minNumEach){
		//watch out for cases where the min can't be met
		int numChunks = s.length/minNumEach;
		if (numChunks == 0) return new Object[][]{s};
		
		double numLeftOver = (double)s.length % (double)minNumEach;
		
		int[] numInEach = new int[numChunks];
		for (int i=0; i< numChunks; i++) numInEach[i] = minNumEach;
		
		while (numLeftOver > 0){
			for (int i=0; i< numChunks; i++) {
				numInEach[i]++;
				numLeftOver--;
				if (numLeftOver == 0) break;
			}
		}
		//build chunk array
		Object[][] chunks = new Object[numChunks][];
		int index = 0;
		//for each chunk
		for (int i=0; i< numChunks; i++){
			//create container and fill it
			Object[] sub = new Object[numInEach[i]];
			for (int j=0; j< sub.length; j++) sub[j] = s[index++];
			chunks[i] = sub;
		}
		return chunks;
	}
	
	/**Splits an object[] into chunks containing a max number in each. Note this is by reference, the array is not copied. */
	public static Object[][] chunkWithMax (Object[] s, int maxInEach){
		//watch out for cases where the min can't be met
		int numChunks = s.length/maxInEach;
		if (numChunks == 0) return new Object[][]{s};

		double numLeftOver = (double)s.length % (double)maxInEach;
		//IO.pl(numChunks + " num chunks init");
		//IO.pl(numLeftOver +" num leftOver init");

		if (numLeftOver != 0) numChunks++;
		//IO.pl(numChunks + " num chunks post");

		//load with max
		int[] numInEach = new int[numChunks];
		for (int i=0; i< numChunks; i++) numInEach[i] = maxInEach;
		int numToSubtract = (numChunks*maxInEach) - s.length;
		//IO.pl(numToSubtract + " num to subtract");

		if (numToSubtract !=0) {
			int numSubtracted = 0;
			boolean loop = true;
			while (loop) {
				for (int i=0; i< numInEach.length; i++) {
					numInEach[i]--;
					numSubtracted++;
					if (numSubtracted == numToSubtract) {
						loop = false;
						break;
					}
				}
			}
		}
		//int totalObjects = 0;
		//for (int i=0; i< numInEach.length; i++) totalObjects+= numInEach[i];
		//IO.pl(totalObjects + " total objects\n");

		//build chunk array
		Object[][] chunks = new Object[numChunks][];
		int index = 0;
		//for each chunk
		for (int i=0; i< numChunks; i++){
			//create container and fill it
			Object[] sub = new Object[numInEach[i]];
			for (int j=0; j< sub.length; j++) sub[j] = s[index++];
			chunks[i] = sub;
		}
		return chunks;
	}


	/**Trims all the strings in the array String.trim()*/
	public static void trim(String[] s){
		for (int i=0; i< s.length; i++) s[i] = s[i].trim();
	}

	/**Removes the random digit_ proceeding a text, (ie 2432344322_fileName.txt to fileName.txt)*/
	public static String clipUnderScore(String string){
		int pos = string.indexOf("_");
		if (pos == -1) return string;
		return string.substring(pos+1);
	}

	/**Truncates a String adding three periods if too long.*/
	public static String truncate(String string, int cutOff){		
		if (string.length()>cutOff){
			StringBuffer sb = new StringBuffer (string.substring(0,cutOff-3));
			sb.append("...");
			return sb.toString();
		}
		return string;
	}
	
	/**Returns a nicely formated time, 15 May 2004 21:53.
	 * (Note, all I can say is that the GC DateFormat Date classes are so convoluted as to be utterly useless. Shame!)*/
	public static String getDateTime(){
		String[] months = {"Jan","Feb","Mar","Apr","May","June","July", "Aug","Sept","Oct","Nov","Dec"};
		GregorianCalendar c = new GregorianCalendar();
		int minutes = c.get(Calendar.MINUTE);
		String min;
		if (minutes < 10) min = "0"+minutes;
		else min = ""+minutes;
		return c.get(Calendar.DAY_OF_MONTH)+" "+months[c.get(Calendar.MONTH)]+" "+ c.get(Calendar.YEAR)+" "+c.get(Calendar.HOUR_OF_DAY)+":"+min;
	}

	/**Returns a nicely formated time (15 May 2004).
	 * (Note, all I can say is that the GC DateFormat Date classes are so convoluted as to be utterly useless. Shame!)*/
	public static String getDate(){
		String[] months = {"Jan","Feb","Mar","Apr","May","June","July", "Aug","Sept","Oct","Nov","Dec"};
		GregorianCalendar c = new GregorianCalendar();
		return c.get(Calendar.DAY_OF_MONTH)+ " "+months[c.get(Calendar.MONTH)]+ " "+ c.get(Calendar.YEAR);
	}
	/**Returns a nicely formated time (15May2004).
	 * (Note, all I can say is that the GC DateFormat Date classes are so convoluted as to be utterly useless. Shame!)*/
	public static String getDateNoSpaces(){
		String[] months = {"Jan","Feb","Mar","Apr","May","June","July", "Aug","Sept","Oct","Nov","Dec"};
		GregorianCalendar c = new GregorianCalendar();
		return c.get(Calendar.DAY_OF_MONTH)+months[c.get(Calendar.MONTH)]+ c.get(Calendar.YEAR);
	}
	
	/**Returns a formatted time 2:15:00 given 2.25*/
	public static String getFormattedTimeFromFraction (double fractionalTime){
		int whole = new Double (fractionalTime).intValue();
		double fraction = fractionalTime - whole;
		int min = (int)(fraction * 60);
		if (min<10) return whole+":0"+min+":00";
		return whole+":"+min+":00";
	}

	/**Takes a String like testSlide:33328:testslide desc...:1 or test: gip: 2,
	 * splits on the : and looks for ID in the last String[x], if not found it returns zero, otherwise
	 * it returns a parsed int.*/			
	public static int extractIdFromConcat(String concat){
		if (concat==null) return 0;
		String[] cells = concat.split(":");
		int numCells = cells.length;
		if (numCells!=0){
			String id = cells[numCells-1].trim();
			//check that a number is present
			try{
				return Integer.parseInt(id);
			} catch (NumberFormatException e){	
				System.err.println("Returning zero from extractIdFromConcat(), concat -> "+concat);
				return 0;
			}
		}
		return 0;
	}

	/**Takes a String like 1: testSlide: 33328:testslide desc etc or 2: test: gip,
	 * splits on the : and looks for ID in the first String[x], if not found it returns zero, otherwise
	 * it returns a parsed int.*/			
	public static int extractIdFromFrontOfConcat(String concat){
		if (concat==null) return 0;
		String[] cells = concat.split(":");
		int numCells = cells.length;
		if (numCells!=0){
			String id = cells[0].trim();
			//check that a number is present
			try{
				return Integer.parseInt(id);
			} catch (NumberFormatException e){
				System.err.println("Returning zero from extractIdFromFrontOfConcat(), concat -> "+concat);
				return 0;
			}
		}
		return 0;
	}

	/**Ditto to extractIDFromConcat but works on String[] returning IDs as a String[].
	 * Pulls off an int from the first cel of a colon delimited String.*/			
	public static String[] extractIDsFromConcats(String[] concats){
		int len = concats.length;
		String[] IDs = new String[len];
		for (int i=0; i<len; i++){
			IDs[i]=concats[i].split(":")[0].trim();
		} 
		return IDs;
	}

	/**Returns a String[2] split on first colon, if no colon returns null.*/	
	public static String[] splitOnFirstColon(String idColonName){
		int f = idColonName.indexOf(":");
		if (f==-1) return null;
		else return new String[]{idColonName.substring(0,f).trim(), idColonName.substring(f+1).trim()};
		
	}
	
	/**Method to check if a text is null or empty after trimming.*/
	public static boolean isEmpty(String text) {
		if (text == null || text.trim().equals("")) return true;
		return false;
	}

	/**Method to check if a text is not null and contains more than "".*/
	public static boolean isNotEmpty(String text) {
		if (text == null || text.trim().equals("")) return false;
		return true;
	}		
	/**Converts an ArrayList of Integer to int[]. Returns null if empty.*/
	public static int[] integerArrayListToIntArray(ArrayList integerAL){
		if (integerAL == null) return null;
		int num = integerAL.size();
		if (num ==0) return null;
		int[] ints = new int[num];
		try {
			for (int i=0; i< num; i++){
				ints[i]= ((Integer)integerAL.get(i)).intValue();
			}
			return ints;
		}
		catch(Exception e){
			System.err.println("Returning null from integerArrayListToIntArray(), ArrayList -> "+integerAL);
			return null;
		}
	}
	/**Returns a String[] given an ArrayList of Strings.*/
	public static String[] stringArrayListToStringArray(ArrayList stringAL){
		String[] s = new String[stringAL.size()];
		stringAL.toArray(s);
		return s;
	}
	
	/**Returns a String given an ArrayList of Character.*/
	public static String charArrayListToString(ArrayList<Character> al){
		char[] s = new char[al.size()];
		for (int i=0; i< s.length; i++) s[i] = al.get(i).charValue();
		return new String(s);
	}

	
	
	/**Returns a String separated by separator for each int[i].
	 * Returns "" if null or uninitialized.*/
	public static String intArrayToString(int[] i, String separator){
		if (i==null) return "";
		int len = i.length;
		if (len==1) return i[0]+"";
		if (len==0) return "";
		StringBuffer sb = new StringBuffer(i[0]+"");
		for (int j=1; j<len; j++){
			sb.append(separator);
			sb.append(i[j]);
		}
		return sb.toString();
	}		
	
	/**Returns a String separated by separator for each float[i].
	 * Returns "" if null or uninitialized.*/
	public static String floatArrayToString(float[] i, String separator){
		if (i==null) return "";
		int len = i.length;
		if (len==1) return i[0]+"";
		if (len==0) return "";
		StringBuffer sb = new StringBuffer(i[0]+"");
		for (int j=1; j<len; j++){
			sb.append(separator);
			sb.append(i[j]);
		}
		return sb.toString();
	}	
	
	/**Returns a String separated by separator for each double[i].
	 * Returns "" if null or uninitialized.*/
	public static String doubleArrayToString(double[] i, String separator){
		if (i==null) return "";
		int len = i.length;
		if (len==1) return i[0]+"";
		if (len==0) return "";
		StringBuilder sb = new StringBuilder(i[0]+"");
		for (int j=1; j<len; j++){
			sb.append(separator);
			sb.append(i[j]);
		}
		return sb.toString();
	}	


	/**Creates a LinkedHashMap from a String[] where the array is ordered such as key,value,key,value...
	 * Direction 0 forward, 1 = reverse of making the Map. Returns null if odd number of Strings*/
	public static LinkedHashMap createLinkedHashMap(String[] s, int direction){
		int len = s.length;
		LinkedHashMap hash= null;
		if (len % 2 !=0) return hash;
		hash = new LinkedHashMap(len/2);
		if (direction == 0){
			for (int i=len-1; i>=0; i= i-2){
				hash.put(s[i-1], s[i]);
			}
		}
		else {
			for (int i=len-1; i>=0; i= i-2){
				hash.put(s[i], s[i-1]);
			}			
		}
		return hash;
	}
	/**Prints a HashSet to a file. Each entry on a new line.*/
	public static boolean printHashSetToFile(HashSet hash, String file){
		try{
			PrintWriter out = new PrintWriter(new FileWriter(file));
			Iterator it = hash.iterator();
			while(it.hasNext()){
				out.println(it.next());
			}
			out.close();
			return true;
		}catch(Exception e){
			System.err.println("Problem with printHashSetToFile()");
			e.printStackTrace();
			return false;
		}
	}
	
	/**Converts an Object[] to an ArrayList*/
	public static ArrayList objectArrayToArrayList(Object[] ob){
		int len = ob.length;
		ArrayList al = new ArrayList(len);
		for (int i=0; i<len; i++){
			al.add(ob[i]);
		}
		return al;
	}
			
	/**Print to System.out a double[][]*/
	public static void printArray(double[][] d){
		int numColumns = d.length;
		int numRows = d[0].length;
		for (int i=0; i<numColumns; i++){
			System.out.println("\n\nColumn: "+i);
			for (int j=0; j<numRows; j++){
				System.out.println("Row: "+j+"\tValue: "+d[i][j]);
			}
		}
	}
	
	/**Print to System.out a double[][]*/
	public static void printArrayByRow(double[][] d){
		for (int i=0; i< d.length; i++){
			for (int j=0; j<d[i].length; j++){
				System.out.print(d[i][j]+"\t");
			}
			System.out.println();
		}
	}

	/**Print to System.out a float[][]*/
	public static void printArray(float[][] d){
		int numColumns = d.length;
		int numRows = d[0].length;
		for (int i=0; i<numColumns; i++){
			System.out.println("\n\nColumn: "+i);
			for (int j=0; j<numRows; j++){
				System.out.println("Row: "+j+"\tValue: "+d[i][j]);
			}
		}
	}
	
	/**Print to System.out a int[][]*/
	public static void printArray(int[][] d){
		if (d==null || d[0] == null){
			System.out.println("null");
			return;
		}
		int numColumns = d.length;
		for (int i=0; i<numColumns; i++){
			System.out.println("\nColumn: "+i);
			for (int j=0; j<d[i].length; j++){
				System.out.println("Row: "+j+"\tValue: "+d[i][j]);
			}
		}
	}	
	
	/**Returns an ArrayList of String given a String[]*/
	public static ArrayList stringArrayToArrayList(String[] s){
		int num = s.length;
		ArrayList al = new ArrayList(num);
		for (int i=0; i<num; i++){
			al.add(s[i]);
		}
		return al;
	}

	/**Returns a String separated by the separator given an ArrayList of String.*/
	public static String stringArrayListToString(ArrayList stringAL, String separator){
		int len = stringAL.size();
		if (len==0) return "";
		if (len==1) return (String)stringAL.get(0);
		StringBuffer sb = new StringBuffer((String)stringAL.get(0));
		for (int i=1; i<len; i++){
			sb.append(separator);
			sb.append((String)stringAL.get(i));
		}
		return sb.toString();
	}

	/**Returns a String separated by commas for each bin.*/
	public static String stringArrayToString(String[] s, String separator){
		if (s==null) return "";
		int len = s.length;
		if (len==1) return s[0];
		if (len==0) return "";
		StringBuffer sb = new StringBuffer(s[0]);
		for (int i=1; i<len; i++){
			sb.append(separator);
			sb.append(s[i]);
		}
		return sb.toString();
	}
	
	/**Returns a String separated by commas for each bin, skips nulls.*/
	public static String stringArrayToStringSkipNulls(String[] s, String separator){
		if (s==null) return "";
		int len = s.length;
		if (len==1) return s[0];
		if (len==0) return "";
		StringBuffer sb = new StringBuffer(s[0]);
		for (int i=1; i<len; i++){
			if (s[i] != null){
				sb.append(separator);
				sb.append(s[i]);
			}
		}
		return sb.toString();
	}
	
	/**Creates a HashMap from a String[] where the array is ordered such as key,value,key,value...
	 * Returns null if odd number of Strings*/
	public static HashMap createHashMap(String[] s){
		int len = s.length;
		HashMap hash= null;
		if (len % 2 !=0) return hash;
		hash = new HashMap(len/2);
		for (int i=len-1; i>=0; i= i-2){
			hash.put(s[i-1], s[i]);
		}
		return hash;
	}
    /**Returns true if a String contains just white space or is empty.*/
	public static boolean isStringEmpty(String string){
		if (string==null || string.trim().equals("")) return true;
		return false;
	}

	/**Short cut for debugging*/
	public static void p(String s) {
		System.err.println(s);
	}

	public static void invertArray(Object[] obj) {
		int lenthMinusOne = obj.length -1;
		int stop = obj.length/2;
		for (int i=0; i<stop; i++){
			int switchIndex = lenthMinusOne - i;
			Object first = obj[i];
			Object last = obj[switchIndex];
			obj[i] = last;
			obj[switchIndex] = first;
		}
		
	}

	/**Copies and merges.*/
	public static String[] copyAndMerge(String[] a, String[] b) {
		String[] c = new String[a.length+ b.length];
		System.arraycopy(a, 0, c, 0, a.length);
		System.arraycopy(b, 0, c, a.length, b.length);
		return c;
	}	
}
