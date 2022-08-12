package util.gen;

import java.io.File;
import java.util.ArrayList;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.text.similarity.LevenshteinDistance;

import edu.utah.seq.run.PlatformGenderInfo;
 
public class Delme {
 
	public static void main(String[] args) {
		
		String x= "TL-18-843E9B_XT.V1_2018-10-26_gz_Neeraj_Agarwal_F.json.txt.gz";
		
		IO.pl(x);
		IO.pl(removeExtension(x));

		
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
	
	
}
	
	