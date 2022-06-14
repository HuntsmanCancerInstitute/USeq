package util.gen;

import java.io.File;
import java.util.ArrayList;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.text.similarity.LevenshteinDistance;
 
public class Delme {
 
	public static void main(String[] args) {
		
		String x= "2021-07-27 14:28:33       4236 result_json/result_bbaa32a7-4fb2-4b11-b805-deea248aa231.json";
		x=x.trim();
		Pattern jsonPat = Pattern.compile(".+\\sresult_json/.+\\.json");
		Matcher mat = pat.matcher(x);
		if (mat.matches()) IO.pl("Match ");

		
	}
}
	
	