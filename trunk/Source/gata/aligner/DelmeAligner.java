package gata.aligner;

import java.util.regex.*;
public class DelmeAligner {

	public static void main(String[] args) {
		Pattern pat = Pattern.compile("\\d+");
		Matcher mat = pat.matcher("dme.xxy_321.2-4445_gipers_555968");
		if (mat.find()) System.out.println(mat.group());
		else System.out.println("no find");
	}
}
