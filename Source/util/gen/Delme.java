package util.gen;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Delme {

	public static void main(String[] args) throws IOException {
		String x = "DP=1134;TI=NM_003482.3;GI=KMT2D;FC=Nonsense;PC=Q3812*;DC=c.11434C>T;LO=EXON;EXON=39;CI=";
		Pattern pat = Pattern.compile(".+TI=(NM_[\\d\\.]+).+");
		Matcher mat = pat.matcher(x);
		if (mat.matches()) IO.pl("Found:"+mat.group(1));
		else IO.pl("NoMatch");

	}

}
