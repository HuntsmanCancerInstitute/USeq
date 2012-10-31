package trans.misc;
import java.io.*;
import java.util.*;
import java.util.regex.*;

import util.gen.*;
import util.bio.annotation.*;

public class Delme {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Pattern pat = Pattern.compile("=");
		String x = "chrom=chr19";
		
		String[] t = pat.split(x);
		
		
		
		System.out.println(t[1]);
	}
	
	
}
