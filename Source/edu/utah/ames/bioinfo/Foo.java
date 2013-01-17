package edu.utah.ames.bioinfo;

import java.text.ParseException;
import java.text.SimpleDateFormat;

public class Foo {
	public static void main(String[] args) throws ParseException {
		long time = System.currentTimeMillis();
		System.out.println(time);
		
		long fiveDaysAgo = time - (5L * 24 * 60 * 60 * 1000);
		long sevenDaysAgo = time - (7L * 24 * 60 * 60 * 1000);
		
		System.out.println("5 days ago: " + fiveDaysAgo);
		System.out.println("7 days ago: " + sevenDaysAgo);
		
		SimpleDateFormat formatter = new SimpleDateFormat("yyyy-MM-dd");
		
		long d1 = formatter.parse("2001-1-1").getTime();
		long d2 = formatter.parse("2012-1-1").getTime();
		
		System.out.println(Math.abs((d1-d2)/(1000*60*60*24)));
	}
}
