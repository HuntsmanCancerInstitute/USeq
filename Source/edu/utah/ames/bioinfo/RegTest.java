package edu.utah.ames.bioinfo;

/*
 * Tests regular expressions against strings via console input. Exit
 * by entering a blank string and selecting "N/n".
 */

import java.util.regex.*;
import java.util.Scanner;
/**
 * 
 * @author darren.ames@hci.utah.edu
 *
 */

public final class RegTest {

	static String r, s;
	static Pattern pattern;
	static Matcher matcher;
	static boolean match, validRegex, doneMatching;
	
	private static Scanner sc = new Scanner(System.in);
	
	public static void main(String[] args) {
		
		System.out.println("Welcome to the regex tester\n");
		
		do {
			do {
				System.out.print("Enter regex: ");
				r = sc.nextLine();
				validRegex = true;
				try {
					pattern = Pattern.compile(r);
				}
				catch (Exception e) {
					System.out.println(e.getMessage());
					validRegex = false;
				}
			}
			while (!validRegex);
			doneMatching = false;
			
			while (!doneMatching) {
				System.out.print("Enter a string to test against: ");
				s = sc.nextLine();
				if (s.length() == 0) 
					doneMatching = true;
				else {
					matcher = pattern.matcher(s);
					if (matcher.matches()) 
						System.out.println("Match.");
					else
						System.out.println("Not a match.");
				}
			}
		}
		while (askAgain());
	}
	
	private static boolean askAgain() {
		System.out.println("Another? (Y or N)");
		String reply = sc.nextLine();
		if (reply.equalsIgnoreCase("Y"))
			return true;
		return false;
	}
}
