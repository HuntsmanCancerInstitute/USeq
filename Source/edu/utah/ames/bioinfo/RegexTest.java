package edu.utah.ames.bioinfo;

import java.util.regex.*;
import java.util.Scanner;
/**
 * 
 * @author darren.ames@hci.utah.edu
 *
 */

public class RegexTest {

	static String r;
	static String s;
	static Pattern p;
	static Matcher m;
	static boolean match;
	static boolean validRegex;
	static boolean doneMatching;

	private static Scanner sc = new Scanner(System.in);

	public static void main(String[] args) {

		System.out.println("Welcome to the regex tester\n");

		do {
			do {
				System.out.print("Enter regex: ");
				r = sc.nextLine();
				validRegex = true;
				try {
					p = Pattern.compile(r);
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
					m = p.matcher(s);
					if (m.find()) { 
						//System.out.println("Match.");
						System.out.println(m.group(0));
						System.out.println(m.group(1));
					}
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
	/**
	String data="300x250,468x60,300x400v(480x320,768x1024,100x100),400x300v,640x480v(200x200,728x90)";

	List<String> tokens=new ArrayList<>();
	StringBuilder buffer=new StringBuilder();

	int parenthesesCounter=0;

	for (char c : data.toCharArray()){
	    if (c=='(') parenthesesCounter++;
	    if (c==')') parenthesesCounter--;
	    if (c==',' && parenthesesCounter==0){
	        //lets add token inside buffer to our tokens
	        tokens.add(buffer.toString());
	        //now we need to clear buffer  
	        buffer.delete(0, buffer.length());
	    }
	    else 
	        buffer.append(c);
	}
	//lets not forget about part after last comma
	tokens.add(buffer.toString());

	String[] splitedArray=tokens.toArray(new String[tokens.size()]);

	//lets test what is inside our array
	for (String s : splitedArray)
	    System.out.println(s);
	    */
}

