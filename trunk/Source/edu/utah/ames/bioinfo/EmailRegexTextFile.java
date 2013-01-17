package edu.utah.ames.bioinfo;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * 
 * @author darren.ames@hci.utah.edu
 */
public class EmailRegexTextFile {
	
	//fields
	static String EMAIL = "^#e\\s+(.+@.+)";
	
	public static void main(String[] args) throws IOException {

		// expression to find a valid e-mail address within a file
		Pattern pattern = Pattern.compile(EMAIL, Pattern.CASE_INSENSITIVE);

		// read file
		FileInfo fi = new FileInfo();
		File file = fi.getCmdFile();

		BufferedReader br = new BufferedReader(new FileReader(file));
		int lines = 0;
		int matches = 0;

		// search for email address String
		for (String line = br.readLine(); line != null; line = br.readLine()) {
			lines++;
			Matcher matcher = pattern.matcher(line);
			if (matcher.find()) {
				// System.out.println(matcher.group(1));
				matches++;

				// set email address in EmailAddress class
				EmailAddress email = new EmailAddress(matcher.group(1));
				System.out.println(email.getEmailAddress());
			}
		}

		// output summary
		if (matches == 0) {
			System.out.println("No matches for valid email address" + " in "
					+ file);

		} else {
			System.out.println("\n" + matches + " match(es) in " + file);
		}
	}
}
