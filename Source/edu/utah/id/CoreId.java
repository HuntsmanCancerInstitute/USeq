package edu.utah.id;

import java.util.Random;
import java.util.regex.Pattern;

import util.gen.IO;

/** Generates and tests coreIds
 * an 8 letter number string in the following pattern, LLDLLDLL, L=A-Za-z but no OoIiLl, D=2-9, 262 billion combinations */
public class CoreId {
	
	
	//fields
	public static final String CORE_ID_DESCRIPTION = "LLDLLDLL, L=A-Za-z but no OoIiLl, D=2-9";
	public static final Pattern CORE_ID_Pattern = Pattern.compile("[A-Za-z&&[^OoIilL]]{2}[2-9][A-Za-z&&[^OoIilL]]{2}[2-9][A-Za-z&&[^OoIilL]]{2}");
	private static final String LETTERS = "ABCDEFGHJKMNPQRSTUVWXYZabcdefghjkmnpqrstuvwxyz";
	private static final String NUMBERS = "23456789";
	private Random random = null;
	
	public CoreId () {
		random = new Random();
	}
	
	public static boolean isCoreId(String testId) {
		return CORE_ID_Pattern.matcher(testId).matches();
	}
	
	public String createCoreId() {
		String sixLetters = getRandomString(6,LETTERS);
		String twoNumbers = getRandomString(2,NUMBERS);
		StringBuilder sb = new StringBuilder();
		sb.append(sixLetters.substring(0, 2));
		sb.append(twoNumbers.substring(0,1));
		sb.append(sixLetters.substring(2, 4));
		sb.append(twoNumbers.substring(1,2));
		sb.append(sixLetters.substring(4, 6));
		return sb.toString();
	}
	
	private String getRandomString(int length, String saltChars) {
        StringBuilder salt = new StringBuilder();
        while (salt.length() < length) {
            int index = (int) (random.nextFloat() * saltChars.length());
            salt.append(saltChars.charAt(index));
        }
        return salt.toString();
    }
	
	public static void main (String[] args) {
		CoreId coreId = new CoreId();
		for (int i=0; i< 10; i++) {
			String id = coreId.createCoreId();
			IO.pl(id+" "+CoreId.isCoreId(id));
		}
	}

}
