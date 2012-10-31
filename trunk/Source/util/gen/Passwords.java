package util.gen;
import java.util.*;
/**
 * Methods for generating one time passwords.
 * @author nix
 *
 */
public class Passwords {
	//fields
	private String date;
	private String challengeConcatinate;
	private String htmlChallengeTable;
	
	public static String[] numberLetters = {
		"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", 
		"P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "1", "2", "3", "4", 
		"5", "6", "7", "8", "9", "0"};
	//alphabet minus 28 abiguous characters
	public static String[] nonAmbiguousLetters = {"A","B","C","D","E","F","G","H","J","K","L","M","N",
			"P","Q","R","T","U","V","W","X","Y","3","4","6","7","8","9"};		
	
	public void makeChallengeCard (){
		date = Misc.getDate();
		
		//start challenge text
		String [] randomWords = createRandomWords(nonAmbiguousLetters, 4, 100);
		StringBuffer challengeSB = new StringBuffer();
		challengeSB.append(date);
		int counter = 0;
		
		//make letter line
		StringBuffer ll = new StringBuffer("<TR ID='g'><TD>");
			for (int i=0; i<10; i++) {
				ll.append("<TD>");
				ll.append(numberLetters[i]);
			  }
			ll.append("<TR>\n");
		
		//start HTML Table	
		StringBuffer htmlTable = new StringBuffer(
			"<html><head><title>"+date+
			"</title><STYLE TYPE='text/css'>#g{background-color:rgb(200,200,200); font-weight:bold;} " +
			"#l{background-color:rgb(200,200,200);}"+
			"TD {font:arial; font-size:12; text-align:center;} BODY {font:arial; font-size:14;}</STYLE></head>\n<body>" +
			"BioRoot Challenge Card: &nbsp; &nbsp; " + date+ 
			"<table>" +ll);

		//loop through random words assigning to challenge text and HTML table
		for (int i=0; i<10; i++){
			htmlTable.append("<TR><TD ID='g'> &nbsp; "+(i+1)+" &nbsp; ");
			for (int j=0; j<10; j++){
				htmlTable.append("<TD>"+randomWords[counter]);
				challengeSB.append(":"+numberLetters[j]+(i+1)+"@"+randomWords[counter]);
				counter++;
			}
			htmlTable.append("</TR>");
		}
		htmlTable.append("</table></body></html>");	
		htmlChallengeTable = htmlTable.toString();
		challengeConcatinate = challengeSB.toString();
	}
	
	/**Creates pseudorandom Strings derived from an alphabet of String[] using the
	 * java.util.Random class.  Indicate how long you want a particular word and
	 * the number of words.*/
	public static String[] createRandomWords(String[] alphabet,int lengthOfWord,int numberOfWords) {
		ArrayList words = new ArrayList();
		Random r = new Random();
		int len = alphabet.length;
		for (int i = 0; i < numberOfWords; i++) {
			StringBuffer w = new StringBuffer();
			for (int j = 0; j < lengthOfWord; j++) {
				w.append(alphabet[r.nextInt(len)]);
			}
			words.add(w.toString());
		}
		String[] w = new String[words.size()];
		words.toArray(w);
		return w;
	}
	
	/**Returns a random word using nonambiguous alphabet.  Don't use this method for creating more than one word!*/
	public static String createRandowWord(int lengthOfWord){
		return createRandomWords(nonAmbiguousLetters, lengthOfWord,1)[0];
	}
	
	/**An example in how to make a huge random word table.*/
	public static void makeHugeTable(){		
	
		//text of file to write the html table
		String file = "randomWordTable.html";
		
		//make 5200 random 5 mers
		String[] randomWords = createRandomWords(nonAmbiguousLetters, 4, 5200);
		
		//make call and response A1-Z200 String[] and print table as html
		String[] callRes = new String[5200];
		String[] AthruZ = {"A","B","C","D","E","F","G","H","I","J","K","L","M","N",
		"O","P","Q","R","S","T","U","V","W","X","Y","Z"};
		int counter = 0;
		
		//make letter line
		StringBuffer ll = new StringBuffer("<TR ID='g'><TD>");
		int len = AthruZ.length;
		for (int i=0; i<len; i++) {
			ll.append("<TD>");
			ll.append(AthruZ[i]);
		}
		ll.append("<TR>\n");
		String letterLine = ll.toString();
		
		//print out html page
		Date date = new Date();
		StringBuffer htmlTable = new StringBuffer(
			"<html><head><title>"+date+
			"</title><STYLE TYPE='text/css'>#g{background-color:rgb(200,200,200);}"+
			"TD {font:arial; font-size:10; text-align:center;}</STYLE></head>\n<body>"+
			date+"<table>"+letterLine		
		);
		
		int trNum = 1;
		int hrNum = 1;
		for (int j=1; j<201; j++){
			if (hrNum==26){
				htmlTable.append("<TR><TD>&nbsp;</TD></TR>");
				htmlTable.append(letterLine);
				hrNum=1;
				j--;
			}
			else{
				hrNum++;
			
				if (trNum ==5){
					htmlTable.append("<TR ID='g'><TD>");
					trNum =1;
				}
				else {
					trNum++;
					htmlTable.append("<TR><TD>");
				} 
				htmlTable.append(j);
				for (int i=0; i<26; i++){
					htmlTable.append("<TD>");
					htmlTable.append(randomWords[counter]);
					callRes[counter]= AthruZ[i]+j+":"+randomWords[counter];
					counter++;
				}
				htmlTable.append("<TD>"+j+"</TR>\n");	
			}
		}
		htmlTable.append("</TABLE></body></html>");
		
		IO.writeString(htmlTable.toString(), file);
	}
	public static String[] getNonAmbiguousLetters() {
		return nonAmbiguousLetters;
	}

	public static String[] getNumberLetters() {
		return numberLetters;
	}

	public String getChallengeConcatinate() {
		return challengeConcatinate;
	}

	public String getDate() {
		return date;
	}

	public String getHtmlChallengeTable() {
		return htmlChallengeTable;
	}

}

