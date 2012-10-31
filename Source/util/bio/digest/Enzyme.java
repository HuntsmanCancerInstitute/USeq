package util.bio.digest;
import java.util.regex.*;


/**
 * Restriction enzyme information.
 */
public class Enzyme implements Comparable{
	
	//Fields	
	private String name;
	private String recogSeq;
	private String strippedRecogSeq;
	
	private String methSens;
	private String avail;
	private String iso;
	
	private char overHangs;
	private int  recogLength;
	private int  cutTop;
	private int  cutBottom;
	private int[] topCuts;	//relative to a particular digested sequence
	private int	numberOfCuts;
	
	//Constructor
	public Enzyme (String aName, String aRecogSeq, String aMethSens, String aAvail, String aIsoz) {		
		name = aName;			
		recogSeq = aRecogSeq;	
		methSens = aMethSens;
		avail = aAvail;
		iso = aIsoz;
		
		//Calculate length of recogSeq
		//Strip out anything but ACGTRYMKSWBDHVN
		Pattern p = Pattern.compile("[^ACGTRYMKSWBDHVN]", Pattern.CASE_INSENSITIVE);
		Matcher m = p.matcher ( recogSeq );
		strippedRecogSeq = m.replaceAll("");
		recogLength = strippedRecogSeq.length();
		
		//Determine if the enzyme cut leaves a B=blunt, 3= 3' overhang 5= 5' overhang
		//Check to see if there is a ^ symbol, indicates palindrome
		double carrot = (double)recogSeq.indexOf("^"); 	//cast int to double
		if (carrot >= 0) {
			cutTop = (int)carrot;
			cutBottom = recogLength - (int)carrot;
		}    
		
		//Check to see if there is a (5/9), if two then it matches the second one, this is really cheating
		//(5/9) indicates a non palindrome cutter, some like BaeI (10/15)ACNNNNNGTAYC(12/7) cut twice
		//but this class will only take the second cut
		else {
			p = Pattern.compile("\\w+\\((-*\\d+)/(-*\\d{1,2})", Pattern.CASE_INSENSITIVE);
			m = p.matcher ( recogSeq );
			m.find();
			Integer i = new Integer (m.group(1));
			Integer j = new Integer (m.group(2));
			cutTop = recogLength + i.intValue();
			cutBottom = recogLength + j.intValue();
		}
		//Set the type of cut it leaves               
		if (cutTop > cutBottom) {
			overHangs = '3';
		}
		else if (cutTop < cutBottom) {
			overHangs = '5';
		}
		else {
			overHangs = 'B';
		}
	}	
	
	//Getter methods
	public String getName() {
		return name;
	}
	public String getRecogSeq() {	//NEB's format for recog seq  ie GAATCCC(3/5)
		return recogSeq;
	}
	public String getMethSens() {
		return methSens;
	}
	public String getAvail() {
		return avail;
	}
	public String getIso() {
		return iso;
	}	
	public int getCutTop() {	//this number referres to the base preceeding the cut  2  =  ga^tatc when
		return cutTop;		//the first base is 1
	}	
	public int getCutBottom() {	
		return cutBottom;
	}				
	public int getLength() {	//length of the recog seq GGGCCC = 6 
		return recogLength;
	}
	public String getRegEx() {	//ready to be put into Pattern.compile
		strippedRecogSeq = strippedRecogSeq.toUpperCase();
		StringBuffer sb = new StringBuffer();
		char base;
		String appender;
		for (int i=0; i< strippedRecogSeq.length(); i++) {
			base = strippedRecogSeq.charAt(i);
			switch (base) {
			case 'G':
			case 'A':
			case 'T':
			case 'C': sb.append(base); break;
			case 'N': sb.append("[GATC]"); break;
			case 'R': sb.append("[AG]"); break;
			case 'Y': sb.append("[CT]"); break;
			case 'M': sb.append("[AC]"); break;
			case 'K': sb.append("[GT]"); break;
			case 'S': sb.append("[CG]"); break;
			case 'W': sb.append("[AT]"); break;
			case 'B': sb.append("[CGT]"); break;
			case 'D': sb.append("[AGT]"); break;
			case 'H': sb.append("[ACT]"); break;
			}	
		}
		return new String(sb);
	}
	public char getOverHangs() {		
		return overHangs;
	}		
	
	public void printEnz() {
		System.out.println( "Enzyme text       : " + getName() + "\n" +
				"Recognition Seq   : " + getRecogSeq() + "\n" +
				"Stripped RecogSeq : " + strippedRecogSeq + "\n"+
				"Recognition Length: " + getLength() + "\n"+
				"Regular Expression: " + getRegEx() + "\n" +
				"Cuts Top Strand   : " + getCutTop() + "\n" +
				"Cuts Bot Strand   : " + getCutBottom() + "\n"+
				"Type of Overhang  : " + getOverHangs() + "\n"+
				"Isoschizomers     : " + getIso() + "\n"+
				"Methylation Sens  : " + getMethSens() + "\n"+
				"Availibility      : " + getAvail() + "\n");
	}	
	public int[] getTopCuts() {
		return topCuts;
	}
	public void setTopCuts(int[] topCuts) {
		this.topCuts = topCuts;
		if (topCuts!= null) numberOfCuts = topCuts.length;
		else numberOfCuts = 0;
	}

	public int compareTo(Object obj) {
		Enzyme other = (Enzyme) obj;
		//sort by number of cuts
		if (other.numberOfCuts<numberOfCuts) return 1;
		if (other.numberOfCuts>numberOfCuts) return -1;
		//sort by position
		if (numberOfCuts!=0){
			if (other.getTopCuts()[0]< topCuts[0]) return 1;
			if (other.getTopCuts()[0]> topCuts[0]) return -1;
		}
		//sort by alphabetically
		return other.name.compareTo(name)*-1;
	}
	public int getNumberOfCuts() {
		return numberOfCuts;
	}
	public String getStrippedRecogSeq() {
		return strippedRecogSeq;
	}
}

