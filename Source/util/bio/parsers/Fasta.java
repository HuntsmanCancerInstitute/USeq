package util.bio.parsers;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Fasta implements Comparable{
	
	//fields
	private String seq;
	private String name;
	private int number =0;
	private final Pattern numberPat = Pattern.compile("\\d+$");

	//constructor
	public Fasta (String name, String seq){
		this.seq = seq;
		this.name = name;
		//attempt to extract a number
		Matcher mat = numberPat.matcher(name);
		if (mat.find()) number = Integer.parseInt(mat.group());
	}
	
	//for sorting by extracted number
	public int compareTo (Object obj){
		Fasta other = (Fasta)obj;
		if (other.number < number) return 1;
		if (other.number > number) return -1;
		return 0;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public int getNumber() {
		return number;
	}

	public void setNumber(int number) {
		this.number = number;
	}

	public String getSeq() {
		return seq;
	}

	public void setSeq(String seq) {
		this.seq = seq;
	}
}
