package trans.cel;

/**Class for holding info about a text cel file line, enables sorting after rotation.*/
public class CelLine implements Comparable{
	
	private int x;
	private int y;
	private String line;
	
	public CelLine (int x, int y, String line){
		this.x = x;
		this.y = y;
		this.line = line;
	}
	/**Sort first by y coordinate, then x.*/
	public int compareTo(Object obj){
		CelLine other = (CelLine) obj;
		if (other.y < y ) return 1;
		if (other.y > y ) return -1;
		if (other.x > x ) return -1;
		if (other.x < x ) return 1;




		//same
		return 0;
	}
	
	public String getLine(){
		return line;
	}

}
