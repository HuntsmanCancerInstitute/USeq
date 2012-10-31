package meme;

/**Holds data related to a particular match of a matrix to a sequence. Interbase coordinates.*/
public class MotifHit implements Comparable{
    private int start;  //the first base  0
    private int stop;    //not included
    private double score;
    private String seq;	
    private int orientation; //0 is forward, 1 reverseComplement
    
    
    public MotifHit(int st, int sp, double sc, String s, int o){
        start = st;
        stop = sp;
        score = sc;
        seq = s;
        orientation = o;
    }
    public int compareTo(Object otherObject){
        MotifHit other = (MotifHit) otherObject;
        //sort first by score
        if (score > other.score) return -1;
        if (score < other.score) return 1;
        // if scores same, sort by postion
        if (start > other.start) return -1;
        if (start < other.start) return 1;
        return 0;
    }
    public String toString(){
        return
        "** MotifHit ** Seq: "+seq+ "  Score: "+score + "  Start: "+start+"  Stop : "+stop+"  Ori : "+orientation;
    }
	public double getScore() {
		return score;
	}
	public void setScore(double score) {
		this.score = score;
	}
	public String getSeq() {
		return seq;
	}
	public void setSeq(String seq) {
		this.seq = seq;
	}
	public int getStart() {
		return start;
	}
	public void setStart(int start) {
		this.start = start;
	}
	public int getStop() {
		return stop;
	}
	public void setStop(int stop) {
		this.stop = stop;
	}
    
	public int getOrientation() {
		return orientation;
	}
}

