package util.bio.parsers;
import java.util.*;
import util.gen.*;

public class BaseCount {
	
	//fields
	private int position;
	private char base;
	private int g;
	private int a;
	private int t;
	private int c;
	private int n;
	private int deletions;
	private ArrayList <Integer> insertionLengths = new ArrayList <Integer> ();
	
	
	public BaseCount (int position, char base){
		this.position = position;
		this.base = base;
		countBase(base);
	}
	
	/**Base must be lower case, if not gatc then n is incremented.*/
	public void countBase(char b){
		if (b == 'g') g++;
		else if (b == 'a') a++;
		else if (b == 't') t++;
		else if (b == 'c') c++;
		else n++;
	}
	
	//Position Base G A T C N Del Int
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append(position); sb.append("\t");
		sb.append(base); sb.append("\t");
		sb.append(g); sb.append("\t");
		sb.append(a); sb.append("\t");
		sb.append(t); sb.append("\t");
		sb.append(c); sb.append("\t");
		sb.append(n); sb.append("\t");
		sb.append(deletions); sb.append("\t");
		if (insertionLengths.size() != 0){
			//int[] il = Misc.integerArrayListToIntArray(insertionLengths);
			//sb.append(Misc.intArrayToString(il, ","));
			sb.append(insertionLengths.size());
		}
		else sb.append("0");
		
		return sb.toString();
	}

	public int getA() {
		return a;
	}

	public void setA(int a) {
		this.a = a;
	}

	public char getBase() {
		return base;
	}

	public void setBase(char base) {
		this.base = base;
	}

	public int getC() {
		return c;
	}

	public void setC(int c) {
		this.c = c;
	}

	public int getDeletions() {
		return deletions;
	}

	public void setDeletions(int deletions) {
		this.deletions = deletions;
	}

	public int getG() {
		return g;
	}

	public void setG(int g) {
		this.g = g;
	}

	public ArrayList <Integer> getInsertionLengths() {
		return insertionLengths;
	}

	public void setInsertionLengths(ArrayList<Integer> insertionLengths) {
		this.insertionLengths = insertionLengths;
	}

	public int getN() {
		return n;
	}

	public void setN(int n) {
		this.n = n;
	}

	public int getPosition() {
		return position;
	}

	public void setPosition(int position) {
		this.position = position;
	}

	public int getT() {
		return t;
	}

	public void setT(int t) {
		this.t = t;
	}
	
 
}
