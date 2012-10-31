package trans.cel;
import java.io.*;

public class ControlStats implements Serializable{
	
	private static final long serialVersionUID = 0;
	private double[] stats;
	private String[] names; 
	
	public boolean ordered(){
		for (int x=0; x< stats.length-1; x++) {
			if (stats[x]<= stats[x+1] == false) return false;
		}
		return true;
	}
	
	public String toString (){
		StringBuffer sb = new StringBuffer();
		sb.append(stats[0]);
		for (int x=1; x< stats.length; x++) {
			sb.append("\t");
			sb.append(stats[x]);
		}
		return sb.toString();
	}

	public String[] getNames() {
		return names;
	}

	public void setNames(String[] names) {
		this.names = names;
	}

	public double[] getStats() {
		return stats;
	}

	public void setStats(double[] stats) {
		this.stats = stats;
	}

}
