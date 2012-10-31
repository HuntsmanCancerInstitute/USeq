package util.bio.wrappers;
import java.io.*;
import util.gen.*;
import java.util.*;

/**Container for sequence used in picking primers with Primer3.*/
public class Primer3Pick {

	private String sequence;
	private String line;
	private String resultSummaryLine = "";
	public static final String header = "PRIMER_LEFT_SEQUENCE\tPRIMER_RIGHT_SEQUENCE\tPRIMER_LEFT_TM\tPRIMER_RIGHT_TM\tPRIMER_PRODUCT_SIZE\tSequence\t";
	
	public Primer3Pick (String line){
		String[] tokens = line.split("\\s");
		sequence = tokens[0];
		this.line = line;
	}
	
	/**Makes a tab delimited results line:
	 * PRIMER_LEFT_SEQUENCE
	 * PRIMER_RIGHT_SEQUENCE
	 * PRIMER_LEFT_TM
	 * PRIMER_RIGHT_TM
	 * PRIMER_PRODUCT_SIZE
	 * The user entered line
	 * 
	 * If no results then leaves no picks message.
	 * 
	 * */
	public void loadResults( HashMap results){
		StringBuffer sb = new StringBuffer ();
		//sucess?
		if (results.containsKey("PRIMER_LEFT_SEQUENCE") == false){
			sb.append("No picks\t\t\t\t\t");
			sb.append(line);
		}
		else {
			sb.append(Misc.fetchExit("PRIMER_LEFT_SEQUENCE", results));
			sb.append("\t");
			sb.append(Misc.fetchExit("PRIMER_RIGHT_SEQUENCE", results));
			sb.append("\t");
			sb.append(Misc.fetchExit("PRIMER_LEFT_TM", results));
			sb.append("\t");
			sb.append(Misc.fetchExit("PRIMER_RIGHT_TM", results));
			sb.append("\t");
			sb.append(Misc.fetchExit("PRIMER_PRODUCT_SIZE", results));
			sb.append("\t");
			sb.append(line);
		}
		resultSummaryLine = sb.toString();
	}
	
	/**Parses a tab delimited file where the sequence to pick primers from is in column 0.*/
	public static Primer3Pick[] parseFile(File file){
		ArrayList al = new ArrayList();
		Primer3Pick[] p = null;
		try{
			BufferedReader in = new BufferedReader(new FileReader(file));
			String line;
			while ((line = in.readLine()) != null){
				line = line.trim();
				if (line.startsWith("#") || line.length()==0) continue;
				al.add(new Primer3Pick(line));
			}
			in.close();
			p = new Primer3Pick[al.size()];
			al.toArray(p);
		}catch (Exception e){
			e.printStackTrace();
		}
		return p;
	}

	public String getLine() {
		return line;
	}

	public String getSequence() {
		return sequence;
	}

	public static String getHeader() {
		return header;
	}

	public String getResultSummaryLine() {
		return resultSummaryLine;
	}

	public void setResultSummaryLine(String resultSummaryLine) {
		this.resultSummaryLine = resultSummaryLine;
	}

	public void setLine(String line) {
		this.line = line;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

}
