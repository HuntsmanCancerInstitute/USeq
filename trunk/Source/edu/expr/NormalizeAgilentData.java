package edu.expr;
import java.io.*;
import java.util.*;
import util.gen.*;

/**Takes two replicas, parses out those probes that are not 1's, median centers, then averages
 * ProbeName	ProbeName	LogRatio	gIsPosAndSignif	rIsPosAndSignif	LogRatio	gIsPosAndSignif	rIsPosAndSignif
A_51_P100021	A_51_P100021	-0.087965	1	1	0.344989	0	1
A_51_P100034	A_51_P100034	0.038444	1	1	0.034538	1	1
A_51_P100052	A_51_P100052	-0.146955	1	0	0	0	0*/
public class NormalizeAgilentData {

	public NormalizeAgilentData(String[] args){
		
		try {
			//split data into flagged and non flagged
			File file = new File(args[0]);
			File save = new File (Misc.removeExtension(file.getCanonicalPath())+"Norm.xls");
			BufferedReader in = new BufferedReader(new FileReader(file));
			PrintWriter out = new PrintWriter (new FileWriter(save));
			ArrayList<DataLine> good = new ArrayList<DataLine>(); 
			ArrayList<DataLine> bad = new ArrayList<DataLine>();
			String line;
			while ((line = in.readLine()) != null){
				line = line.trim();
				if (line.length() == 0 || line.startsWith("Probe")) continue;
				DataLine data = new DataLine(line);
				if (data.flagged) bad.add(data);
				else good.add(data);
			}
			in.close();
			//convert to arrays
			DataLine[] goodData = new DataLine[good.size()];
			good.toArray(goodData);
			DataLine[] badData = new DataLine[bad.size()];
			bad.toArray(badData);
			
			//calculate median of good data
			double[] scoresA = new double[goodData.length];
			double[] scoresB = new double[goodData.length];
			double[] average = new double[goodData.length];
			for (int i=0; i< goodData.length; i++){
				scoresA[i] = goodData[i].scoreA;
				scoresB[i] = goodData[i].scoreB;
			}
			Arrays.sort(scoresA);
			Arrays.sort(scoresB);
			double medianA = Num.median(scoresA);
			double medianB = Num.median(scoresB);
			
			//median scale
			double scalarA = 1/medianA;
			double scalarB = 1/medianB;
			for (int i=0; i< goodData.length; i++){
				//scale
				goodData[i].scoreA *= scalarA;
				goodData[i].scoreB *= scalarB;
				//save to check
				scoresA[i] = goodData[i].scoreA;
				scoresB[i] = goodData[i].scoreB;
				//average 
				average[i] = (goodData[i].scoreA + goodData[i].scoreB)/2;
				goodData[i].log2Ave = Num.log2(average[i]);
			}
			
			Arrays.sort(scoresA);
			Arrays.sort(scoresB);
			Arrays.sort(average);
			double medianAPost = Num.median(scoresA);
			double medianBPost = Num.median(scoresB);
			double medianAve = Num.median(average);
			
			//stats
			String stats = 
				goodData.length + "\t# Good\n"+
				badData.length  + "\t# Bad\n"+
				medianA + "\tMedian A\n"+
				medianB + "\tMedian B\n"+
				scalarA + "\tScalar A\n"+
				scalarB + "\tScalar B\n"+
				medianAPost + "\tMedian A Post\n"+
				medianBPost + "\tMedian B Post\n"+
				medianAve + "\tMedian Ave\n";
			System.out.println("\n"+stats);
			
			//combine
			ArrayList<DataLine> combine = new ArrayList<DataLine>();
			for (int i=0; i< goodData.length; i++) combine.add(goodData[i]);
			for (int i=0; i< badData.length; i++) combine.add(badData[i]);
			DataLine[] data = new DataLine[combine.size()];
			combine.toArray(data);
			Arrays.sort(data);
			
			//print 
			for (int i=0; i< data.length; i++) out.println(data[i].name+"\t"+data[i].log2Ave);

			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		new NormalizeAgilentData (args);

	}
	
	private class DataLine implements Comparable {
		String name;
		double scoreA;
		double scoreB;
		double log2Ave;
		boolean flagged;
		public DataLine(String line){
			//load
			String[] tokens = line.split("\\t+");
			name = tokens[0];
			scoreA = Double.parseDouble(tokens[2]);
			scoreB = Double.parseDouble(tokens[5]);
			if (scoreA == 0 || scoreB == 0) {
				flagged = true;
				log2Ave = 0;
			}
			else flagged = false;
			//unlog scores
			scoreA = Num.antiLog(scoreA, 10);
			scoreB = Num.antiLog(scoreB, 10);
		}
		public int compareTo(Object obj){
			DataLine other = (DataLine) obj;
			return name.compareTo(other.name);
		}
	}

}
