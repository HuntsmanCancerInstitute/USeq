package trans.main;
import java.io.*;
import java.util.regex.*;
import util.gen.*;

public class ExportIntervalData {
	
	//fields
	private File[] intervalFiles;
	private Interval[] intervals;
	private float[][] treatments; //[replica index][intensities]
	private float[][] controls;
	private int numTreatmentReplicas;
	private int numControlReplicas;
	
	public ExportIntervalData(String[] args){
		processArgs(args);
		//for each interval file
		for (int i=0; i< intervalFiles.length; i++){
			try {
				//load interval array and check for oligos
				if (loadIntervals(intervalFiles[i]) == false) Misc.printExit("\nError: cannot load an interval array from "+intervalFiles[i].getName()+"\n");
				if (intervals[0].getOligos() == null) {
					System.out.println("\nError: cannot extract intensities for "+intervalFiles[i].getName()+". Need oligo info.  Run LoadIntervalOligoInfo app first.\n");
					continue;
				}
				numTreatmentReplicas = intervals[0].getNumberTreatmentIntensities();
				numControlReplicas = intervals[0].getNumberControlIntensities();
				
				//print header
				File fileOut = new File (intervalFiles[i].getCanonicalPath()+".sum");
				PrintWriter out = new PrintWriter (new FileWriter (fileOut));
				//out.println(getHeader());
				
				//for each interval
				for (int j=0; j<intervals.length; j++){
					//load best window oligo intensities
					loadBestWindowOligoIntensities (intervals[j]);
					//use pseudoMedian to summarize intensities
					double[] treatmentSummaries = Num.pseudoMedian(treatments);
					double[] controlSummaries = Num.pseudoMedian(controls);
					//print
					StringBuffer sb = new StringBuffer();
					sb.append(j);
					for (int k=0; k<numTreatmentReplicas; k++) {
						sb.append("\t");
						sb.append(treatmentSummaries[k]);
						
					}
					for (int k=0; k<numControlReplicas; k++) {
						sb.append("\t");
						sb.append(controlSummaries[k]);
					}
					out.println(sb);
				}
				out.close();
			} catch (IOException e){
				e.printStackTrace();
			}
		}
	}
	
	public static String getHeader(int numTreatmentReplicas, int numControlReplicas){
		StringBuffer sb = new StringBuffer("#Name");
		for (int j=0; j<numTreatmentReplicas; j++) {
			sb.append("\t");
			sb.append("T");
			sb.append(j);
		}
		for (int j=0; j<numControlReplicas; j++) {
			sb.append("\t");
			sb.append("C");
			sb.append(j);
		}
		return sb.toString();
	}
	
	
	/**Loads the treatment and control intensity arrays for a particular interval.*/
	public void loadBestWindowOligoIntensities(Interval interval){
		//for each Interval
		Window best = interval.getBestWindow();
		Oligo[] oligos = interval.extractSubRegionOligos(best.getStart1stOligo(), best.getStartLastOligo());
		int numOligos = oligos.length;
		//get intensities[replica index][intensities]
		float[][] treatments = new float[numTreatmentReplicas][numOligos];
		float[][] controls = new float[numControlReplicas][numOligos];
		//for each oligo
		for (int i=0; i<numOligos; i++){
			float[] t = oligos[i].getTreatmentIntensities(interval.getNumberTreatmentIntensities());
			float[] c = oligos[i].getControlIntensities(interval.getNumberControlIntensities());	
			//for each replica
			for (int j=0; j< numTreatmentReplicas; j++) treatments[j][i] = t[j];
			for (int j=0; j< numControlReplicas; j++) controls[j][i] = t[j];
		}
	}
	
	public boolean loadIntervals( File intervalFile){
		try {
			intervals = (Interval[])IO.fetchObject(intervalFile);
		} catch (Exception e){
			return false;
		}
		return true;
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}	
		new ExportIntervalData(args);
	}		
	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'i': intervalFiles = IO.extractFiles(new File (args[i+1])); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (intervalFiles == null || intervalFiles.length ==0) Misc.printExit("\nError: cannot find your interval file(s)!\n");
	}	
	
	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Aug  2006                                **\n" +
				"**************************************************************************************\n" +

		"**************************************************************************************\n");		
	}
}