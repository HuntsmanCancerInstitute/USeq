package edu.expr;
import java.io.*;
import util.gen.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ParseMiRNATable {
	
	//fields
	private File sangerMiRNAsTable;
	private File specialMiRNANamesFile;
	private HashMap<String,String> specialMiRNANames;
	private HashMap miRNAs;
	private HashMap targetScores;
	private HashMap<String,String> ensblTransID2GeneName;
	private HashMap<String,HashSet<String>> geneTargetMiRNAs;
	private String[] specialMiRNANamesNotFound;
	private Pattern colonNeg = Pattern.compile(":-");
	
	public ParseMiRNATable (String[] args){
		sangerMiRNAsTable = new File (args[0]);
		//create hash of miRNAs (text:MiRNA w/targets)
		parseMiRNATable();
		//load list of select miRNAs from file
		specialMiRNANamesFile = new File (args[1]);
		specialMiRNANames = IO.loadFileIntoHashMap(specialMiRNANamesFile);
		//load gene text look up hash
		ensblTransID2GeneName = IO.loadFileIntoHashMap(new File(args[2]));
		//sum target scores
		sumTargetScores();
		//make result lines saving in ArrayList
		ArrayList<ScoreLine> resultLines = new ArrayList<ScoreLine>();
		Iterator<String> it = targetScores.keySet().iterator();
		while (it.hasNext()){
			String geneTargetName = it.next();
			float sumScore = ((Float)(targetScores.get(geneTargetName))).floatValue();
			HashSet<String> mis = geneTargetMiRNAs.get(geneTargetName);
			//any alternative gene text present?
			if (ensblTransID2GeneName.containsKey(geneTargetName)){
				geneTargetName = geneTargetName+"\t"+ensblTransID2GeneName.get(geneTargetName);
			}
			else geneTargetName = geneTargetName+"\t";
			//split miRNA's by score sign
			String[] plusMinus = splitByScoreSign(Misc.hashSetToStringArray(mis));
			String summaryLine = geneTargetName+"\t"+sumScore+"\t"+plusMinus[0]+ "\t"+ plusMinus[1];
			resultLines.add(new ScoreLine(sumScore, summaryLine));
		}
		//convert results AL into an array and sort
		ScoreLine[] sl = new ScoreLine[resultLines.size()];
		resultLines.toArray(sl);
		Arrays.sort(sl);
		//print results
		System.out.println("EnsemblTransName\tGeneName\tSumScore\tMiRNAs:PlusFoldChange\tMiRNAs:MinusFoldChange");
		for (int i=0; i< sl.length; i++) System.out.println(sl[i].line);
		System.out.println("\nmiRNAs not found: "+Misc.stringArrayToString(specialMiRNANamesNotFound, ", "));
	}
	
	private class ScoreLine implements Comparable{
		float score;
		String line;
		
		public ScoreLine (float score, String line){
			this.score = score;
			this.line = line;
		}
		
		public int compareTo (Object obj){
			ScoreLine other = (ScoreLine)obj;
			if (other.score > score) return 1;
			if (other.score < score) return -1;
			return 0;
		}
	}
	
	public String[] splitByScoreSign(String[] items){
		ArrayList<String> plus = new ArrayList<String>();
		ArrayList<String> minus = new ArrayList<String>();
		for (int i=0; i< items.length; i++){
			Matcher mat = colonNeg.matcher(items[i]);
			if (mat.find()) minus.add(items[i]);
			else plus.add(items[i]);
		}
		return new String[]{Misc.stringArrayListToString(plus, ", "), Misc.stringArrayListToString(minus, ", ")};
	}

	public void sumTargetScores (){
		targetScores = new HashMap();
		ArrayList notFound = new ArrayList();
		geneTargetMiRNAs = new HashMap<String,HashSet<String>>();
		//for each miRNA in users hit list
		Iterator<String> it = specialMiRNANames.keySet().iterator();
		while (it.hasNext()){
		//for (int i=0; i< specialMiRNANames.length; i++){
			String namePickedMiRNA = it.next();
			//is it in the global list
			if (miRNAs.containsKey(namePickedMiRNA)) {
				MiRNA miRNA = (MiRNA)miRNAs.get(namePickedMiRNA);
				//for each gene target
				MiRNATarget[] targets = miRNA.getTargets();
				Float sum = new Float(0);
				for (int j=0; j< targets.length; j++){
					//is it in the hash?
					String targetName = targets[j].getName();
					if (targetScores.containsKey(targetName) == false){
						sum = new Float(targets[j].getScore());
					}
					else {
						sum = (Float)targetScores.get(targetName);
						sum = new Float (sum.floatValue() + targets[j].getScore());
					}
					targetScores.put(targetName, sum);
					//save miRNA hit in targets
					if (geneTargetMiRNAs.containsKey(targetName)){
						HashSet<String> mis = geneTargetMiRNAs.get(targetName);
						mis.add(namePickedMiRNA+":"+specialMiRNANames.get(namePickedMiRNA));
					}
					else {
						HashSet<String> mis = new HashSet<String>();
						mis.add(namePickedMiRNA+":"+specialMiRNANames.get(namePickedMiRNA));
						geneTargetMiRNAs.put(targetName, mis);
					}
				}
			}
			else {
				notFound.add(namePickedMiRNA);
			}
		}
		specialMiRNANamesNotFound = new String[notFound.size()];
		notFound.toArray(specialMiRNANamesNotFound);
	}
	
	public void parseMiRNATable(){
		//make hash of miRNAName and ArrayList of hits
		HashMap map = hashMiRNATable(sangerMiRNAsTable);
		//make MiRNAs
		MiRNA[] miRNAsArray = makeMiRNAs(map);
		//make hash based on text
		miRNAs = new HashMap();
		for (int i=0 ;i<miRNAsArray.length; i++) {
			String name = miRNAsArray[i].getName();
			if (miRNAs.containsKey(name)) Misc.printExit("\nError: duplicate miRNA text! -> "+name );
			miRNAs.put(name, miRNAsArray[i]);
		}
	}
	

	
	public static MiRNA[] makeMiRNAs(HashMap map){
		MiRNA[] miRNAs = new MiRNA[map.size()];
		Iterator it = map.keySet().iterator();
		int index = 0;
		while (it.hasNext()){
			String miRNAName = (String) it.next();
			ArrayList hits = (ArrayList) map.get(miRNAName);
			MiRNATarget[] targets = new MiRNATarget[hits.size()];
			hits.toArray(targets);
			miRNAs[index++] = new MiRNA (miRNAName, targets);
		}
		return miRNAs;
	}
	
	public static HashMap hashMiRNATable(File txtTable){
		//make hash map to hold objects
		HashMap map = new HashMap();
		//run through lines
		try {
			String line;
			String[] tokens;
			ArrayList al;
			BufferedReader in = new BufferedReader (new FileReader (txtTable));
			while ((line=in.readLine()) != null){
				line = line.trim();
				if (line.startsWith("#") || line.length() ==0) continue;
				tokens = line.split("\\t");
				//System.out.println(line+"\t"+tokens.length);
				String miRNAName = tokens[1];
				//look to see if already exists
				if (map.containsKey(miRNAName)){
					 al = (ArrayList) map.get(miRNAName);
				}
				else {
					al = new ArrayList();
					map.put(miRNAName, al);
				}
				//make a hit and put in ArrayList, String text, String alias, String score, String pvalue
				String alias = "";
				if (tokens.length > 12) alias = tokens[12];
				MiRNATarget hit = new MiRNATarget(tokens[11], alias, tokens[9], tokens[10]);
				al.add(hit);
			}
			in.close();
			return map;
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}

	}
	public static void main (String[] args){
		if (args.length == 0) Misc.printExit("\nSangerMiRNATableDBDumpFile SelectTargetsFile(text score) EnsTransID2GeneNameFile");
		new ParseMiRNATable (args);
	}

}
