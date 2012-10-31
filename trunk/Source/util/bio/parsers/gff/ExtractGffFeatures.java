package util.bio.parsers.gff;
import java.io.*;
import java.util.regex.*;
import util.gen.*;
import java.util.*;

public class ExtractGffFeatures {

	private File[] gffFiles;

	public ExtractGffFeatures (String[] args){
		processArgs(args);
		extractGffFeatures();
	}
	
	public void extractGffFeatures(){
		//for each file
		for (int i=0; i< gffFiles.length; i++){
			
			//parse file setting types to extract
			Gff3Parser parser = new Gff3Parser();
			parser.setRegExTypes("5'-UTR|3'-UTR");
			parser.setRelax(true);
			if (parser.parseIt(gffFiles[i]) == false) Misc.printExit("Problem parsing gff! "+gffFiles[i]);
			//subtract one from coordinates to put into zero base
			//parser.subtractOneFromFeatures();
			
			Gff3Feature[] gffs = parser.getFeatures();
			System.out.println("num "+gffs.length);
			for (int j=0; j< gffs.length; j++) System.out.println(gffs[j]);
			System.exit(0);
			
			/*HashSet uni = new HashSet();
			for (int x=0; x< gffs.length; x++){
				String[] t = gffs[x].getAttributes().split("\\s");
				if (t.length > 1) {
					System.out.println(t[0]);
					if (uni.contains(t[0])) System.out.println("*************************** Dup "+t[0]);
					else uni.add(t[0]);
				}
			}
			System.exit(0);
			*/
			
			//group by gene text
			LinkedHashMap map = new LinkedHashMap();
			//Pattern text = Pattern.compile(".*RNA\\s([SP|chrM].+)");
			ArrayList al = new ArrayList();
			for (int x=0; x< gffs.length; x++){
				//String nameAtts = gffs[x].getAttributes().split(";")[0].trim();
				//Matcher mat = text.matcher(nameAtts);
				//if (mat.matches() == false) {
				//	System.out.println("Error: matching "+gffs[x].getAttributes()+" with "+text.pattern()+ " skipping!");
				//	continue;
				//}
				//String geneName = mat.group(1);
				//gffs[x].setName(geneName);
				String geneName = gffs[x].getAttributes().split("\\s")[0];
				gffs[x].setName(geneName);
				//add to map
				if (map.containsKey(geneName)) al = (ArrayList) map.get(geneName);
				else {
					al = new ArrayList();
					map.put(geneName, al);
				}
				al.add(gffs[x]);
			}
			
			//print each cluster
			//printClusters(map);
			
			//print in ucsc table format
			printUCSCFormat (map);



		}

	}
	
	/**Prints all the gff lines.*/
	public void printClusters(LinkedHashMap map){
		Iterator it = map.keySet().iterator();
		while (it.hasNext()){
			ArrayList gffAL = (ArrayList) map.get(it.next());
			int num = gffAL.size();
			for (int k=0; k< num; k++){
				Gff3Feature gffFeature = (Gff3Feature) gffAL.get(k);
				System.out.println(gffFeature.toStringNoAttributes()+"\tName="+gffFeature.getName());
			}
			
		}
	}
	
	/**	Print each cluster in ucsc table format
		refGene.name	refGene.chrom	refGene.strand	refGene.txStart	refGene.txEnd	refGene.cdsStart	refGene.cdsEnd	refGene.exonCount	refGene.exonStarts	refGene.exonEnds	refLink.product	geneName.name	refSeqSummary.summary
		 NM_198576	chr1	+	995569	1031415	995619	1030284	36	995569,997647,1010723,1016111,1016522,1016827,1017258,1018541,1018840,1019125,1019411,1019636,1020463,1020661,1021035,1021266,1021462,1021699,1022122,1022629,1022875,1023078,1023314,1024169,1024538,1024868,1025205,1025535,1025729,1026028,1026555,1026755,1027030,1029055,1029750,1030126,	995820,997909,1010771,1016327,1016747,1017052,1017465,1018760,1019035,1019326,1019560,1019742,1020580,1020826,1021179,1021391,1021568,1022038,1022260,1022757,1022990,1023198,1023668,1024362,1024754,1025098,1025340,1025632,1025894,1026140,1026672,1026948,1027118,1029280,1029854,1031415,	agrin	AGRIN	Agrin is a neuronal aggregating factor that induces the aggregation of ...
	*/
	public void printUCSCFormat(LinkedHashMap map){
		Iterator it = map.keySet().iterator();
		while (it.hasNext()){
			ArrayList gffAL = (ArrayList) map.get(it.next());
			int num = gffAL.size();
			//collect exon starts and ends
			int[] starts = new int[num];
			int[] ends = new int[num];
			Gff3Feature gffFeature = null;
			int smallest = 1000000000;
			int largest = -1;
			for (int k=0; k< num; k++){
				gffFeature = (Gff3Feature) gffAL.get(k);
				starts[k] = gffFeature.getStart();
				if (starts[k] < smallest) smallest = starts[k];
//adding 1 to stop of each feature to get to interbase
				ends[k] = gffFeature.getEnd()+1;
				if (ends[k] > largest) largest = ends[k];
			}
			System.out.println(
					gffFeature.getName()+"\t"+ 
					gffFeature.getSeqId()+"\t"+ 
					gffFeature.getStrand()+"\t"+ 
					smallest+"\t"+ 
					largest+"\t"+
					smallest+"\t"+ 
					largest+"\t"+
					starts.length+"\t"+ 
					Misc.intArrayToString(starts, ",")+"\t"+ 
					Misc.intArrayToString(ends, ","));//+"\t"+
					//"x\tx\tx"
					//);
			
			
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ExtractGffFeatures(args);
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
					case 'g': gffFiles = IO.extractFiles(new File (args[i+1])); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (gffFiles == null || gffFiles.length ==0) Misc.printExit("\nError: cannot find your gff file(s)!\n");
	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Extract Gff Features    May 2007                     **\n" +
				"**************************************************************************************\n" +
				"EGF is a workbench for parsing Gff files.  Grab a java programmer and play.\n\n"+
				
				"Parameters:\n"+
				"-g Full path file text for a gff file or directory containing such.\n\n"+
				
		        "**************************************************************************************\n");		
	}
}
