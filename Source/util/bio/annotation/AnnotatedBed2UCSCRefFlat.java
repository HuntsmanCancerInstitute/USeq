package util.bio.annotation;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;
import util.bio.parsers.*;

public class AnnotatedBed2UCSCRefFlat {

	private File bedFile;
	private HashMap<String, UCSCGeneLine[]> chromGenes = null;
	private int minNumExons = 1;
	private String requiredAnnoType = null;
	
	public AnnotatedBed2UCSCRefFlat (String[] args) throws FileNotFoundException, IOException{
		processArgs(args);
		parseBed();
		saveUCSC();
	}


	private void saveUCSC() throws FileNotFoundException, IOException {
		//write out the ucsc file
		String name = Misc.removeExtension(bedFile.getName());
		File ucscFile = new File (bedFile.getParentFile(), name+".ucsc.gz");
		Gzipper out = new Gzipper( ucscFile);
		for (String chr: chromGenes.keySet()) {
			for (UCSCGeneLine ugl: chromGenes.get(chr)) {
				out.println(ugl.toUCSC());
			}
		}
		out.close();
		IO.pl("\nSaved "+ucscFile);
		
	}


	private void parseBed() {
		//split bed by geneName, assuming the name looks like "NSG00000124260.7_MAGEA10_3pUTR,ENSG00000266560.1_RP11-1007I13.4_Intron"
		Bed[] bed = Bed.parseFile(bedFile, 0, 0);
		HashMap<String, ArrayList<Bed>> geneBed = new HashMap<String, ArrayList<Bed>>();
		for (Bed b: bed){
			if (b.getName() == null || b.getName().equals(".")) Misc.printErrAndExit("ERROR; this bed line doesn't contain a gene name "+b.toString());
			String[] splitName = Misc.COMMA.split(b.getName());
			for (String sn: splitName) {
				//appropriate anno type
				if (requiredAnnoType !=null && sn.contains(requiredAnnoType)==false) continue;
				int lastUnderscore = sn.lastIndexOf('_');
				String geneName = sn.substring(0, lastUnderscore);
				ArrayList<Bed> beds = geneBed.get(geneName);
				if (beds == null) {
					beds = new ArrayList<Bed>();
					geneBed.put(geneName, beds);
				}
				beds.add(b);
			}
		}
		IO.pl("\t"+geneBed.size()+ "\tNumber of Genes ");

		//build genes for those with two or more exons
		ArrayList<UCSCGeneLine> genesAL = new ArrayList<UCSCGeneLine>();
		for (String name: geneBed.keySet()) {
			ArrayList<Bed> beds = geneBed.get(name);
			if (beds.size() < minNumExons) continue;
			Bed[] toSort = new Bed[beds.size()];
			beds.toArray(toSort);
			Arrays.sort(toSort);
			UCSCGeneLine ugl = new UCSCGeneLine(toSort, name);
			genesAL.add(ugl);
		}
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader();
		UCSCGeneLine[] genes = new UCSCGeneLine[genesAL.size()];
		genesAL.toArray(genes);
		reader.setGeneLines(genes);

		IO.pl("\t"+ genes.length + "\tPassing criteria");

		if (genes == null || genes.length == 0) Misc.printExit("\nProblem loading your USCS gene model table or bed file? No genes/ regions?\n");
		//check ordering
		if (reader.checkStartStopOrder() == false) Misc.printExit("\nOne of your regions's coordinates are reversed. Check that each start is less than the stop.\n");
		//check gene name is unique
		if (reader.uniqueGeneNames() == false) Misc.printExit("\nDuplicate gene names were found in your gene / bed file, these must be unique.\n");
		//check that genes are stranded
		if (reader.checkStrand() == false) Misc.printExit("\nError: your bed file doesn't appear to be stranded?\n");
		chromGenes = reader.getChromSpecificGeneLines();

	}


	public static void main(String[] args) throws FileNotFoundException, IOException {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AnnotatedBed2UCSCRefFlat(args); 
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bedFile = new File (args[i+1]); i++; break;
					case 'm': minNumExons = Integer.parseInt(args[i+1]); i++; break;
					case 'a': requiredAnnoType = args[i+1]; i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (bedFile == null || bedFile.canRead() == false) Misc.printExit("\nError: cannot find your bed file!\n");	
	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Annotated Bed 2 UCSC RefFlat   June 2015                **\n" +
				"**************************************************************************************\n" +
				"Run AnnotateBedWithGenes using the addGeneFeature option then use this tool to convert\n"+
				"the annotated bed to genes for those with a gene name. To convert to GTF, use the UCSC\n"+
				"genePredToGtf tool.\n\n"+

				"Parameters:\n"+
				"-b Annotated bed file to convert\n"+
				"-m Minimum number exons for conversion, defaults to 1\n"+
				"-a Required anno type, defaults to none\n"+

				"\nExample: java -Xmx2G -jar pathTo/USeq/Apps/AnnotatedBed2UCSCRefFlat -b anno.bed.gz\n"+
				"   -a 3pUTR -m 2\n\n"+

				"**************************************************************************************\n");		
	}		



}

