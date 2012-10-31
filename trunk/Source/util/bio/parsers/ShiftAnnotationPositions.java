package util.bio.parsers;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.*;
import util.bio.annotation.*;
import util.bio.seq.ConcatinateFastas;

/**Shifts the bp position of UCSCGeneLines matching particular chromosomes. Good for shifting coordinates after making a composite mock chromosome from 
 * unassembled sequences/ scaffolds.*/
public class ShiftAnnotationPositions {

	private File ucscGeneFile;
	private File shifterFile;
	private String chromName;
	private File shiftedGeneFile;
	private UCSCGeneModelTableReader reader;
	UCSCGeneLine[] lines;
	private File bedFile;


	public ShiftAnnotationPositions(String[] args){

		processArgs(args);

		loadGeneModels();

		HashMap<String, Integer> nameNumber = IO.loadFileIntoHashMapStringInteger(shifterFile);
		int numShifted = 0;

		try {
			PrintWriter out = new PrintWriter (new FileWriter (shiftedGeneFile));
			for (int i=0; i< lines.length; i++){
				String testChrom = lines[i].getChrom();
				//shift it?
				if (nameNumber.containsKey(testChrom)){
					numShifted++;
					lines[i].setChrom(chromName);
					int add = nameNumber.get(testChrom).intValue();
					lines[i].setTxStart(lines[i].getTxStart() + add);
					lines[i].setTxEnd(lines[i].getTxEnd() + add);
					if (ucscGeneFile != null){
						lines[i].setTss(lines[i].getTss() + add);
						lines[i].setCdsStart(lines[i].getCdsStart() + add);
						lines[i].setCdsEnd(lines[i].getCdsEnd() + add);
						ExonIntron[] ex = lines[i].getExons();
						for (int x=0; x< ex.length; x++){
							ex[x].setStart(ex[x].getStart() +add);
							ex[x].setEnd(ex[x].getEnd() +add);
						}
					}
				}
				if (bedFile == null) out.println(lines[i]);
				else out.println(lines[i].toStringBedFormat());
			}
			out.close();
			System.out.print("Shifted "+numShifted+" annotations. Each assigned to chromosome '"+chromName+"'.\n\n");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public void loadGeneModels(){
		//load gene models from refFlat for refSeq UCSC gene table
		UCSCGeneModelTableReader reader = null;
		if (ucscGeneFile != null) {
			reader = new UCSCGeneModelTableReader(ucscGeneFile, 0);
			lines = reader.getGeneLines();
		}

		//or from bed file
		else if (bedFile != null) {
			Bed[] bed = Bed.parseFile(bedFile, 0, 0);
			if (bed != null) {
				lines = new UCSCGeneLine[bed.length];
				for (int i=0; i< bed.length; i++){
					lines[i] = new UCSCGeneLine(bed[i]);
				}
				reader = new UCSCGeneModelTableReader();
				reader.setGeneLines(lines);
			}
		}
		if (lines == null || lines.length == 0) Misc.printExit("\nProblem loading your USCS gene model table or bed file, no annotations parsed?\n");
		//check ordering
		if (reader.checkStartStopOrder() == false) Misc.printExit("\nOne of your regions's coordinates are reversed. Check that each start is less than the stop.\n");
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ShiftAnnotationPositions(args);
	}		

	/**This method will process each argument and assign new variables*/
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
					case 'u': ucscGeneFile = new File(args[++i]); break;
					case 's': shifterFile = new File(args[++i]); break;
					case 'b': bedFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for files
		if (bedFile == null && ucscGeneFile == null) Misc.printErrAndExit("\nError: please enter either a UCSC or bed annotation file to shift.\n");
		if (bedFile != null && ucscGeneFile != null) Misc.printErrAndExit("\nError: please enter either a UCSC or bed annotation file to shift BUT not both!\n");
		if (ucscGeneFile != null && ucscGeneFile.exists() == false) Misc.printErrAndExit("\nError: cannot find your ucsc annotation file?\n");
		if (bedFile != null && bedFile.exists()== false) Misc.printErrAndExit("\nError: cannot find your bed annotation file?\n");
		if (shifterFile != null && shifterFile.exists()== false) Misc.printErrAndExit("\nError: cannot find your xxx.shifter file?\n");
		if (shifterFile.getName().endsWith(".shifter.txt") == false) Misc.printErrAndExit("\nError: your -s xxx.shifter.txt file does not end with .shifter.txt?!\n");

		//set chromosome name
		chromName = shifterFile.getName().replace(".shifter.txt", "");

		String name;
		File parentFile;
		if (ucscGeneFile != null) {
			name = ucscGeneFile.getName();
			parentFile = ucscGeneFile.getParentFile();
		}
		else {
			name = bedFile.getName();
			parentFile = bedFile.getParentFile();
		}
		name = name.replace(".gz", "");
		name = name.replace(".txt", "");
		name = name.replace(".bed", "");
		String extension = ".txt";
		if (bedFile != null) extension = ".bed";
		shiftedGeneFile = new File (parentFile, name+"_"+chromName+"Shifted"+extension);

	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         Shift Annotation Positions: Oct 2010                     **\n" +
				"**************************************************************************************\n" +
				"Uses the information in an xxx.shifter.txt file from the ConcatinateFastas app to\n" +
				"shift the annotation to match the coordinates of the concatinated sequence. Good for\n" +
				"working with poorly assembled genomes. Run this multiple times with different shifter\n" +
				"files. All files are assumed to use interbase coordinates.\n\n" +

				"Options:\n"+
				"-b Full path file name for a xxx.bed formatted annotation file.\n"+
				"-u (OR) Full path file name for a UCSC refflat/ refseq formatted gene table.\n"+
				"-s Full path file name for the xxx.shifter.txt file from the ConcatinateFastas app.\n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/ShiftAnnotationPositions \n" +
				"    -u /zv8/ucscRefSeq.txt -f /zv8/BadFastas/chrScaffold.shifter.txt\n\n" +

		"**************************************************************************************\n");

	}

}
