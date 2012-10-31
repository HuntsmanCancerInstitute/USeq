package trans.misc;
import util.gen.*;
import edu.utah.seq.parsers.*;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import trans.roc.*;

public class Wig2Bar {

	private File[] files;
	private String genomeVersion;
	//don't change from -66666666666.0f
	private float skipValue = -66666666666.0f;

	public Wig2Bar(String[] args) {
		//check for args 
		processArgs(args);

		//for each file
		try{
			System.out.println("\nParsing...");
			for (int x=0; x< files.length; x++){
				System.out.println("\t"+files[x].getName());
				//what kind of wig file?
				String type =  parseWigFileType(files[x]);
				if (type == null) Misc.printExit("\nCould not parse the wig file type from this file, aborting. -> "+files[x]);
				//parse file
				GrGraph[] grs = null;
				if (type.equals("variableStep")) grs = parseVariableStepWigFile(files[x], skipValue);
				else if (type.equals("fixedStep")) grs = parseFixedStepWigFile(files[x], skipValue);
				else Misc.printExit("\nCannot parse this wig file type! -> "+type);

				//make a save directory
				String fileName = files[x].getName();
				fileName = fileName.replace(".gz", "");
				fileName = fileName.replace(".zip", "");
				fileName = Misc.removeExtension(fileName);
				String firstLetter = fileName.substring(0, 1).toUpperCase();
				fileName = firstLetter + fileName.substring(1);
				File barDir = new File (files[x].getParentFile(), fileName);
				barDir.mkdir();

				//for each GrGraph, write bar file
				BarParser bp = new BarParser();
				HashMap<String,String> tagValues = new HashMap<String,String>();
				tagValues.put(BarParser.SOURCE_TAG, type+" wig file "+files[x].getName());
				tagValues.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
				for (int i=0; i< grs.length; i++){
					File bar = new File (barDir, grs[i].getChromosome()+".bar");
					bp.writeBarFile(bar, grs[i].getChromosome(), genomeVersion, '.', grs[i].getBasePositions(), grs[i].getValues(), tagValues);
				}
			}
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	/**Looks for the first line that begins with 'fixedStep' or 'variableStep', if none found returns null, otherwise returns one of those
	 * two names.*/
	public static String parseWigFileType(File wigFile) throws Exception{
		BufferedReader in = IO.fetchBufferedReader(wigFile);
		String line;
		String type = null;
		while ((line=in.readLine())!=null){
			line = line.trim();
			if (line.length()==0) continue;
			if (line.startsWith("fixedStep")) {
				type = "fixedStep";
				break;
			}
			if (line.startsWith("variableStep")) {
				type = "variableStep";
				break;
			}
		}
		in.close();
		return type;
	}

	/**Parses a variableStep wig file, skips any track type data. Span is assumed to be 1.
	 * Should look something like:
	 * 	variableStep chrom=chr19 span=1
	 * 	59304701 10.0
	 * 	59304901 12.5
	 * 	59305401 15.0
	 * 	59305601 17.5
	 * Subtracts 1 from each position to convert from 1-relative coordinates specified from UCSC Wig fromat
	 * to interbase coordinates. 
	 * Skips lines with designated skipValue*/
	public static GrGraph[] parseVariableStepWigFile(File wigFile, float skipValue) throws Exception{
		ArrayList<GrGraph> graphs = new ArrayList<GrGraph>();
		//load file
		String line;
		String[] tokens = null;
		String chromosome = null;
		ArrayList<Gr> grs = new ArrayList<Gr>();
		BufferedReader in = IO.fetchBufferedReader(wigFile);
		Pattern number = Pattern.compile("^\\d");
		Pattern space = Pattern.compile("\\s");
		Pattern equal = Pattern.compile("=");
		Matcher mat;
		//for each line
		while ((line=in.readLine())!=null){
			line = line.trim();
			if (line.length()==0) continue;
			//start with a number?
			mat = number.matcher(line);
			if (mat.find()){
				tokens = space.split(line);
				if (tokens.length !=2) throw new Exception("Problem with parsing position:value from "+wigFile+" line -> "+line);
				float value = Float.parseFloat(tokens[1]);
				if (value != -66666666666.0f) grs.add(new Gr(Integer.parseInt(tokens[0])-1, value));
			}
			//variableStep
			else if (line.startsWith("variableStep")){
				//parse chrom
				tokens = space.split(line);
				tokens = equal.split(tokens[1]);
				if (tokens.length !=2) throw new Exception ("Problem parsing chromosome from"+wigFile+" line -> "+line); 
				//first one or old
				if (chromosome != null) {
					Gr[] g = new Gr[grs.size()];
					grs.toArray(g);
					graphs.add(new GrGraph(chromosome, null, g));
					grs.clear();
					g = null;
				}
				chromosome = tokens[1];
			}
		}

		if (chromosome == null) throw new Exception ("No 'variableStep chrom=...' line found in "+wigFile);
		//save last chromosome
		Gr[] g = new Gr[grs.size()];
		grs.toArray(g);
		graphs.add(new GrGraph(chromosome, null, g));
		grs.clear();
		g = null;
		in.close();

		GrGraph[] gg = new GrGraph[graphs.size()];
		graphs.toArray(gg);
		return gg;
	}

	/**Parses a fixedStep wig file, skips any track type data. 
	 * Should look something like:
	 * 	fixedStep chrom=chrY start=668 step=1
	 *	0.012
	 *	0.021
	 *	0.028
	 *	0.033
	 *	0.036
	 * Subtracts 1 from each position to convert from 1-relative coordinates specified from UCSC Wig format
	 * to interbase coordinates. */
	public static GrGraph[] parseFixedStepWigFile(File wigFile, float skipValue) throws Exception{
		ArrayList<GrGraph> graphs = new ArrayList<GrGraph>();
		//load file
		String line;
		String[] tokens = null;
		String chromosome = null;
		ArrayList<Gr> grs = new ArrayList<Gr>();
		BufferedReader in = IO.fetchBufferedReader(wigFile);

		Pattern number = Pattern.compile("^\\d");
		Pattern space = Pattern.compile("\\s");
		Pattern equal = Pattern.compile("=");
		Matcher mat;
		HashSet<String> chroms = new HashSet<String>();
		
		int startPosition = 0;
		int stepSize = 0;
		
		//for each line
		while ((line=in.readLine())!=null){
			line = line.trim();
			//empty?
			if (line.length() == 0) continue;
			//start with a number?
			mat = number.matcher(line);
			if (mat.find()){
				float value = Float.parseFloat(line);
				if (value != skipValue) {
					grs.add(new Gr(startPosition, value));
					//System.out.println("\tGR "+chromosome+" "+startPosition+" "+value);	
				}
				//else System.out.println("Skipping "+value);
				//increment position
				startPosition+= stepSize;
			}
			//fixedStep?
			else if (line.startsWith("fixedStep")){
				//split line and check 'fixedStep chrom=chrY start=668 step=1'
				tokens = space.split(line);
				if (tokens.length !=4) throw new Exception("Problem with parsing fixedStep line from "+wigFile+" line -> "+line);
				//parse chrom
				String[] chromTokens = equal.split(tokens[1]);
				if (chromTokens.length !=2) throw new Exception ("Problem parsing chromosome from"+wigFile+" line -> "+line); 
				//first one or old
				if (chromosome == null) chromosome = chromTokens[1];
				if (chromosome != null) {
					//different chromosome?
					if (chromosome.equals(chromTokens[1]) == false){
						//close old and start new
						Gr[] g = new Gr[grs.size()];
						grs.toArray(g);
						if (chroms.contains(chromosome) == false) chroms.add(chromosome);
						else Misc.printExit("\nWig file is not sorted by chromosome! Aborting.\n");
						graphs.add(new GrGraph(chromosome, null, g));
						grs.clear();
						chromosome = chromTokens[1];
					}
				}
				//set start
				String[] startTokens = equal.split(tokens[2]);
				if (startTokens.length !=2) throw new Exception ("Problem parsing start position from"+wigFile+" line -> "+line);
				startPosition = Integer.parseInt(startTokens[1]) -1;
				//set step
				String[] stepTokens = equal.split(tokens[3]);
				if (stepTokens.length !=2) throw new Exception ("Problem parsing start position from"+wigFile+" line -> "+line);
				stepSize = Integer.parseInt(stepTokens[1]);
				//System.out.println(chromosome+" "+startPosition+" "+stepSize);
			}
		}

		if (chromosome == null) throw new Exception ("No 'fixedStep chrom=...' line found in "+wigFile);
		//save last chromosome
		Gr[] g = new Gr[grs.size()];
		grs.toArray(g);
		graphs.add(new GrGraph(chromosome, null, g));
		grs.clear();
		g = null;
		in.close();

		GrGraph[] gg = new GrGraph[graphs.size()];
		graphs.toArray(gg);
		return gg;
	}

	public static void main(String[] args) {
		if (args.length==0){
			printDocs();
			System.exit(0);
		}	
		new Wig2Bar(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File file = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': file = new File(args[i+1]); i++; break;
					case 'v': genomeVersion = args[i+1]; i++; break;
					case 's': skipValue = Float.parseFloat(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nError: unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//parse files and genome version
		if (file == null) Misc.printExit("\nError: cannot find your xxx.wig file(s)?");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(file,".wig");
		tot[1] = IO.extractFiles(file,".wig.zip");
		tot[2] = IO.extractFiles(file,".wig.gz");

		files = IO.collapseFileArray(tot);
		if (files == null || files.length == 0) Misc.printExit("\nError: cannot find your xxx.wig file(s)?");
		if (genomeVersion == null) Misc.printExit("\nError: you must supply a genome version. Goto http://genome.ucsc.edu/cgi-" +
		"bin/hgGateway load your organism to find the associated genome version.\n");
	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Wig2Bar: Oct 2009                                 **\n" +
				"**************************************************************************************\n" +
				"Converts variable step and fixed step xxx.wig(.zip/.gz OK) files to chrom specific\n" +
				"bar files.\n\n" +

				"-f The full path directory/file text for your xxx.wig(.gz/.zip OK) file(s).\n" +
				"-v Genome version (ie H_sapiens_Mar_2006), get from UCSC Browser,\n" +
				"      http://genome.ucsc.edu/FAQ/FAQreleases\n"+ 
				"-s Skip wig lines with designated value/score.\n"+

				"\nExample: java -Xmx1500M -jar pathTo/Apps/Wig2Bar -f /WigFiles/ -v hg18 -s 0.0 \n\n" +

		"**************************************************************************************\n");		
	}

}
