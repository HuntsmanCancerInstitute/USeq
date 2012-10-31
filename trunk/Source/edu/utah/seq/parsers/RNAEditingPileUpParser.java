package edu.utah.seq.parsers;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.data.ComparatorPointPosition;
import edu.utah.seq.data.Info;
import edu.utah.seq.data.Point;
import edu.utah.seq.data.PointData;
import util.gen.*;


/**@author davidnix*/
public class RNAEditingPileUpParser {

	//user defined fields
	private File[] pileupFiles;
	private float minimumReadCoverage = 5.0f;
	private float minimumFractionEditedBase = 0.01f;
	private String versionedGenome;
	
	//internal fields
	private Pattern space = Pattern.compile("\\t");
	private Pattern GBase = Pattern.compile("G");
	private Pattern cBase = Pattern.compile("c");
	private int chromIndex = 0;
	private int positionIndex = 1;
	private int refseqIndex = 2;
	private int readDepthIndex = 3;
	private int baseCallIndex = 4;
	//private int baseScoreIndex = 5;
	private File dataDirectory;
	private ArrayList<File> dataFiles;

	public RNAEditingPileUpParser(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		
		System.out.println("Processing pileup files:");
		
		for (File pf : pileupFiles){
			System.out.println("\t" + pf.getName());
			
			parseFile(pf);
			
			convertToPointData();
		}
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}
	
	public void convertToPointData(){
		// for each data file
		for (File dataFile: dataFiles){
			String fileName = dataFile.getName();
			int fileLength = fileName.length();
			
			String strand = fileName.substring(fileLength-1);
			String chromosome = fileName.substring(0, fileLength-1);
			
			//parse binary data: position fraction
			DataInputStream dis = null;
			ArrayList<Point> ptAL = new ArrayList<Point>();
			
			try {
				dis = new DataInputStream(new BufferedInputStream(new FileInputStream(dataFile)));
				
				while (true){
					int position = dis.readInt();
					float fraction = dis.readFloat();
					ptAL.add(new Point (position, fraction));
				}
			} catch (EOFException eof){	
				if (ptAL.size() != 0){
					PointData pd = Point.extractPositionScores(ptAL);
					Info info = new Info("ParsedPileUp", versionedGenome, chromosome, strand, 1, null);
					pd.setInfo(info);
					//save
					pd.writePointData(dataDirectory);
				}
			}
			catch (Exception e){
				e.printStackTrace();
				Misc.printErrAndExit("\nError encountered in parsing this binary file? "+dataFile);
			} finally {
				if (dis != null) {
					try {
						dis.close();
					} catch (IOException ignored) {}
				}
			}
		}
	}
	
	public void parseFile(File pileupFile){
		try {
			//input streams
			BufferedReader in = IO.fetchBufferedReader(pileupFile);
			String baseName = Misc.removeExtension(pileupFile.getCanonicalPath());
			dataDirectory = new File (baseName+ "_" +(int)minimumReadCoverage +"RC"+minimumFractionEditedBase+"FE");
			dataDirectory.mkdir();
			dataFiles = new ArrayList<File>();
			
			//data output streams
			String chromosome = "";
			DataOutputStream dosPlus = null;
			DataOutputStream dosMinus = null;
			
			//for each line in the file
			String line;
			System.out.print("\t\t");
			while ((line = in.readLine()) != null){
				String[] tokens = space.split(line);
				if (tokens.length < 5){
					System.err.println("\nMalformed pileup line, skipping -> "+line+"\n\t"+tokens.length+" Tokens");
					continue;
				}
				//refseq base an A?
				String refSeqBase = tokens[refseqIndex].toUpperCase();
				
				//minimum coverage
				float coverage = Float.parseFloat(tokens[readDepthIndex]);
				if (coverage < minimumReadCoverage) continue;
				
				//new chromosome?
				if (tokens[chromIndex].equals(chromosome) == false){
					chromosome = tokens[chromIndex];
					System.out.print(chromosome+" ");
					//close old streams
					if (dosPlus != null) {
						dosPlus.close();
						dosMinus.close();
					}
					//make new ones
					File f = new File (dataDirectory, chromosome+"+");
					dosPlus = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(f)));
					dataFiles.add(f);
					f.deleteOnExit();
					f = new File (dataDirectory, chromosome+"-");
					dosMinus = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(f)));
					dataFiles.add(f);
					f.deleteOnExit();
				}
				
				//plus strand
				if (refSeqBase.equals("A")){
					//look for upper case G's
					int numGs = countGs(tokens[baseCallIndex]);
					if (numGs !=0){
						float fraction = ((float)numGs) / coverage;
						if (fraction < minimumFractionEditedBase) continue;
						int position = Integer.parseInt(tokens[positionIndex]) -1;
						dosPlus.writeInt(position);
						dosPlus.writeFloat(fraction);
					}
				}
				//minus strand
				else if (refSeqBase.equals("T")){
					//look for lower case c's
					int numCs = countCs(tokens[baseCallIndex]);
					if (numCs !=0){
						float fraction = ((float)numCs) / Float.parseFloat(tokens[readDepthIndex]);
						if (fraction < minimumFractionEditedBase) continue;
						int position = Integer.parseInt(tokens[positionIndex]) -1;
						dosMinus.writeInt(position);
						dosMinus.writeFloat(fraction);
					}
				}
				
			}
			System.out.println();
			
			//close streams
			if (dosPlus != null) {
				dosPlus.close();
				dosMinus.close();
			}
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public int countGs(String baseCalls){
		Matcher mat = GBase.matcher(baseCalls);
		int num = 0;
		while (mat.find()) num++;
		return num;
	}
	
	public int countCs(String baseCalls){
		Matcher mat = cBase.matcher(baseCalls);
		int num = 0;
		while (mat.find()) num++;
		return num;
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new RNAEditingPileUpParser(args);
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
					case 'p': pileupFiles = IO.extractFiles( new File(args[++i])); break;
					case 'r': minimumReadCoverage = Float.parseFloat(args[++i]); break;
					case 'f': minimumFractionEditedBase = Float.parseFloat(args[++i]); break;
					case 'v': versionedGenome = args[i+1]; i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for bam files
		if (pileupFiles == null || pileupFiles.length == 0 || pileupFiles[0].canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your pileup file(s)?\n");
		if (versionedGenome == null) Misc.printExit("\nPlease enter a genome version recognized by UCSC, see http://genome.ucsc.edu/FAQ/FAQreleases.\n");
		
		System.out.println(minimumReadCoverage+"\tMinimum Read Coverage");
		System.out.println(minimumFractionEditedBase+"\tMinimum Fraction Edited Bases\n");
		
		
	}	
	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            RNA Editing PileUp Parser: Jan 2012                       **\n" +
				"**************************************************************************************\n" +
				"Parses a SAMTools mpileup output file for refseq A bases that show evidence of\n" +
				"RNA editing via conversion to Gs, stranded. Bases passing the thresholds are\n" +
				"written to file in bar format for viewing in IGB and subsequent clustering with\n" +
				"the RNAEditingScanSeqs app.\n\n"+

				"Options:\n"+
				"-p Path to pileup file (.gz or.zip OK) or directory containing such. Multiple files\n" +
				"      are processed independently.\n"+
				"-v Versioned Genome (ie H_sapiens_Mar_2006), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +
				"-r Minimum read coverage, defaults to 5.\n"+
				"-f Minimum fraction edited base, defaults to 0.01\n"+

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/RNAEditingPileUpParser -p /Pileups/\n" +
				"      -v C_elegans_May_2008\n\n" +

		"**************************************************************************************\n");

	}

}
