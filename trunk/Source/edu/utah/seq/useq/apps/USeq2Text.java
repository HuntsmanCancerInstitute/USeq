package edu.utah.seq.useq.apps;
import java.io.*;
import java.util.regex.*;
import java.util.zip.*;
import java.util.*;

import util.gen.IO;

import edu.utah.seq.useq.*;
import edu.utah.seq.useq.data.*;


/**Converts USeq binary archives to minimal native output, wig or bed format.*/
public class USeq2Text {

	//fields
	private File[] useqArchives;
	private boolean printBedFormat = false;
	private boolean printWigFormat = false;
	private boolean skipZeroBlockBedGraphs = true;
	private boolean convertScoresToBedFormat = false;

	public USeq2Text(String[] args){
		processArgs(args);

		try {
			
		System.out.println("Processing:");
		//for each zip archive
		for (int i=0; i< useqArchives.length; i++){

			System.out.println("\t"+useqArchives[i].getName());

			if (printWigFormat){
				//is it stranded
				USeqArchive ua = new USeqArchive(useqArchives[i]);
				if (ua.isStranded()){
					
					File wigFile = new File (useqArchives[i].getParentFile(), USeqUtilities.removeExtension(useqArchives[i].getName())+"Plus.wig");
					print2WigFile(useqArchives[i], wigFile, "+");
					wigFile = new File (useqArchives[i].getParentFile(), USeqUtilities.removeExtension(useqArchives[i].getName())+"Minus.wig");
					print2WigFile(useqArchives[i], wigFile, "-");
				}
				else {
					File wigFile = new File (useqArchives[i].getParentFile(), USeqUtilities.removeExtension(useqArchives[i].getName())+".wig");
					print2WigFile(useqArchives[i], wigFile, null);
				}
			}

			else {
				String extension = ".txt";
				if (printBedFormat) extension = ".bed";
				File txtFile = new File (useqArchives[i].getParentFile(), USeqUtilities.removeExtension(useqArchives[i].getName())+extension);
				print2TextFile(useqArchives[i], txtFile, printBedFormat, convertScoresToBedFormat);
			}
		}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("\nDone!");
	}
	
	public USeq2Text(){};

	@SuppressWarnings("unchecked")
	/**Prints a UseqArchive binary file to a text file in either minimal native format (see comment at beginning of file for a description of the columns) or in bed6 or bed12 format.*/
	public static void print2TextFile (File inputUseqArchive, File outputTextFile, boolean printOuputInBedFormat, boolean fixBedScores){
		try {
			PrintWriter out = new PrintWriter (new FileWriter (outputTextFile));
			ZipFile zf = new ZipFile(inputUseqArchive);
			Enumeration<ZipEntry> e = (Enumeration<ZipEntry>) zf.entries();

			//make an ArchiveInfo object on the first element in the zip archive
			ZipEntry ze = e.nextElement();
			if (ze.getName().equals(ArchiveInfo.ARCHIVE_README_NAME) == false) throw new IOException("The first zip entry -> "+ze.getName()+", is not the "+ArchiveInfo.ARCHIVE_README_NAME+"! Aborting.");
			ArchiveInfo ai = new ArchiveInfo(zf.getInputStream(ze), false);

			//write out ai info as comments
			ai.appendCommentedKeyValues(out);

			//load data slices
			while(e.hasMoreElements()) {
				ze = e.nextElement();
				//make a SliceInfo object
				SliceInfo si = new SliceInfo(ze.getName());
				DataInputStream dis = new DataInputStream( new BufferedInputStream(zf.getInputStream(ze)));
				String extension = si.getBinaryType();

				//print bed format
				if (printOuputInBedFormat){
					//call appropriate maker
					//Position
					if (USeqUtilities.POSITION.matcher(extension).matches()) new PositionData (dis, si).writeBed(out);
					//PositionScore
					else if (USeqUtilities.POSITION_SCORE.matcher(extension).matches()) new PositionScoreData (dis, si).writeBed(out, fixBedScores);
					//PositionText
					else if (USeqUtilities.POSITION_TEXT.matcher(extension).matches()) new PositionTextData (dis, si).writeBed(out);
					//PositionScoreText
					else if (USeqUtilities.POSITION_SCORE_TEXT.matcher(extension).matches()) new PositionScoreTextData (dis, si).writeBed(out, fixBedScores);
					//Region
					else if (USeqUtilities.REGION.matcher(extension).matches()) new RegionData (dis, si).writeBed(out);
					//RegionScore
					else if (USeqUtilities.REGION_SCORE.matcher(extension).matches()) new RegionScoreData (dis, si).writeBed(out, fixBedScores);
					//RegionText
					else if (USeqUtilities.REGION_TEXT.matcher(extension).matches())  new RegionTextData (dis, si).writeBed(out);
					//RegionScoreText
					else if (USeqUtilities.REGION_SCORE_TEXT.matcher(extension).matches()) new RegionScoreTextData (dis, si).writeBed(out, fixBedScores);
					else  throw new IOException("\nFailed to recognize the binary file extension! "+ze.getName());
				}

				//print native minimal format
				else {
					//call appropriate maker
					//Position
					if (USeqUtilities.POSITION.matcher(extension).matches()) new PositionData (dis, si).writeNative(out);
					//PositionScore
					else if (USeqUtilities.POSITION_SCORE.matcher(extension).matches()) new PositionScoreData (dis, si).writeNative(out);
					//PositionText
					else if (USeqUtilities.POSITION_TEXT.matcher(extension).matches()) new PositionTextData (dis, si).writeNative(out);
					//PositionScoreText
					else if (USeqUtilities.POSITION_SCORE_TEXT.matcher(extension).matches()) new PositionScoreTextData (dis, si).writeNative(out);
					//Region
					else if (USeqUtilities.REGION.matcher(extension).matches()) new RegionData (dis, si).writeNative(out);
					//RegionScore
					else if (USeqUtilities.REGION_SCORE.matcher(extension).matches()) new RegionScoreData (dis, si).writeNative(out);
					//RegionText
					else if (USeqUtilities.REGION_TEXT.matcher(extension).matches())  new RegionTextData (dis, si).writeNative(out);
					//RegionScoreText
					else if (USeqUtilities.REGION_SCORE_TEXT.matcher(extension).matches()) new RegionScoreTextData (dis, si).writeNative(out);
					else  throw new IOException("\nFailed to recognize the binary file extension! "+ze.getName());
				}
				dis.close();
			}
			out.close();
		} catch (IOException e) {
			System.err.println("\nError, could not process binary archive!");
			e.printStackTrace();
		}
	}

	public String buildWigHeader(File inputUseqArchive, ArchiveInfo ai) throws IOException{

		//what kind of data?
		//add required header
		StringBuilder header = new StringBuilder();
		String graphType = ai.getValue(ArchiveInfo.GRAPH_STYLE_KEY);
		if (graphType == null) throw new IOException("\nFailed to identify a graph type.\n");
		if (graphType.equals(ArchiveInfo.GRAPH_STYLE_VALUE_HEATMAP) || graphType.equals(ArchiveInfo.GRAPH_STYLE_VALUE_STAIRSTEP)) header.append("track type=bedGraph");
		else header.append("track type=wiggle_0");

		//name
		String name = inputUseqArchive.getName().replaceAll(USeqUtilities.USEQ_EXTENSION_WITH_PERIOD, "");
		header.append(" name=\"");
		header.append(name);
		header.append("\"");

		//description
		header.append(" description=\"");
		String desc = ai.getValue(ArchiveInfo.DESCRIPTION_KEY);
		if (desc != null) {
			header.append(desc);
			header.append(" - ");
		}
		header.append(ai.getVersionedGenome());
		header.append("\"");

		//visibility
		header.append(" visibility=full");
		
		//zero graph
		header.append(" alwaysZero=on");

		//any color?
		String hexColor = ai.getValue(ArchiveInfo.COLOR_KEY);
		if (hexColor!=null) {
			String rgb = USeqUtilities.convertHexadecimal2RGB(hexColor, ",");
			if (rgb == null) throw new IOException("\nFailed to convert the hex color code '"+hexColor+"' to rgb. \n");
			header.append(" color=");
			header.append(rgb);
		}

		//dot graph type?
		if (graphType.equals(ArchiveInfo.GRAPH_STYLE_VALUE_DOT)) header.append(" graphType=points");

		return header.toString();
	}

	@SuppressWarnings("unchecked")
	/**Prints a UseqArchive binary file to a variable step wig file. Set strand to null for all or use + or - .*/
	public void print2WigFile (File inputUseqArchive, File outputTextFile, String strand){
		try {
			
			PrintWriter out = new PrintWriter (new FileWriter (outputTextFile));
			ZipFile zf = new ZipFile(inputUseqArchive);
			Enumeration<ZipEntry> e = (Enumeration<ZipEntry>) zf.entries();

			//make an ArchiveInfo object on the first element in the zip archive
			ZipEntry ze = e.nextElement();
			if (ze.getName().equals(ArchiveInfo.ARCHIVE_README_NAME) == false) throw new IOException("The first zip entry -> "+ze.getName()+", is not the "+ArchiveInfo.ARCHIVE_README_NAME+"! Aborting.");
			ArchiveInfo ai = new ArchiveInfo(zf.getInputStream(ze), false);

			//graph data?
			if (ai.isRegionData()) throw new IOException("\nThis USeq archive looks like it contains region data, not graph data.  Use the native text or bed file output option. \n");

			//write out ai info as comments
			out.println("# This wig file was generated by converting the '"+inputUseqArchive+"' archive with the USeq2Text application.");
			ai.appendCommentedKeyValues(out);

			//print header
			out.println(buildWigHeader(inputUseqArchive, ai));

			//should this be saved as a bedGraph?
			boolean bedGraphFormat = false;
			String graphType = ai.getValue(ArchiveInfo.GRAPH_STYLE_KEY);
			if (graphType.equals(ArchiveInfo.GRAPH_STYLE_VALUE_STAIRSTEP) || graphType.equals(ArchiveInfo.GRAPH_STYLE_VALUE_HEATMAP)) bedGraphFormat = true;
			//write data slices
			if (bedGraphFormat) writeBedGraph(zf, e, out, strand);

			else {
				String chromosome = "";
				
				while(e.hasMoreElements()) {
					ze = e.nextElement();
					//make a SliceInfo object
					SliceInfo si = new SliceInfo(ze.getName());
					DataInputStream dis = new DataInputStream( new BufferedInputStream(zf.getInputStream(ze)));
					String extension = si.getBinaryType();
					//correct strand?
					if (strand != null && strand.equals(si.getStrand()) == false) {
						continue;
					}
					//add new chromosome line?
					if (si.getChromosome().equals(chromosome) == false) {
						chromosome = si.getChromosome();
						out.println("variableStep chrom="+chromosome);

					}
					//call appropriate maker
					//Position
					if (USeqUtilities.POSITION.matcher(extension).matches()) {
						//System.out.println("POSITION");
						new PositionData (dis, si).writePositionScore(out);
					}
					//PositionScore
					else if (USeqUtilities.POSITION_SCORE.matcher(extension).matches()) {
						//System.out.println("POSITION_SCORE");
						new PositionScoreData (dis, si).writePositionScore(out);
					}
					//PositionText
					else if (USeqUtilities.POSITION_TEXT.matcher(extension).matches()) {
						//System.out.println("POSITION_TEXT");
						new PositionTextData (dis, si).writePositionScore(out);
					}
					//PositionScoreText
					else if (USeqUtilities.POSITION_SCORE_TEXT.matcher(extension).matches()) {
						//System.out.println("POSITION_SCORE_TEXT");
						new PositionScoreTextData (dis, si).writePositionScore(out);
					}
					else  throw new IOException("\nThis USeq archive looks like it contains region data, not graph data.  Use the native text or bed file output option. \n");
					dis.close();
				}
			}
			//close the PrintWriter
			out.close();
		} catch (IOException e) {
			System.err.println("\nError, could not process your binary archive!");
			outputTextFile.delete();
			e.printStackTrace();
		}
	}
	public void writeBedGraph(ZipFile zf, Enumeration<ZipEntry> e, PrintWriter out, String strand) throws IOException {
		
		ZipEntry ze;
		String tab = "\t";
		String chromosome = null;
		ArrayList<PositionScoreData> psAL = new ArrayList<PositionScoreData>();
		while(e.hasMoreElements()) {
			ze = e.nextElement();
			//make a SliceInfo object
			SliceInfo si = new SliceInfo(ze.getName());
			DataInputStream dis = new DataInputStream( new BufferedInputStream(zf.getInputStream(ze)));
			String extension = si.getBinaryType();
			//correct strand?
			if (strand != null && strand.equals(si.getStrand()) == false) {
				continue;
				
			}
			//new chromosome line?
			if (chromosome == null) {
				chromosome = si.getChromosome();
			}
			
			else if (si.getChromosome().equals(chromosome) == false) {
				//merge
				PositionScoreData merged = PositionScoreData.merge(psAL);
				int[] positions = merged.getBasePositions();
				float[] scores = merged.getBaseScores();
				
				//write paired blocks that have the same score
				int lastPosition = -1;
				for (int i=0; i< scores.length; i++){
					int next = i+1;
					if (next == scores.length) break;
					//same score?
					if (scores[i] == scores[next]){
						//zero?
						if (skipZeroBlockBedGraphs && scores[i] == 0) continue;
						out.print(chromosome); out.print(tab);
						
						//check position
						int pos = positions[i];
						if (pos < lastPosition) pos = lastPosition;
						else lastPosition = pos;
						out.print(pos); out.print(tab); 
						
						pos = positions[next]+1;
						if (pos < lastPosition) pos = lastPosition;
						else lastPosition = pos;
						
						out.print(pos); out.print(tab); 
						out.println(scores[i]);
						i++;
					}
					//different scores
					else {
						//look to see if next position is one off
						if ((positions[next]-1) == positions[i]){
							//zero?
							if (skipZeroBlockBedGraphs && scores[i] == 0) continue;
							out.print(chromosome); out.print(tab); 
							
							//check position
							int pos = positions[i];
							if (pos < lastPosition) pos = lastPosition;
							else lastPosition = pos;
							out.print(pos); out.print(tab); 
							
							pos = positions[next]+1;
							if (pos < lastPosition) pos = lastPosition;
							else lastPosition = pos;
							
							out.print(pos); out.print(tab); 
							out.println(scores[i]);
						}
					}
				}
				//clear
				psAL.clear();
				chromosome = si.getChromosome();
			}
			
			//fetch PositionScore data and add to ArrayList
			//PositionScore
			if (USeqUtilities.POSITION_SCORE.matcher(extension).matches()) {
				psAL.add(new PositionScoreData (dis, si));
			}
			//PositionScoreText
			else if (USeqUtilities.POSITION_SCORE_TEXT.matcher(extension).matches()) {
				//new PositionScoreTextData (dis, si).writePositionScore(out);
				PositionScoreTextData p = new PositionScoreTextData (dis, si);
				psAL.add(new PositionScoreData (p.getBasePositions(), p.getBaseScores(), si));
			}
			else  throw new IOException("\nThis USeq archive lacks score information thus it cannot be made into a bed graph!  Use the native text or bed file output option. \n");
			
			
			dis.close();
		}
		
		//write last chromosome, might be the first if only one slice!
		//merge
		PositionScoreData merged = PositionScoreData.merge(psAL);
		int[] positions = merged.getBasePositions();
		float[] scores = merged.getBaseScores();

		//write paired blocks that have the same score
		int lastPosition = -1;
		for (int i=0; i< scores.length; i++){
			int next = i+1;
			if (next == scores.length) break;
			//same score?
			if (scores[i] == scores[next]){
				//zero?
				if (skipZeroBlockBedGraphs && scores[i] == 0) continue;
				out.print(chromosome); out.print(tab); 
				//check position
				int pos = positions[i];
				if (pos < lastPosition) pos = lastPosition;
				else lastPosition = pos;
				out.print(pos); out.print(tab); 
				
				
				pos = positions[next]+1;
				if (pos < lastPosition) pos = lastPosition;
				else lastPosition = pos;
				
				out.print(pos); out.print(tab); 
				out.println(scores[i]);
			} else {
				//look to see if next position is one off
				if ((positions[next]-1) == positions[i]){
					//zero?
					if (skipZeroBlockBedGraphs && scores[i] == 0) continue;
					out.print(chromosome); out.print(tab); 
					
					//check position
					int pos = positions[i];
					if (pos < lastPosition) pos = lastPosition;
					else lastPosition = pos;
					out.print(pos); out.print(tab); 
					
					pos = positions[next]+1;
					if (pos < lastPosition) pos = lastPosition;
					else lastPosition = pos;
					
					out.print(pos); out.print(tab); 
					out.println(scores[i]);
				}
			}
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new USeq2Text(args);
	}

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+USeqUtilities.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': useqArchives = USeqUtilities.extractFiles(new File(args[++i]), USeqUtilities.USEQ_EXTENSION_NO_PERIOD); break;
					case 'b': printBedFormat = true; break;
					case 'c': convertScoresToBedFormat = true; break;
					case 'w': printWigFormat = true; printBedFormat = false; break;
					case 'h': printDocs(); System.exit(0); break;
					default: USeqUtilities.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					USeqUtilities.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//pull files
		if (useqArchives == null || useqArchives.length == 0) USeqUtilities.printExit("\nCannot find any xxx."+USeqUtilities.USEQ_EXTENSION_NO_PERIOD+" USeq archives?\n");

	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                USeq 2 Text: Oct 2012                             **\n" +
				"**************************************************************************************\n" +
				"Converts USeq archives to text either as minimal native, bed, or wig graph format. \n" +
				"\n" +

				"\nOptions:\n"+
				"-f Full path file/directory containing xxx."+USeqUtilities.USEQ_EXTENSION_NO_PERIOD+" files.\n" +
				"-b Print bed format, defaults to native text format.\n"+
				"-c Convert scores to bed format 0-1000.\n"+
				"-w Print wig graph format (var step or bed graph), defaults to native format.\n\n"+

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/USeq2Text -f\n" +
				"      /AnalysisResults/USeqDataArchives/ \n\n" +

		"**************************************************************************************\n");

	}

	public File[] getUseqArchives() {
		return useqArchives;
	}

	public void setUseqArchives(File[] useqArchives) {
		this.useqArchives = useqArchives;
	}

	public boolean isPrintBedFormat() {
		return printBedFormat;
	}

	public void setPrintBedFormat(boolean printBedFormat) {
		this.printBedFormat = printBedFormat;
	}

	public boolean isPrintWigFormat() {
		return printWigFormat;
	}

	public void setPrintWigFormat(boolean printWigFormat) {
		this.printWigFormat = printWigFormat;
	}

	public boolean isSkipZeroBlockBedGraphs() {
		return skipZeroBlockBedGraphs;
	}

	public void setSkipZeroBlockBedGraphs(boolean skipZeroBlockBedGraphs) {
		this.skipZeroBlockBedGraphs = skipZeroBlockBedGraphs;
	}


}
