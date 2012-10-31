package edu.utah.seq.useq.data;
import edu.utah.seq.useq.*;
import edu.utah.seq.useq.apps.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

/**Container for a sorted PositionScoreText[].
 * @author david.nix@hci.utah.edu*/
public class PositionScoreTextData extends USeqData{

	//fields
	private PositionScoreText[] sortedPositionScoreTexts;
	private int[] basePositions;
	private float[] scores;

	//constructors
	public PositionScoreTextData(){}

	/**Note, be sure to sort the PositionScoreText[].*/
	public PositionScoreTextData(PositionScoreText[] sortedPositionScoreTexts, SliceInfo sliceInfo){
		this.sortedPositionScoreTexts = sortedPositionScoreTexts;
		this.sliceInfo = sliceInfo;
	}
	public PositionScoreTextData(File binaryFile) throws IOException{
		sliceInfo = new SliceInfo(binaryFile.getName());
		read (binaryFile);
	}
	public PositionScoreTextData(DataInputStream dis, SliceInfo sliceInfo){
		this.sliceInfo = sliceInfo;
		read (dis);
	}

	//methods
	/**Updates the SliceInfo setting just the FirstStartPosition, LastStartPosition, and NumberRecords.*/
	public static void updateSliceInfo (PositionScoreText[] sortedPositionScoreTexts, SliceInfo sliceInfo){
		sliceInfo.setFirstStartPosition(sortedPositionScoreTexts[0].position);
		sliceInfo.setLastStartPosition(sortedPositionScoreTexts[sortedPositionScoreTexts.length-1].position);
		sliceInfo.setNumberRecords(sortedPositionScoreTexts.length);
	}
	/**Returns the position of the last position in the sortedPositionScoreTexts array.*/
	public int fetchLastBase(){
		return sortedPositionScoreTexts[sortedPositionScoreTexts.length-1].position;
	}
	/**Writes 6 or 12 column xxx.bed formatted lines to the PrintWriter.*/
	public void writeBed (PrintWriter out, boolean fixScore){
		String chrom = sliceInfo.getChromosome();
		String strand = sliceInfo.getStrand();
		for (int i=0; i< sortedPositionScoreTexts.length; i++){
			//chrom start stop name score strand
			//bed12?
			String[] tokens = Text2USeq.PATTERN_TAB.split(sortedPositionScoreTexts[i].text);
			if (fixScore){
				int score = USeqUtilities.fixBedScore(sortedPositionScoreTexts[i].score);
				if (tokens.length == 7) out.println(chrom+"\t"+sortedPositionScoreTexts[i].position+"\t"+(sortedPositionScoreTexts[i].position + 1)+"\t"+ tokens[0] +"\t"+score+"\t"+strand+"\t"+tokens[1]+"\t"+tokens[2]+"\t"+tokens[3]+"\t"+tokens[4]+"\t"+tokens[5]+"\t"+tokens[6]);
				else out.println(chrom+"\t"+sortedPositionScoreTexts[i].position+"\t"+(sortedPositionScoreTexts[i].position + 1)+"\t"+ sortedPositionScoreTexts[i].text +"\t"+score+"\t"+strand);
			}
			else {
				if (tokens.length == 7) out.println(chrom+"\t"+sortedPositionScoreTexts[i].position+"\t"+(sortedPositionScoreTexts[i].position + 1)+"\t"+ tokens[0] +"\t"+sortedPositionScoreTexts[i].score+"\t"+strand+"\t"+tokens[1]+"\t"+tokens[2]+"\t"+tokens[3]+"\t"+tokens[4]+"\t"+tokens[5]+"\t"+tokens[6]);
				else out.println(chrom+"\t"+sortedPositionScoreTexts[i].position+"\t"+(sortedPositionScoreTexts[i].position + 1)+"\t"+ sortedPositionScoreTexts[i].text +"\t"+sortedPositionScoreTexts[i].score+"\t"+strand);
			}
		}
	}

	/**Writes native format to the PrintWriter*/
	public void writeNative (PrintWriter out){
		String chrom = sliceInfo.getChromosome();
		String strand = sliceInfo.getStrand();
		if (strand.equals(".")){
			out.println("#Chr\tPosition\tScore\tText(s)");
			for (int i=0; i< sortedPositionScoreTexts.length; i++) out.println(chrom+"\t"+sortedPositionScoreTexts[i].position+"\t"+sortedPositionScoreTexts[i].score+"\t"+sortedPositionScoreTexts[i].text);
		}
		else {
			out.println("#Chr\tPosition\tScore\tText(s)\tStrand");
			for (int i=0; i< sortedPositionScoreTexts.length; i++){
				//chrom start stop name score strand
				out.println(chrom+"\t"+sortedPositionScoreTexts[i].position+"\t"+sortedPositionScoreTexts[i].score+"\t"+sortedPositionScoreTexts[i].text+"\t"+strand);
			}
		}
	}

	/**Writes position score format to the PrintWriter, 1bp coor*/
	public void writePositionScore (PrintWriter out){
		int prior = -1; 
		for (int i=0; i< sortedPositionScoreTexts.length; i++){
			if (prior != sortedPositionScoreTexts[i].position) {
				out.println((sortedPositionScoreTexts[i].position +1) +"\t"+sortedPositionScoreTexts[i].score);
				prior = sortedPositionScoreTexts[i].position;
			}
		}
	}

	//methods

	/**Writes the PositionScoreText[] to a binary file.
	 * @param saveDirectory, the binary file will be written using the chromStrandStartBP-StopBP.extension notation to this directory
	 * @param attemptToSaveAsShort, scans to see if the offsets exceed 65536 bp, a bit slower to write but potentially a considerable size reduction, set to false for max speed
	 * @return the binaryFile written to the saveDirectory, null if something bad happened
	 * */
	public File write (File saveDirectory, boolean attemptToSaveAsShort) {
		//check to see if this can be saved using shorts instead of ints?
		boolean useShort = false;
		if (attemptToSaveAsShort){			
			int bp = sortedPositionScoreTexts[0].position;
			useShort = true;
			for (int i=1; i< sortedPositionScoreTexts.length; i++){
				int currentStart = sortedPositionScoreTexts[i].position;
				int diff = currentStart - bp;
				if (diff > 65536) {
					useShort = false;
					break;
				}
				bp = currentStart;
			}
		}

		//make and put file type/extension in header
		String fileType;
		if (useShort) fileType = USeqUtilities.SHORT + USeqUtilities.FLOAT;
		else fileType = USeqUtilities.INT + USeqUtilities.FLOAT;
		fileType = fileType + USeqUtilities.TEXT;
		sliceInfo.setBinaryType(fileType);
		binaryFile = new File(saveDirectory, sliceInfo.getSliceName());

		FileOutputStream workingFOS = null;
		DataOutputStream workingDOS = null;
		try {
			//make IO
			workingFOS = new FileOutputStream(binaryFile);
			workingDOS = new DataOutputStream( new BufferedOutputStream (workingFOS));

			//write String header, currently this isn't used
			workingDOS.writeUTF(header);

			//write first position score, always an int
			workingDOS.writeInt(sortedPositionScoreTexts[0].position);
			workingDOS.writeFloat(sortedPositionScoreTexts[0].score);
			workingDOS.writeUTF(sortedPositionScoreTexts[0].text);

			//write shorts?
			if (useShort) {			
				int bp = sortedPositionScoreTexts[0].position;
				for (int i=1; i< sortedPositionScoreTexts.length; i++){
					int currentStart = sortedPositionScoreTexts[i].position;
					//subtract 32768 to extend range of short (-32768 to 32768)
					int diff = currentStart - bp - 32768;
					workingDOS.writeShort((short)(diff));
					workingDOS.writeFloat(sortedPositionScoreTexts[i].score);
					workingDOS.writeUTF(sortedPositionScoreTexts[i].text);
					bp = currentStart;
				}
			}

			//no, write ints
			else {
				int bp = sortedPositionScoreTexts[0].position;
				for (int i=1; i< sortedPositionScoreTexts.length; i++){
					int currentStart = sortedPositionScoreTexts[i].position;
					int diff = currentStart - bp;
					workingDOS.writeInt(diff);
					workingDOS.writeFloat(sortedPositionScoreTexts[i].score);
					workingDOS.writeUTF(sortedPositionScoreTexts[i].text);
					bp = currentStart;
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			binaryFile = null;
		} finally {
			USeqUtilities.safeClose(workingDOS);
			USeqUtilities.safeClose(workingFOS);
		}
		return binaryFile;
	}

	/**Writes the PositionScoreText[] to a ZipOutputStream.
	 * @param	attemptToSaveAsShort	if true, scans to see if the offsets exceed 65536 bp, a bit slower to write but potentially a considerable size reduction, set to false for max speed
	 * */
	public void write (ZipOutputStream out, DataOutputStream dos, boolean attemptToSaveAsShort) {
		//check to see if this can be saved using shorts instead of ints?
		boolean useShort = false;
		if (attemptToSaveAsShort){			
			int bp = sortedPositionScoreTexts[0].position;
			useShort = true;
			for (int i=1; i< sortedPositionScoreTexts.length; i++){
				int currentStart = sortedPositionScoreTexts[i].position;
				int diff = currentStart - bp;
				if (diff > 65536) {
					useShort = false;
					break;
				}
				bp = currentStart;
			}
		}

		//make and put file type/extension in header
		String fileType;
		if (useShort) fileType = USeqUtilities.SHORT + USeqUtilities.FLOAT;
		else fileType = USeqUtilities.INT + USeqUtilities.FLOAT;
		fileType = fileType + USeqUtilities.TEXT;
		sliceInfo.setBinaryType(fileType);
		binaryFile = null;

		try {
			//make new ZipEntry
			out.putNextEntry(new ZipEntry(sliceInfo.getSliceName()));

			//write String header, currently this isn't used
			dos.writeUTF(header);

			//write first position score, always an int
			dos.writeInt(sortedPositionScoreTexts[0].position);
			dos.writeFloat(sortedPositionScoreTexts[0].score);
			dos.writeUTF(sortedPositionScoreTexts[0].text);

			//write shorts?
			if (useShort) {			
				int bp = sortedPositionScoreTexts[0].position;
				for (int i=1; i< sortedPositionScoreTexts.length; i++){
					int currentStart = sortedPositionScoreTexts[i].position;
					//subtract 32768 to extend range of short (-32768 to 32768)
					int diff = currentStart - bp - 32768;
					dos.writeShort((short)(diff));
					dos.writeFloat(sortedPositionScoreTexts[i].score);
					dos.writeUTF(sortedPositionScoreTexts[i].text);
					bp = currentStart;
				}
			}

			//no, write ints
			else {
				int bp = sortedPositionScoreTexts[0].position;
				for (int i=1; i< sortedPositionScoreTexts.length; i++){
					int currentStart = sortedPositionScoreTexts[i].position;
					int diff = currentStart - bp;
					dos.writeInt(diff);
					dos.writeFloat(sortedPositionScoreTexts[i].score);
					dos.writeUTF(sortedPositionScoreTexts[i].text);
					bp = currentStart;
				}
			}

			//close ZipEntry but not streams!
			out.closeEntry();
		} catch (IOException e) {
			e.printStackTrace();
			USeqUtilities.safeClose(out);
			USeqUtilities.safeClose(dos);
		} 
	}

	/**Assumes all are of the same chromosome and strand! Sorts PositionScoreTextData prior to merging*/
	public static PositionScoreTextData merge (ArrayList<PositionScoreTextData> pdAL){
		//convert to arrays and sort
		PositionScoreTextData[] pdArray = new PositionScoreTextData[pdAL.size()];
		pdAL.toArray(pdArray);
		Arrays.sort(pdArray);
		//fetch total size of PositionScore[]
		int num = 0;
		for (int i=0; i< pdArray.length; i++) num += pdArray[i].sortedPositionScoreTexts.length;
		//concatinate
		PositionScoreText[] concatinate = new PositionScoreText[num];
		int index = 0;
		for (int i=0; i< pdArray.length; i++){
			PositionScoreText[] slice = pdArray[i].sortedPositionScoreTexts;
			System.arraycopy(slice, 0, concatinate, index, slice.length);
			index += slice.length;
		}
		//get and modify header
		SliceInfo sliceInfo = pdArray[0].sliceInfo;
		PositionScoreTextData.updateSliceInfo(concatinate, sliceInfo);
		//return new PositionData
		return new PositionScoreTextData(concatinate, sliceInfo);
	}

	public static PositionScoreTextData mergeUSeqData(ArrayList<USeqData> useqDataAL) {
		int num = useqDataAL.size();
		//convert ArrayList
		ArrayList<PositionScoreTextData> a = new ArrayList<PositionScoreTextData>(num);
		for (int i=0; i< num; i++) a.add((PositionScoreTextData) useqDataAL.get(i));
		return merge (a);
	}

	/**Reads a DataInputStream into this PositionScoreData.*/
	public void read (DataInputStream dis) {
		try {
			//read text header, currently not used
			header = dis.readUTF();	

			//make array
			int numberPositions = sliceInfo.getNumberRecords();
			sortedPositionScoreTexts = new PositionScoreText[numberPositions];

			//make first position
			sortedPositionScoreTexts[0] = new PositionScoreText(dis.readInt(), dis.readFloat(), dis.readUTF());

			//what kind of data to follow? 
			String fileType = sliceInfo.getBinaryType();

			//ints?
			if (USeqUtilities.POSITION_SCORE_TEXT_INT_FLOAT_TEXT.matcher(fileType).matches()){	
				//read and resolve offsets to real bps
				for (int i=1; i< numberPositions; i++){
					sortedPositionScoreTexts[i] = new PositionScoreText(sortedPositionScoreTexts[i-1].position + dis.readInt(), dis.readFloat(), dis.readUTF());	
				}
			}
			//shorts?
			else if (USeqUtilities.POSITION_SCORE_TEXT_SHORT_FLOAT_TEXT.matcher(fileType).matches()){
				//read and resolve offsets to real bps
				for (int i=1; i< numberPositions; i++){
					sortedPositionScoreTexts[i] = new PositionScoreText(sortedPositionScoreTexts[i-1].position + dis.readShort() + 32768,  dis.readFloat(), dis.readUTF());	
				}
			}
			//unknown!
			else {
				throw new IOException ("Incorrect file type for creating a PositionScoreText[] -> '"+fileType+"' in "+binaryFile +"\n");
			}
		} catch (IOException e) {
			e.printStackTrace();
			USeqUtilities.safeClose(dis);
		}
	}

	public PositionScoreText[] getPositionScoreTexts() {
		return sortedPositionScoreTexts;
	}

	public void setPositionScoreTexts(PositionScoreText[] sortedPositionScoreTexts) {
		this.sortedPositionScoreTexts = sortedPositionScoreTexts;
		updateSliceInfo(sortedPositionScoreTexts, sliceInfo);
	}
	/**Returns whether data remains.*/
	public boolean trim(int beginningBP, int endingBP) {
		ArrayList<PositionScoreText> al = new ArrayList<PositionScoreText>();
		for (int i=0; i< sortedPositionScoreTexts.length; i++){
			if (sortedPositionScoreTexts[i].isContainedBy(beginningBP, endingBP)) al.add(sortedPositionScoreTexts[i]);
		}
		if (al.size() == 0) return false;
		sortedPositionScoreTexts = new PositionScoreText[al.size()];
		al.toArray(sortedPositionScoreTexts);
		updateSliceInfo(sortedPositionScoreTexts, sliceInfo);
		return true;
	}

	public int[] getBasePositions(){
		if (basePositions == null){
			basePositions = new int[sortedPositionScoreTexts.length];
			scores = new float[sortedPositionScoreTexts.length];
			for (int i=0; i<basePositions.length; i++) {
				basePositions[i] = sortedPositionScoreTexts[i].position;
				scores[i] = sortedPositionScoreTexts[i].score;
			}
		}
		return basePositions;
	}
	public float[] getBaseScores(){
		if (scores== null) getBasePositions();
		return scores;
	}
}
