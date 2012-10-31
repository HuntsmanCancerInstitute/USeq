package edu.utah.seq.useq.data;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;
import edu.utah.seq.useq.*;
import edu.utah.seq.useq.apps.*;

/**Container for a sorted PositionText[].
 * @author david.nix@hci.utah.edu*/
public class PositionTextData extends USeqData{

	//fields
	private PositionText[] sortedPositionTexts;

	//constructors
	public PositionTextData(){}

	/**Note, be sure to sort the PositionText[].*/
	public PositionTextData(PositionText[] sortedPositionTexts, SliceInfo sliceInfo){
		this.sortedPositionTexts = sortedPositionTexts;
		this.sliceInfo = sliceInfo;
	}
	public PositionTextData(File binaryFile) throws Exception{
		sliceInfo = new SliceInfo(binaryFile.getName());
		read (binaryFile);
	}
	public PositionTextData(DataInputStream dis, SliceInfo sliceInfo) {
		this.sliceInfo = sliceInfo;
		read (dis);
	}

	//methods
	/**Updates the SliceInfo setting just the FirstStartPosition, LastStartPosition, and NumberRecords.*/
	public static void updateSliceInfo (PositionText[] sortedPositionTexts, SliceInfo sliceInfo){
		sliceInfo.setFirstStartPosition(sortedPositionTexts[0].position);
		sliceInfo.setLastStartPosition(sortedPositionTexts[sortedPositionTexts.length-1].position);
		sliceInfo.setNumberRecords(sortedPositionTexts.length);
	}
	
	/**Returns the position of the last position in the sortedPositionTexts array.*/
	public int fetchLastBase(){
		return sortedPositionTexts[sortedPositionTexts.length-1].position;
	}
	/**Writes 6 or 12 column xxx.bed formatted lines to the PrintWriter*/
	public void writeBed (PrintWriter out){
		String chrom = sliceInfo.getChromosome();
		String strand = sliceInfo.getStrand();
		for (int i=0; i< sortedPositionTexts.length; i++){
			//chrom start stop name score strand
			//bed12?
			String[] tokens = Text2USeq.PATTERN_TAB.split(sortedPositionTexts[i].text);
			if (tokens.length == 7) out.println(chrom+"\t"+sortedPositionTexts[i].position+"\t"+(sortedPositionTexts[i].position + 1)+"\t"+tokens[0]+"\t0\t"+strand+"\t"+tokens[1]+"\t"+tokens[2]+"\t"+tokens[3]+"\t"+tokens[4]+"\t"+tokens[5]+"\t"+tokens[6]);
			else out.println(chrom+"\t"+sortedPositionTexts[i].position+"\t"+(sortedPositionTexts[i].position + 1)+"\t"+sortedPositionTexts[i].text+"\t0\t"+strand);
		}
	}
	
	/**Writes native format to the PrintWriter*/
	public void writeNative (PrintWriter out){
		String chrom = sliceInfo.getChromosome();
		String strand = sliceInfo.getStrand();
		if (strand.equals(".")){
			out.println("#Chr\tPosition\tText(s)");
			for (int i=0; i< sortedPositionTexts.length; i++) out.println(chrom+"\t"+sortedPositionTexts[i].position+"\t"+sortedPositionTexts[i].text);
		}
		else {
			out.println("#Chr\tPosition\tText(s)\tStrand");
			for (int i=0; i< sortedPositionTexts.length; i++){
				//chrom start stop name score strand
				out.println(chrom+"\t"+sortedPositionTexts[i].position+"\t"+sortedPositionTexts[i].text+"\t"+strand);
			}
		}
	}
	
	/**Writes position score format to the PrintWriter, 1bp coor*/
	public void writePositionScore (PrintWriter out){
		int prior = -1;
		for (int i=0; i< sortedPositionTexts.length; i++){
				if (prior != sortedPositionTexts[i].position) {
					out.println((sortedPositionTexts[i].position +1) +"\t0");
					prior = sortedPositionTexts[i].position;
				}
		}
	}
	
	
	/**Writes the PositionText[] to a binary file.
	 * @param saveDirectory, the binary file will be written using the chromStrandStartBP-StopBP.extension notation to this directory
	 * @param attemptToSaveAsShort, scans to see if the offsets exceed 65536 bp, a bit slower to write but potentially a considerable size reduction, set to false for max speed
	 * @return the binaryFile written to the saveDirectory
	 * */
	public File write (File saveDirectory, boolean attemptToSaveAsShort) {
		//check to see if this can be saved using shorts instead of ints?
		boolean useShort = false;
		if (attemptToSaveAsShort){			
			int bp = sortedPositionTexts[0].position;
			useShort = true;
			for (int i=1; i< sortedPositionTexts.length; i++){
				int currentStart = sortedPositionTexts[i].position;
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
		if (useShort) fileType = USeqUtilities.SHORT + USeqUtilities.TEXT;
		else fileType = USeqUtilities.INT + USeqUtilities.TEXT;
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
			workingDOS.writeInt(sortedPositionTexts[0].position);
			workingDOS.writeUTF(sortedPositionTexts[0].text);

			//write shorts?
			if (useShort) {			
				int bp = sortedPositionTexts[0].position;
				for (int i=1; i< sortedPositionTexts.length; i++){
					int currentStart = sortedPositionTexts[i].position;
					//subtract 32768 to extend range of short (-32768 to 32768)
					int diff = currentStart - bp - 32768;
					workingDOS.writeShort((short)(diff));
					workingDOS.writeUTF(sortedPositionTexts[i].text);
					bp = currentStart;
				}
			}

			//no, write ints
			else {
				int bp = sortedPositionTexts[0].position;
				for (int i=1; i< sortedPositionTexts.length; i++){
					int currentStart = sortedPositionTexts[i].position;
					int diff = currentStart - bp;
					workingDOS.writeInt(diff);
					workingDOS.writeUTF(sortedPositionTexts[i].text);
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
	
	/**Assumes all are of the same chromosome and strand! Sorts PositionTextData prior to merging*/
	public static PositionTextData merge (ArrayList<PositionTextData> pdAL){
		//convert to arrays and sort
		PositionTextData[] pdArray = new PositionTextData[pdAL.size()];
		pdAL.toArray(pdArray);
		Arrays.sort(pdArray);
		//fetch total size of PositionText[]
		int num = 0;
		for (int i=0; i< pdArray.length; i++) num += pdArray[i].sortedPositionTexts.length;
		//concatinate
		PositionText[] concatinate = new PositionText[num];
		int index = 0;
		for (int i=0; i< pdArray.length; i++){
			PositionText[] slice = pdArray[i].sortedPositionTexts;
			System.arraycopy(slice, 0, concatinate, index, slice.length);
			index += slice.length;
		}
		//get and modify header
		SliceInfo sliceInfo = pdArray[0].sliceInfo;
		PositionTextData.updateSliceInfo(concatinate, sliceInfo);
		//return new PositionTextData
		return new PositionTextData(concatinate, sliceInfo);
	}
	
	public static PositionTextData mergeUSeqData(ArrayList<USeqData> useqDataAL) {
		int num = useqDataAL.size();
		//convert ArrayList
		ArrayList<PositionTextData> a = new ArrayList<PositionTextData>(num);
		for (int i=0; i< num; i++) a.add((PositionTextData) useqDataAL.get(i));
		return merge (a);
	}

	/**Writes the PositionText[] to a ZipOutputStream.
	 * @param	attemptToSaveAsShort	if true, scans to see if the offsets exceed 65536 bp, a bit slower to write but potentially a considerable size reduction, set to false for max speed
	 * */
	public void write (ZipOutputStream out, DataOutputStream dos, boolean attemptToSaveAsShort) {
		//check to see if this can be saved using shorts instead of ints?
		boolean useShort = false;
		if (attemptToSaveAsShort){			
			int bp = sortedPositionTexts[0].position;
			useShort = true;
			for (int i=1; i< sortedPositionTexts.length; i++){
				int currentStart = sortedPositionTexts[i].position;
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
		if (useShort) fileType = USeqUtilities.SHORT + USeqUtilities.TEXT;
		else fileType = USeqUtilities.INT + USeqUtilities.TEXT;
		sliceInfo.setBinaryType(fileType);
		binaryFile = null;

		try {
			//make new ZipEntry
			out.putNextEntry(new ZipEntry(sliceInfo.getSliceName()));

			//write String header, currently this isn't used
			dos.writeUTF(header);

			//write first position score, always an int
			dos.writeInt(sortedPositionTexts[0].position);
			dos.writeUTF(sortedPositionTexts[0].text);

			//write shorts?
			if (useShort) {			
				int bp = sortedPositionTexts[0].position;
				for (int i=1; i< sortedPositionTexts.length; i++){
					int currentStart = sortedPositionTexts[i].position;
					//subtract 32768 to extend range of short (-32768 to 32768)
					int diff = currentStart - bp - 32768;
					dos.writeShort((short)(diff));
					dos.writeUTF(sortedPositionTexts[i].text);
					bp = currentStart;
				}
			}

			//no, write ints
			else {
				int bp = sortedPositionTexts[0].position;
				for (int i=1; i< sortedPositionTexts.length; i++){
					int currentStart = sortedPositionTexts[i].position;
					int diff = currentStart - bp;
					dos.writeInt(diff);
					dos.writeUTF(sortedPositionTexts[i].text);
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

	/**Reads a DataInputStream into this PositionTextData.*/
	public void read (DataInputStream dis) {
		try {
			//read text header, currently not used
			header = dis.readUTF();

			//make array
			int numberPositions = sliceInfo.getNumberRecords();
			sortedPositionTexts = new PositionText[numberPositions];

			//make first position
			sortedPositionTexts[0] = new PositionText(dis.readInt(), dis.readUTF());

			//what kind of data to follow? 
			String fileType = sliceInfo.getBinaryType();

			//ints?
			if (USeqUtilities.POSITION_TEXT_INT_TEXT.matcher(fileType).matches()){	
				//read and resolve offsets to real bps
				for (int i=1; i< numberPositions; i++){
					sortedPositionTexts[i] = new PositionText(sortedPositionTexts[i-1].position + dis.readInt(), dis.readUTF());	
				}
			}
			//shorts?
			else if (USeqUtilities.POSITION_TEXT_SHORT_TEXT.matcher(fileType).matches()){
				//read and resolve offsets to real bps
				for (int i=1; i< numberPositions; i++){
					sortedPositionTexts[i] = new PositionText(sortedPositionTexts[i-1].position + dis.readShort() + 32768,  dis.readUTF());	
				}
			}
			//unknown!
			else {
				throw new IOException ("Incorrect file type for creating a PositionText[] -> '"+fileType+"' in "+binaryFile +"\n");
			}
		} catch (IOException e) {
			e.printStackTrace();
			USeqUtilities.safeClose(dis);
		}
	}

	public PositionText[] getPositionTexts() {
		return sortedPositionTexts;
	}

	public void setPositionTexts(PositionText[] sortedPositionTexts) {
		this.sortedPositionTexts = sortedPositionTexts;
		updateSliceInfo(sortedPositionTexts, sliceInfo);
	}

	/**Returns whether data remains.*/
	public boolean trim(int beginningBP, int endingBP) {
		ArrayList<PositionText> al = new ArrayList<PositionText>();
		for (int i=0; i< sortedPositionTexts.length; i++){
			if (sortedPositionTexts[i].isContainedBy(beginningBP, endingBP)) al.add(sortedPositionTexts[i]);
		}
		if (al.size() == 0) return false;
		sortedPositionTexts = new PositionText[al.size()];
		al.toArray(sortedPositionTexts);
		updateSliceInfo(sortedPositionTexts, sliceInfo);
		return true;
	}
}
