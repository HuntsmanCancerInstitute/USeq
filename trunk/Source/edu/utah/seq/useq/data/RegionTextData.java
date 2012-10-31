package edu.utah.seq.useq.data;

import java.io.*;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;
import edu.utah.seq.useq.*;
import edu.utah.seq.useq.apps.*;

/**Container for a sorted RegionText[].
 * @author david.nix@hci.utah.edu*/
public class RegionTextData extends USeqData{

	//fields
	private RegionText[] sortedRegionTexts;

	//constructors
	public RegionTextData(){}

	/**Note, be sure to sort the RegionText[].*/
	public RegionTextData(RegionText[] sortedRegionTexts, SliceInfo sliceInfo){
		this.sortedRegionTexts = sortedRegionTexts;
		this.sliceInfo = sliceInfo;
	}
	public RegionTextData(File binaryFile) throws IOException{
		sliceInfo = new SliceInfo(binaryFile.getName());
		read (binaryFile);
	}
	public RegionTextData(DataInputStream dis, SliceInfo sliceInfo) {
		this.sliceInfo = sliceInfo;
		read (dis);
	}

	//methods
	/**Updates the SliceInfo setting just the FirstStartPosition, LastStartPosition, and NumberRecords.*/
	public static void updateSliceInfo (RegionText[] sortedRegionTexts, SliceInfo sliceInfo){
		sliceInfo.setFirstStartPosition(sortedRegionTexts[0].getStart());
		sliceInfo.setLastStartPosition(sortedRegionTexts[sortedRegionTexts.length-1].start);
		sliceInfo.setNumberRecords(sortedRegionTexts.length);
	}
	/**Returns the bp of the last end position in the array.*/
	public int fetchLastBase(){
		int lastBase = -1;
		for (RegionText r : sortedRegionTexts){
			int end = r.getStop();
			if (end > lastBase) lastBase = end;
		}
		return lastBase;
	}
	/**Writes 6 or 12 column xxx.bed formatted lines to the PrintWriter*/
	public void writeBed (PrintWriter out){
		String chrom = sliceInfo.getChromosome();
		String strand = sliceInfo.getStrand();
		for (int i=0; i< sortedRegionTexts.length; i++){
			//bed12?
			String[] tokens = Text2USeq.PATTERN_TAB.split(sortedRegionTexts[i].text);
			if (tokens.length == 7) out.println(chrom+"\t"+sortedRegionTexts[i].start+"\t"+sortedRegionTexts[i].stop+"\t"+ tokens[0] +"\t0\t"+strand+"\t"+tokens[1]+"\t"+tokens[2]+"\t"+tokens[3]+"\t"+tokens[4]+"\t"+tokens[5]+"\t"+tokens[6]);
			else out.println(chrom+"\t"+sortedRegionTexts[i].start+"\t"+sortedRegionTexts[i].stop+"\t"+ sortedRegionTexts[i].text +"\t0\t"+strand);
		}
	}
	
	/**Writes native format to the PrintWriter*/
	public void writeNative (PrintWriter out){
		String chrom = sliceInfo.getChromosome();
		String strand = sliceInfo.getStrand();
		if (strand.equals(".")){
			out.println("#Chr\tStart\tStop\tText(s)");
			for (int i=0; i< sortedRegionTexts.length; i++) out.println(chrom+"\t"+sortedRegionTexts[i].start+"\t"+sortedRegionTexts[i].stop+"\t"+sortedRegionTexts[i].text);
		}
		else {
			out.println("#Chr\tStart\tStop\tText(s)\tStrand");
			for (int i=0; i< sortedRegionTexts.length; i++) out.println(chrom+"\t"+sortedRegionTexts[i].start+"\t"+sortedRegionTexts[i].stop+"\t"+sortedRegionTexts[i].text+"\t"+strand);
		}
	}

	
	/**Writes the RegionText[] to a binary file.  Each region's start/stop is converted to a running offset/length which are written as either ints or shorts.
	 * @param saveDirectory, the binary file will be written using the chromStrandStartBP-StopBP.extension notation to this directory
	 * @param attemptToSaveAsShort, scans to see if the offsets and region lengths exceed 65536 bp, a bit slower to write but potentially a considerable size reduction, set to false for max speed
	 * @return the binaryFile written to the saveDirectory
	 * */
	public File write (File saveDirectory, boolean attemptToSaveAsShort) {
		//check to see if this can be saved using shorts instead of ints?
		boolean useShortBeginning = false;
		boolean useShortLength = false;
		if (attemptToSaveAsShort){			
			int bp = sortedRegionTexts[0].start;
			useShortBeginning = true;
			for (int i=1; i< sortedRegionTexts.length; i++){
				int currentStart = sortedRegionTexts[i].start;
				int diff = currentStart - bp;
				if (diff > 65536) {
					useShortBeginning = false;
					break;
				}
				bp = currentStart;
			}
			//check to short length
			useShortLength = true;
			for (int i=0; i< sortedRegionTexts.length; i++){
				int diff = sortedRegionTexts[i].stop - sortedRegionTexts[i].start;
				if (diff > 65536) {
					useShortLength = false;
					break;
				}
			}
		}
		//make and put file type/extension in header
		String fileType;
		if (useShortBeginning) fileType = USeqUtilities.SHORT;
		else fileType = USeqUtilities.INT;
		if (useShortLength) fileType = fileType+ USeqUtilities.SHORT;
		else fileType = fileType+ USeqUtilities.INT;
		fileType = fileType+ USeqUtilities.TEXT;
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

			//write first position, always an int
			workingDOS.writeInt(sortedRegionTexts[0].start);

			//write short position?
			int bp = sortedRegionTexts[0].start;
			if (useShortBeginning) {			
				//also short length?
				//no
				if (useShortLength == false){
					//write first record's length
					workingDOS.writeInt(sortedRegionTexts[0].stop- sortedRegionTexts[0].start);
					workingDOS.writeUTF(sortedRegionTexts[0].text);
					for (int i=1; i< sortedRegionTexts.length; i++){
						int currentStart = sortedRegionTexts[i].start;
						//subtract 32768 to extend range of short (-32768 to 32768)
						int diff = currentStart - bp - 32768;
						workingDOS.writeShort((short)(diff));
						workingDOS.writeInt(sortedRegionTexts[i].stop- sortedRegionTexts[i].start);
						workingDOS.writeUTF(sortedRegionTexts[i].text);
						bp = currentStart;
					}
				}
				//yes short length
				else {
					//write first record's length, subtracting 32768 to extent the range of the signed short
					workingDOS.writeShort((short)(sortedRegionTexts[0].stop- sortedRegionTexts[0].start - 32768));
					workingDOS.writeUTF(sortedRegionTexts[0].text);
					for (int i=1; i< sortedRegionTexts.length; i++){
						int currentStart = sortedRegionTexts[i].start;
						//subtract 32768 to extend range of short (-32768 to 32768)
						int diff = currentStart - bp - 32768;
						workingDOS.writeShort((short)(diff));
						workingDOS.writeShort((short)(sortedRegionTexts[i].stop- sortedRegionTexts[i].start - 32768));
						workingDOS.writeUTF(sortedRegionTexts[i].text);
						bp = currentStart;
					}
				}
			}

			//no, write int for position
			else {
				//short length? no
				if (useShortLength == false){
					//write first record's length
					workingDOS.writeInt(sortedRegionTexts[0].stop- sortedRegionTexts[0].start);
					workingDOS.writeUTF(sortedRegionTexts[0].text);
					for (int i=1; i< sortedRegionTexts.length; i++){
						int currentStart = sortedRegionTexts[i].start;
						int diff = currentStart - bp;
						workingDOS.writeInt(diff);
						workingDOS.writeInt(sortedRegionTexts[i].stop- sortedRegionTexts[i].start);
						workingDOS.writeUTF(sortedRegionTexts[i].text);
						bp = currentStart;
					}
				}
				//yes
				else {
					//write first record's length
					workingDOS.writeShort((short)(sortedRegionTexts[0].stop- sortedRegionTexts[0].start - 32768));
					workingDOS.writeUTF(sortedRegionTexts[0].text);
					for (int i=1; i< sortedRegionTexts.length; i++){
						int currentStart = sortedRegionTexts[i].start;
						int diff = currentStart - bp;
						workingDOS.writeInt(diff);
						workingDOS.writeShort((short)(sortedRegionTexts[i].stop- sortedRegionTexts[i].start - 32768));
						workingDOS.writeUTF(sortedRegionTexts[i].text);
						bp = currentStart;
					}
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
	
	/**Assumes all are of the same chromosome and strand!*/
	public static RegionTextData merge (ArrayList<RegionTextData> pdAL){
		//convert to arrays and sort
		RegionTextData[] pdArray = new RegionTextData[pdAL.size()];
		pdAL.toArray(pdArray);
		Arrays.sort(pdArray);
		//fetch total size of RegionText[]
		int num = 0;
		for (int i=0; i< pdArray.length; i++) num += pdArray[i].sortedRegionTexts.length;
		//concatinate
		RegionText[] concatinate = new RegionText[num];
		int index = 0;
		for (int i=0; i< pdArray.length; i++){
			RegionText[] slice = pdArray[i].sortedRegionTexts;
			System.arraycopy(slice, 0, concatinate, index, slice.length);
			index += slice.length;
		}
		//get and modify header
		SliceInfo sliceInfo = pdArray[0].sliceInfo;
		RegionTextData.updateSliceInfo(concatinate, sliceInfo);
		//return new RegionTextData
		return new RegionTextData(concatinate, sliceInfo);
	}
	
	public static RegionTextData mergeUSeqData(ArrayList<USeqData> useqDataAL) {
		int num = useqDataAL.size();
		//convert ArrayList
		ArrayList<RegionTextData> a = new ArrayList<RegionTextData>(num);
		for (int i=0; i< num; i++) a.add((RegionTextData) useqDataAL.get(i));
		return merge (a);
	}

	/**Writes the Region[] to a ZipOutputStream.
	 * @param	attemptToSaveAsShort	if true, scans to see if the offsets exceed 65536 bp, a bit slower to write but potentially a considerable size reduction, set to false for max speed
	 * */
	public void write (ZipOutputStream out, DataOutputStream dos, boolean attemptToSaveAsShort) {
		//check to see if this can be saved using shorts instead of ints?
		boolean useShortBeginning = false;
		boolean useShortLength = false;
		if (attemptToSaveAsShort){			
			int bp = sortedRegionTexts[0].start;
			useShortBeginning = true;
			for (int i=1; i< sortedRegionTexts.length; i++){
				int currentStart = sortedRegionTexts[i].start;
				int diff = currentStart - bp;
				if (diff > 65536) {
					useShortBeginning = false;
					break;
				}
				bp = currentStart;
			}
			//check to short length
			useShortLength = true;
			for (int i=0; i< sortedRegionTexts.length; i++){
				int diff = sortedRegionTexts[i].stop - sortedRegionTexts[i].start;
				if (diff > 65536) {
					useShortLength = false;
					break;
				}
			}
		}
		//make and put file type/extension in header
		String fileType;
		if (useShortBeginning) fileType = USeqUtilities.SHORT;
		else fileType = USeqUtilities.INT;
		if (useShortLength) fileType = fileType+ USeqUtilities.SHORT;
		else fileType = fileType+ USeqUtilities.INT;
		fileType = fileType+ USeqUtilities.TEXT;
		sliceInfo.setBinaryType(fileType);
		binaryFile = null;

		try {
			//make new ZipEntry
			out.putNextEntry(new ZipEntry(sliceInfo.getSliceName()));

			//write String header, currently this isn't used
			dos.writeUTF(header);

			//write first position, always an int
			dos.writeInt(sortedRegionTexts[0].start);

			//write short position?
			int bp = sortedRegionTexts[0].start;
			if (useShortBeginning) {			
				//also short length?
				//no
				if (useShortLength == false){
					//write first record's length
					dos.writeInt(sortedRegionTexts[0].stop- sortedRegionTexts[0].start);
					dos.writeUTF(sortedRegionTexts[0].text);
					for (int i=1; i< sortedRegionTexts.length; i++){
						int currentStart = sortedRegionTexts[i].start;
						//subtract 32768 to extend range of short (-32768 to 32768)
						int diff = currentStart - bp - 32768;
						dos.writeShort((short)(diff));
						dos.writeInt(sortedRegionTexts[i].stop- sortedRegionTexts[i].start);
						dos.writeUTF(sortedRegionTexts[i].text);
						bp = currentStart;
					}
				}
				//yes short length
				else {
					//write first record's length, subtracting 32768 to extent the range of the signed short
					dos.writeShort((short)(sortedRegionTexts[0].stop- sortedRegionTexts[0].start - 32768));
					dos.writeUTF(sortedRegionTexts[0].text);
					for (int i=1; i< sortedRegionTexts.length; i++){
						int currentStart = sortedRegionTexts[i].start;
						//subtract 32768 to extend range of short (-32768 to 32768)
						int diff = currentStart - bp - 32768;
						dos.writeShort((short)(diff));
						dos.writeShort((short)(sortedRegionTexts[i].stop- sortedRegionTexts[i].start - 32768));
						dos.writeUTF(sortedRegionTexts[i].text);
						bp = currentStart;
					}
				}
			}

			//no, write int for position
			else {
				//short length? no
				if (useShortLength == false){
					//write first record's length
					dos.writeInt(sortedRegionTexts[0].stop- sortedRegionTexts[0].start);
					dos.writeUTF(sortedRegionTexts[0].text);
					for (int i=1; i< sortedRegionTexts.length; i++){
						int currentStart = sortedRegionTexts[i].start;
						int diff = currentStart - bp;
						dos.writeInt(diff);
						dos.writeInt(sortedRegionTexts[i].stop- sortedRegionTexts[i].start);
						dos.writeUTF(sortedRegionTexts[i].text);
						bp = currentStart;
					}
				}
				//yes
				else {
					//write first record's length
					dos.writeShort((short)(sortedRegionTexts[0].stop- sortedRegionTexts[0].start - 32768));
					dos.writeUTF(sortedRegionTexts[0].text);
					for (int i=1; i< sortedRegionTexts.length; i++){
						int currentStart = sortedRegionTexts[i].start;
						int diff = currentStart - bp;
						dos.writeInt(diff);
						dos.writeShort((short)(sortedRegionTexts[i].stop- sortedRegionTexts[i].start - 32768));
						dos.writeUTF(sortedRegionTexts[i].text);
						bp = currentStart;
					}
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

	/**Reads a DataInputStream into this RegionScoreData.*/
	public void read (DataInputStream dis) {
		try {
			//read text header, currently not used
			header = dis.readUTF();

			//make array
			int numberRegionTexts = sliceInfo.getNumberRecords();
			sortedRegionTexts = new RegionText[numberRegionTexts];

			//what kind of data to follow? 
			String fileType = sliceInfo.getBinaryType();

			//int Position, int Length
			if (USeqUtilities.REGION_TEXT_INT_INT_TEXT.matcher(fileType).matches()){
				//make first RegionText, position is always an int
				int start = dis.readInt();
				sortedRegionTexts[0] = new RegionText(start, start+dis.readInt(), dis.readUTF());
				//read and resolve offsets to real bps and length to stop
				for (int i=1; i< numberRegionTexts; i++){
					start = sortedRegionTexts[i-1].start + dis.readInt();
					sortedRegionTexts[i] = new RegionText(start, start + dis.readInt(), dis.readUTF());	
				}
			}
			//int Position, short Length
			else if (USeqUtilities.REGION_TEXT_INT_SHORT_TEXT.matcher(fileType).matches()){
				//make first RegionText, position is always an int
				int start = dis.readInt();
				sortedRegionTexts[0] = new RegionText(start, start+ dis.readShort() + 32768, dis.readUTF());
				//read and resolve offsets to real bps and length to stop
				for (int i=1; i< numberRegionTexts; i++){
					start = sortedRegionTexts[i-1].start + dis.readInt();
					sortedRegionTexts[i] = new RegionText(start, start + dis.readShort() + 32768, dis.readUTF());	
				}
			}
			//short Postion, short Length
			else if (USeqUtilities.REGION_TEXT_SHORT_SHORT_TEXT.matcher(fileType).matches()){
				//make first RegionText, position is always an int
				int start = dis.readInt();
				sortedRegionTexts[0] = new RegionText(start, start+ dis.readShort() + 32768, dis.readUTF());
				//read and resolve offsets to real bps and length to stop
				for (int i=1; i< numberRegionTexts; i++){
					start = sortedRegionTexts[i-1].start + dis.readShort() + 32768;
					sortedRegionTexts[i] = new RegionText(start, start + dis.readShort() + 32768, dis.readUTF());
				}
			}
			//short Position, int Length
			else if (USeqUtilities.REGION_TEXT_SHORT_INT_TEXT.matcher(fileType).matches()){
				//make first RegionText, position is always an int
				int start = dis.readInt();
				sortedRegionTexts[0] = new RegionText(start, start+ dis.readInt(), dis.readUTF());
				//read and resolve offsets to real bps and length to stop
				for (int i=1; i< numberRegionTexts; i++){
					start = sortedRegionTexts[i-1].start + dis.readShort() + 32768;
					sortedRegionTexts[i] = new RegionText(start, start + dis.readInt(), dis.readUTF());	
				}
			}
			//unknown!
			else {
				throw new IOException ("Incorrect file type for creating a RegionText[] -> '"+fileType+"' in "+binaryFile +"\n");
			}
		} catch (IOException e) {
			e.printStackTrace();
			USeqUtilities.safeClose(dis);
		}
	}

	public RegionText[] getRegionTexts() {
		return sortedRegionTexts;
	}

	public void setRegionTexts(RegionText[] sortedRegionTexts) {
		this.sortedRegionTexts = sortedRegionTexts;
		updateSliceInfo(sortedRegionTexts, sliceInfo);
	}

	/**Returns whether data remains.*/
	public boolean trim(int beginningBP, int endingBP) {
		ArrayList<RegionText> al = new ArrayList<RegionText>();
		for (int i=0; i< sortedRegionTexts.length; i++){
			if (sortedRegionTexts[i].isContainedBy(beginningBP, endingBP)) al.add(sortedRegionTexts[i]);
		}
		if (al.size() == 0) return false;
		sortedRegionTexts = new RegionText[al.size()];
		al.toArray(sortedRegionTexts);
		updateSliceInfo(sortedRegionTexts, sliceInfo);
		return true;
	}
}
