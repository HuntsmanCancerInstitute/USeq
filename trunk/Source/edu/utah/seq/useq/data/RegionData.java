package edu.utah.seq.useq.data;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;
import edu.utah.seq.useq.*;

/**Container for a sorted Region[] and it's associated SliceInfo.
 * @author david.nix@hci.utah.edu*/
public class RegionData extends USeqData{

	//fields
	private Region[] sortedRegions;

	//constructors
	public RegionData(){}

	/**Note, be sure to sort the Region[].*/
	public RegionData(Region[] sortedRegions, SliceInfo sliceInfo){
		this.sortedRegions = sortedRegions;
		this.sliceInfo = sliceInfo;
	}
	public RegionData(File binaryFile) throws IOException{
		sliceInfo = new SliceInfo(binaryFile.getName());
		read (binaryFile);
	}
	public RegionData(DataInputStream dis, SliceInfo sliceInfo) {
		this.sliceInfo = sliceInfo;
		read (dis);
	}

	//methods
	/**Updates the SliceInfo setting just the FirstStartPosition, LastStartPosition, and NumberRecords.*/
	public static void updateSliceInfo (Region[] sortedRegions, SliceInfo sliceInfo){
		sliceInfo.setFirstStartPosition(sortedRegions[0].getStart());
		sliceInfo.setLastStartPosition(sortedRegions[sortedRegions.length-1].start);
		sliceInfo.setNumberRecords(sortedRegions.length);
	}
	/**Returns the bp of the last end position in the array.*/
	public int fetchLastBase(){
		int lastBase = -1;
		for (Region r : sortedRegions){
			int end = r.getStop();
			if (end > lastBase) lastBase = end;
		}
		return lastBase;
	}

	/**Writes six column xxx.bed formatted lines to the PrintWriter*/
	public void writeBed (PrintWriter out){
		String chrom = sliceInfo.getChromosome();
		String strand = sliceInfo.getStrand();
		for (int i=0; i< sortedRegions.length; i++){
			//chrom start stop name score strand
			out.println(chrom+"\t"+sortedRegions[i].start+"\t"+sortedRegions[i].stop+"\t"+".\t0\t"+strand);
		}
	}
	
	/**Writes native format to the PrintWriter*/
	public void writeNative (PrintWriter out){
		String chrom = sliceInfo.getChromosome();
		String strand = sliceInfo.getStrand();
		if (strand.equals(".")){
			out.println("#Chr\tStart\tStop");
			for (int i=0; i< sortedRegions.length; i++) out.println(chrom+"\t"+sortedRegions[i].start+"\t"+sortedRegions[i].stop);
		}
		else {
			out.println("#Chr\tStart\tStop\tStrand");
			for (int i=0; i< sortedRegions.length; i++) out.println(chrom+"\t"+sortedRegions[i].start+"\t"+sortedRegions[i].stop+"\t"+strand);
		}
	}

	/**Writes the Region[] to a binary file.  Each region's start/stop is converted to a running offset/length which are written as either ints or shorts.
	 * @param saveDirectory, the binary file will be written using the chromStrandStartBP-StopBP.extension notation to this directory
	 * @param attemptToSaveAsShort, scans to see if the offsets and region lengths exceed 65536 bp, a bit slower to write but potentially a considerable size reduction, set to false for max speed
	 * @return the binaryFile written to the saveDirectory
	 * */
	public File write (File saveDirectory, boolean attemptToSaveAsShort) {
		//check to see if this can be saved using shorts instead of ints?
		boolean useShortBeginning = false;
		boolean useShortLength = false;
		if (attemptToSaveAsShort){			
			int bp = sortedRegions[0].start;
			useShortBeginning = true;
			for (int i=1; i< sortedRegions.length; i++){
				int currentStart = sortedRegions[i].start;
				int diff = currentStart - bp;
				if (diff > 65536) {
					useShortBeginning = false;
					break;
				}
				bp = currentStart;
			}
			//check to short length
			useShortLength = true;
			for (int i=0; i< sortedRegions.length; i++){
				int diff = sortedRegions[i].stop - sortedRegions[i].start;
				if (diff > 65536) {
					useShortLength = false;
					break;
				}
			}
		}
		//make and put file type/extension 
		String fileType;
		if (useShortBeginning) fileType = USeqUtilities.SHORT;
		else fileType = USeqUtilities.INT;
		if (useShortLength) fileType = fileType+ USeqUtilities.SHORT;
		else fileType = fileType+ USeqUtilities.INT;
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
			workingDOS.writeInt(sortedRegions[0].start);

			//write short position?
			int bp = sortedRegions[0].start;
			if (useShortBeginning) {			
				//also short length?
				//no
				if (useShortLength == false){
					//write first record's length
					workingDOS.writeInt(sortedRegions[0].stop- sortedRegions[0].start);
					for (int i=1; i< sortedRegions.length; i++){
						int currentStart = sortedRegions[i].start;
						//subtract 32768 to extend range of short (-32768 to 32768)
						int diff = currentStart - bp - 32768;
						workingDOS.writeShort((short)(diff));
						workingDOS.writeInt(sortedRegions[i].stop- sortedRegions[i].start);
						bp = currentStart;
					}
				}
				//yes short length
				else {
					//write first record's length, subtracting 32768 to extent the range of the signed short
					workingDOS.writeShort((short)(sortedRegions[0].stop- sortedRegions[0].start - 32768));
					for (int i=1; i< sortedRegions.length; i++){
						int currentStart = sortedRegions[i].start;
						//subtract 32768 to extend range of short (-32768 to 32768)
						int diff = currentStart - bp - 32768;
						workingDOS.writeShort((short)(diff));
						workingDOS.writeShort((short)(sortedRegions[i].stop- sortedRegions[i].start - 32768));
						bp = currentStart;
					}
				}
			}

			//no, write int for position
			else {
				//short length? no
				if (useShortLength == false){
					//write first record's length
					workingDOS.writeInt(sortedRegions[0].stop- sortedRegions[0].start);
					for (int i=1; i< sortedRegions.length; i++){
						int currentStart = sortedRegions[i].start;
						int diff = currentStart - bp;
						workingDOS.writeInt(diff);
						workingDOS.writeInt(sortedRegions[i].stop- sortedRegions[i].start);
						bp = currentStart;
					}
				}
				//yes
				else {
					//write first record's length
					workingDOS.writeShort((short)(sortedRegions[0].stop- sortedRegions[0].start - 32768));
					for (int i=1; i< sortedRegions.length; i++){
						int currentStart = sortedRegions[i].start;
						int diff = currentStart - bp;
						workingDOS.writeInt(diff);
						workingDOS.writeShort((short)(sortedRegions[i].stop- sortedRegions[i].start - 32768));
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
	public static RegionData merge (ArrayList<RegionData> pdAL){
		//convert to arrays and sort
		RegionData[] pdArray = new RegionData[pdAL.size()];
		pdAL.toArray(pdArray);
		Arrays.sort(pdArray);
		//fetch total size of Region[]
		int num = 0;
		for (int i=0; i< pdArray.length; i++) num += pdArray[i].sortedRegions.length;
		//concatinate
		Region[] concatinate = new Region[num];
		int index = 0;
		for (int i=0; i< pdArray.length; i++){
			Region[] slice = pdArray[i].sortedRegions;
			System.arraycopy(slice, 0, concatinate, index, slice.length);
			index += slice.length;
		}
		//get and modify header
		SliceInfo sliceInfo = pdArray[0].sliceInfo;
		RegionData.updateSliceInfo(concatinate, sliceInfo);
		//return new RegionData
		return new RegionData(concatinate, sliceInfo);
	}
	
	public static RegionData mergeUSeqData(ArrayList<USeqData> useqDataAL) {
		int num = useqDataAL.size();
		//convert ArrayList
		ArrayList<RegionData> a = new ArrayList<RegionData>(num);
		for (int i=0; i< num; i++) a.add((RegionData) useqDataAL.get(i));
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
			int bp = sortedRegions[0].start;
			useShortBeginning = true;
			for (int i=1; i< sortedRegions.length; i++){
				int currentStart = sortedRegions[i].start;
				int diff = currentStart - bp;
				if (diff > 65536) {
					useShortBeginning = false;
					break;
				}
				bp = currentStart;
			}
			//check to short length
			useShortLength = true;
			for (int i=0; i< sortedRegions.length; i++){
				int diff = sortedRegions[i].stop - sortedRegions[i].start;
				if (diff > 65536) {
					useShortLength = false;
					break;
				}
			}
		}
		//make and put file type/extension 
		String fileType;
		if (useShortBeginning) fileType = USeqUtilities.SHORT;
		else fileType = USeqUtilities.INT;
		if (useShortLength) fileType = fileType+ USeqUtilities.SHORT;
		else fileType = fileType+ USeqUtilities.INT;
		sliceInfo.setBinaryType(fileType);
		binaryFile = null;

		try {
			//make new ZipEntry
			out.putNextEntry(new ZipEntry(sliceInfo.getSliceName()));

			//write String header, currently this isn't used
			dos.writeUTF(header);

			//write first position, always an int
			dos.writeInt(sortedRegions[0].start);

			//write short position?
			int bp = sortedRegions[0].start;
			if (useShortBeginning) {			
				//also short length?
				//no
				if (useShortLength == false){
					//write first record's length
					dos.writeInt(sortedRegions[0].stop- sortedRegions[0].start);
					for (int i=1; i< sortedRegions.length; i++){
						int currentStart = sortedRegions[i].start;
						//subtract 32768 to extend range of short (-32768 to 32768)
						int diff = currentStart - bp - 32768;
						dos.writeShort((short)(diff));
						dos.writeInt(sortedRegions[i].stop- sortedRegions[i].start);
						bp = currentStart;
					}
				}
				//yes short length
				else {
					//write first record's length, subtracting 32768 to extent the range of the signed short
					dos.writeShort((short)(sortedRegions[0].stop- sortedRegions[0].start - 32768));
					for (int i=1; i< sortedRegions.length; i++){
						int currentStart = sortedRegions[i].start;
						//subtract 32768 to extend range of short (-32768 to 32768)
						int diff = currentStart - bp - 32768;
						dos.writeShort((short)(diff));
						dos.writeShort((short)(sortedRegions[i].stop- sortedRegions[i].start - 32768));
						bp = currentStart;
					}
				}
			}
			//no, write int for position
			else {
				//short length? no
				if (useShortLength == false){
					//write first record's length
					dos.writeInt(sortedRegions[0].stop- sortedRegions[0].start);
					for (int i=1; i< sortedRegions.length; i++){
						int currentStart = sortedRegions[i].start;
						int diff = currentStart - bp;
						dos.writeInt(diff);
						dos.writeInt(sortedRegions[i].stop- sortedRegions[i].start);
						bp = currentStart;
					}
				}
				//yes
				else {
					//write first record's length
					dos.writeShort((short)(sortedRegions[0].stop- sortedRegions[0].start - 32768));
					for (int i=1; i< sortedRegions.length; i++){
						int currentStart = sortedRegions[i].start;
						int diff = currentStart - bp;
						dos.writeInt(diff);
						dos.writeShort((short)(sortedRegions[i].stop- sortedRegions[i].start - 32768));
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


	/**Reads a DataInputStream into this RegionData.*/
	public void read (DataInputStream dis) {
		try {
			//read text header, currently not used
			header = dis.readUTF();

			//make array
			int numberRegions = sliceInfo.getNumberRecords();
			sortedRegions = new Region[numberRegions];

			//what kind of data to follow? 
			String fileType = sliceInfo.getBinaryType();

			//int Position, int Length
			if (USeqUtilities.REGION_INT_INT.matcher(fileType).matches()){
				//make first Region, position is always an int
				int start = dis.readInt();
				sortedRegions[0] = new Region(start, start+dis.readInt());
				//read and resolve offsets to real bps and length to stop
				for (int i=1; i< numberRegions; i++){
					start = sortedRegions[i-1].start + dis.readInt();
					sortedRegions[i] = new Region(start, start + dis.readInt());	
				}
			}
			//int Position, short Length
			else if (USeqUtilities.REGION_INT_SHORT.matcher(fileType).matches()){
				//make first Region, position is always an int
				int start = dis.readInt();
				sortedRegions[0] = new Region(start, start+ dis.readShort() + 32768);
				//read and resolve offsets to real bps and length to stop
				for (int i=1; i< numberRegions; i++){
					start = sortedRegions[i-1].start + dis.readInt();
					sortedRegions[i] = new Region(start, start + dis.readShort() + 32768);	
				}
			}
			//short Postion, short Length
			else if (USeqUtilities.REGION_SHORT_SHORT.matcher(fileType).matches()){
				//make first Region, position is always an int
				int start = dis.readInt();
				sortedRegions[0] = new Region(start, start+ dis.readShort() + 32768);
				//read and resolve offsets to real bps and length to stop
				for (int i=1; i< numberRegions; i++){
					start = sortedRegions[i-1].start + dis.readShort() + 32768;
					sortedRegions[i] = new Region(start, start + dis.readShort() + 32768);	
				}
			}
			//short Position, int Length
			else if (USeqUtilities.REGION_SHORT_INT.matcher(fileType).matches()){
				//make first Region, position is always an int
				int start = dis.readInt();
				sortedRegions[0] = new Region(start, start+ dis.readInt());
				//read and resolve offsets to real bps and length to stop
				for (int i=1; i< numberRegions; i++){
					start = sortedRegions[i-1].start + dis.readShort() + 32768;
					sortedRegions[i] = new Region(start, start + dis.readInt());	
				}
			}
			//unknown!
			else {
				throw new IOException ("Incorrect file type for creating a Region[] -> '"+fileType+"' in "+binaryFile +"\n");
			}
		} catch (IOException e) {
			e.printStackTrace();
			USeqUtilities.safeClose(dis);
		}
	}

	public Region[] getRegions() {
		return sortedRegions;
	}

	public void setRegions(Region[] sortedRegions) {
		this.sortedRegions = sortedRegions;
		updateSliceInfo(sortedRegions, sliceInfo);
	}

	/**Returns whether data remains.*/
	public boolean trim(int beginningBP, int endingBP) {
		ArrayList<Region> al = new ArrayList<Region>();
		for (int i=0; i< sortedRegions.length; i++){
			if (sortedRegions[i].isContainedBy(beginningBP, endingBP)) al.add(sortedRegions[i]);
		}
		if (al.size() == 0) return false;
		sortedRegions = new Region[al.size()];
		al.toArray(sortedRegions);
		updateSliceInfo(sortedRegions, sliceInfo);
		return true;
	}
}
