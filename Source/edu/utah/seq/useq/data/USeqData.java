package edu.utah.seq.useq.data;

import java.io.*;
import edu.utah.seq.useq.*;

/**Parent Container for a sorted Data[] and it's associated SliceInfo.
 * @author david.nix@hci.utah.edu*/
public class USeqData {

	//fields
	protected SliceInfo sliceInfo;
	protected File binaryFile;
	/**Currently not used by useq archives, will be written to and read from binary file.*/
	protected String header = "";


	//methods
	public SliceInfo getSliceInfo() {
		return sliceInfo;
	}
	public File getBinaryFile() {
		return binaryFile;
	}
	public String getHeader() {
		return header;
	}
	public void setHeader(String header) {
		this.header = header;
	}
	public void setSliceInfo(SliceInfo sliceInfo) {
		this.sliceInfo = sliceInfo;
	}
	public void setBinaryFile(File binaryFile) {
		this.binaryFile = binaryFile;
	}
	/**Reads a binary file this Data object.*/
	public void read (File binaryFile) {
		FileInputStream workingFIS = null;
		DataInputStream workingDIS = null;
		try {
			//open IO
			this.binaryFile = binaryFile;
			workingFIS = new FileInputStream(binaryFile);
			workingDIS = new DataInputStream( new BufferedInputStream(workingFIS ));
			read(workingDIS);
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			USeqUtilities.safeClose(workingFIS);
			USeqUtilities.safeClose(workingDIS);
		}
	}

	/**Reads a DataInputStream into this XXXData object, to be overridden.*/
	public void read (DataInputStream dis) {}
}
