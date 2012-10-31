package util.bio.seq;
import java.util.HashMap;
import util.gen.Misc;
import java.io.*;

/**Can read or write a compressed boolean[] by converting
 * each 8 booleans to a unique text and assigns them to byte.
 * Reads and writes the bytes.*/
public class BinaryBooleanReaderWriter {

	//fields
	private HashMap <String,Byte> trueFalseByte = null;
	private boolean[][] byteBoolean = null;
	private int sizeOriginalBoolean;
	private byte[] bytes;
	private boolean[] booleans;
	private File file;
	
	//constructor
	public BinaryBooleanReaderWriter(){
		makeMaps();
	}
	
	//methods
	/**Saves a boolean[] in a compressed binary format.
	 * @return boolean - true if successful.*/
	public boolean saveBinaryBooleanArray (boolean[] toSave, File file){
		sizeOriginalBoolean = toSave.length;
		booleans = appendToEight (toSave);
		bytes = convert2Bytes (booleans);
		if (booleans.length != sizeOriginalBoolean) booleans = null;
		this.file = file;
		return writeBytes();
	}
	
	/**Loads the boolean[] booleans field from a compressed binary format.
	 * @return boolean[] or null if problem is encountered.*/
	public boolean[] loadBinaryBooleanArray (File file){
		this.file = file;
		return readBytesAndConvertToBoolean();
	}
	
	private boolean writeBytes(){
		try {
			DataOutputStream out = new DataOutputStream(new BufferedOutputStream( new FileOutputStream(file)));
			//write size of original boolean array
			out.writeInt(sizeOriginalBoolean);
			//write out size of byte array
			out.writeInt(bytes.length);
			//write bytes
			out.write(bytes);
			out.close();
			//clean up
			bytes = null;
			return true;
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}
	}
	
	private boolean[] readBytesAndConvertToBoolean(){
		try {
			DataInputStream in = new DataInputStream( new BufferedInputStream( new FileInputStream(file)) );
			//read in size of originalBoolean array
			sizeOriginalBoolean = in.readInt();
			//read in size of byte array
			int sizeByteArray = in.readInt();
			//make boolean array
			boolean[] loadedBoolean = new boolean[sizeByteArray * 8];
			int index = 0;
			for (int i=0; i< sizeByteArray; i++){
				int num = in.readByte() + 128;
				boolean[] toAppend = byteBoolean[num];
				System.arraycopy(toAppend, 0, loadedBoolean, index, 8);
				index += 8;
			}
			//trim array?
			if (sizeOriginalBoolean == loadedBoolean.length) booleans = loadedBoolean;
			else {
				booleans = new boolean[sizeOriginalBoolean];
				System.arraycopy(loadedBoolean, 0, booleans, 0, sizeOriginalBoolean);
				loadedBoolean = null;
			}
			in.close();
			return booleans;
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}
	}

	/**Converts a boolean[] to a 0=false,1=true String.*/
	public static String convert2String (boolean[] b, int startIndex, int stopIndex_Excluded){
		char[] c = new char[stopIndex_Excluded-startIndex];
		int counter = 0;
		for (int i= startIndex; i< stopIndex_Excluded; i++){
			if (b[i] == true) c[counter++] = '1';
			else c[counter++] = '0';
		}
		return new String(c);
	}

	public boolean[] convert2Boolean (byte[] bytes){
		//make the array
		boolean[] bArray = new boolean[bytes.length * 8];
		int index =0;
		//for each byte 
		for (int i=0; i< bytes.length; i++){
			int num = bytes[i] + 128;
			System.arraycopy(byteBoolean[num], 0, bArray, index, 8);
			index += 8;
		}
		return bArray;
	}

	/**Assumes the boolean[] is evenly divisable by 8.*/
	public byte[] convert2Bytes (boolean[] booleans){

		//make byte[]
		byte[] converted = new byte[booleans.length/8];
		//for each 8 booleans
		int index =0;
		for (int i=0; i< booleans.length; i+=8){
			String stringRep = convert2String(booleans, i, i+8);
			converted[index++] = trueFalseByte.get(stringRep).byteValue();
		}
		return converted;
	}

	/**Copies and increases the size of the boolean[] to be divisible by 8.
	 * Returns the original array if it is already divisible by 8.*/
	public static boolean[] appendToEight(boolean[] b){
		int remainder = b.length % 8;
		if (remainder == 0) return b;
		int numToAppend = 8-remainder;
		//make new array
		boolean[] appended = new boolean[b.length+numToAppend];
		System.arraycopy(b, 0, appended, 0, b.length);
		return appended;
	}
	
	public void makeMaps(){
		//make the converter
		//1= true, 0= false
		byteBoolean = new boolean[256][8];
		boolean[] fT = {false, true};
		int index = 0;
		trueFalseByte = new HashMap();
		byte counter = -128;
		for (int a=0; a<2; a++){
			for (int b=0; b<2; b++){
				for (int c=0; c< 2; c++){
					for (int d=0; d<2; d++){
						for (int e=0; e<2; e++){
							for (int f=0; f<2; f++){
								for (int g=0; g<2; g++){
									for (int h=0; h<2; h++){
										trueFalseByte.put(a+""+b+""+c+""+d+""+e+""+f+""+g+""+h, new Byte(counter++));
										byteBoolean[index++] = new boolean[]{fT[a], fT[b], fT[c], fT[d], fT[e], fT[f], fT[g], fT[h]};
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	//main for testing
	/*public static void main(String[] args) {
		boolean[] b = {true,true,true,false,true,true,true,true,false};
		File save = new File ("deleteMe.bb");
		//save it
		BinaryBooleanReaderWriter bb = new BinaryBooleanReaderWriter();
		bb.saveBinaryBooleanArray(b, save);
		bb.loadBinaryBooleanArray(save);
		Misc.printArray(bb.getBooleans());
		save.deleteOnExit();
	}*/
	
	

	public boolean[] getBooleans() {
		return booleans;
	}

}
