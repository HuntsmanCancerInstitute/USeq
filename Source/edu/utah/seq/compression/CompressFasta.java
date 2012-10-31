package edu.utah.seq.compression;
import java.util.zip.*;
import java.util.*;
import util.bio.seq.*;
import util.gen.*;
import java.io.*;

/**Application that compresses or decompresses a multi fasta file.  
 * Names of the individual reads will be converted to 1,2,3...
 * For compression, each of the reads must be 1) the same length 2) on one line 3) only contain gatc characters.
 * Thus this app is best suited for compressing signature sequencer reads.
 * It has a 20% better compression rate that just zipping, gzipping, or bzip2ing the text fasta.*/
public class CompressFasta {
	//fields
	private HashMap seqByte = Seq.makeByte4BaseMap();
	private String[] byteSeq = Seq.all4BaseCombinations;
	private DataOutputStream out;
	private int size = 0;
	private int numberOfSequences;
	private static final String nameVersion = "CompressedFasta1.0";

	public static void main(String[] args) {
		if (args.length ==0)Misc.printExit("\nEnter the file or directory containing same length, single line, DNA sequence xxx.fasta file(s). GATC only.\n");
		new CompressFasta(IO.extractFiles(new File(args[0]), "fasta"));
		//new CompressFasta(new File(""))
	}

	public void decompressIt(File f){

	}

	public CompressFasta(File[] files){
		for (int i=0; i< files.length; i++){
			String name = IO.getFullPathName(files[i]);
			File compressed = new File (name+".binary");
			String notes = "Original file: "+name+" date last modified: "+new Date (files[i].lastModified());
			System.out.print("\nCompressing: "+files[i].getName()+" \t");
			compressIt(files[i],compressed, notes);
			System.out.println("\nDecompressing....");
			printBinarySequence(new File (compressed+".zip"), new File (compressed.getParentFile(),"decom.txt"));
		}
	}

	public boolean findSeqSizeAndNumberReads(File fastaFile){
		try {
			BufferedReader in = new BufferedReader (new FileReader(fastaFile));
			String line;
			//reset size and number of sequences
			size = 0;
			numberOfSequences = 0;
			//calculate size based on first seq
			while ((line = in.readLine()) !=null){
				//advance to header line
				if (line.startsWith(">") == false) continue;
				line = in.readLine().trim();
				//set size of read
				size = line.length();
				//increment number of seqs
				numberOfSequences++;
				break;
			}
			//count number of seqs
			while ((line = in.readLine()) !=null){
				if (line.startsWith(">")) numberOfSequences++;
				else {
					int length = line.trim().length();
					if (length != size && length !=0){
						System.err.println("\n\nSequence sizes differ, all must be the same length ("+size+"), aborting, see '"+line+"' in "+fastaFile+"\n");
						return false;
					}
				}
			}
			in.close();
			//report
			if (size == -1 ) {
				System.err.println("\n\nFailed to determine size of first sequence, check format, aborting. "+fastaFile+"\n");
				return false;
			}
			else if (numberOfSequences == 0 ) {
				System.err.println("\n\nNo sequences found, check format, aborting. "+fastaFile+"\n");
				return false;
			}
			else System.out.println("Size: "+size+" nts\t Number: "+numberOfSequences);
			return true;
		} catch (Exception e) {
			System.err.println("\n\nFailed to determine size and number of sequences, aborting. "+fastaFile+"\n");
			e.printStackTrace();
			return false;
		}
	}

	public boolean compressIt(File fastaFile, File compressedFile, String notes){
		//stat fasta file
		if (findSeqSizeAndNumberReads(fastaFile) == false) return false;
		try {
			BufferedReader in = new BufferedReader (new FileReader(fastaFile));
			out = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(compressedFile)));
			//write format text, size of individual sequence and number of counts
			out.writeUTF(nameVersion);
			out.writeUTF(notes);
			out.writeInt(size);
			out.writeInt(numberOfSequences);
			//run through fasta
			int numSeqs = 0;
			String line;
			StringBuffer sb = new StringBuffer();
			while ((line = in.readLine()) !=null){
				line = line.trim();
				int lineLength = line.length();
				if (line.startsWith(">") || lineLength == 0) continue;
				numSeqs++;
				line = line.toLowerCase();
				sb.append(line);
				int remainder = lineLength % 4;
				if (remainder == 0) {
					writeBytes(sb.toString());
					sb = new StringBuffer();
				}
			}
			//write final seq?
			if (sb.length() !=0) writeFinalBytes (sb.toString());
			out.close();
			in.close();
			//check number of sequences
			if (numSeqs != numberOfSequences) {
				System.err.println("\n\nIncorrect number of sequences found, aborting. "+fastaFile+"\n");
				compressedFile.delete();
				return false;
			}
			//zip it
			if (IO.zipAndDelete(compressedFile) == false){
				System.err.println("\n\nZip compression failed, aborting. "+fastaFile+"\n");
				compressedFile.delete();
				return false;

			}
			return true;
		} catch (Exception e) {
			System.err.print("\n\nFailed to compress fasta file, aborting. "+fastaFile+"\n");
			e.printStackTrace();
			compressedFile.delete();
			return false;
		}
	}

	private void writeBytes(String seq) throws Exception{
		int lengthSeq = seq.length();
		//write 4 bp bytes
		for (int i=0; i< lengthSeq; i+=4){
			String sub = seq.substring(i,i+4);
			Object o = seqByte.get(sub);
			if (o == null){
				throw new Exception("\n\nError converting sequence to byte, check for non gatc characters in -> "+seq+"\n");
			}
			out.writeByte(((Byte)o).byteValue());
		}
	}

	private void writeFinalBytes(String seq) throws Exception{
		int lengthSeq = seq.length();
		//append c's to stop to bring up to a multiple of four
		int basesToAdd = 4 - lengthSeq % 4;			
		if (basesToAdd != 4){
			if (basesToAdd == 1) seq = seq + "c";
			else if (basesToAdd == 2) seq = seq + "cc";
			else seq = seq + "ccc";
		}
		writeBytes(seq);
	}
	
	public boolean printBinarySequence(File compressed, File decompressed){
		try {
			DataInputStream in;
			//unzip it?
			if (compressed.getName().endsWith(".zip")){
				ZipFile zf = new ZipFile(compressed);
				ZipEntry ze = (ZipEntry) zf.entries().nextElement();
				in = new DataInputStream( new BufferedInputStream(zf.getInputStream(ze)));
			}
			else { 
				in = new DataInputStream( new BufferedInputStream (new FileInputStream(compressed) ));
			}
			//read array text
			String fileType = in.readUTF();
			if (fileType.equals(nameVersion) == false) {
				System.err.println("\nCannot read in binary file, wrong file type. Looking for a "+nameVersion+" file, found a "+fileType+"\n");
				return false;
			}
			//make writer
			PrintWriter out = new PrintWriter (new FileWriter(decompressed));
			//read notes, sequence size, and number
			String notes = in.readUTF();
			System.out.println("\t"+notes);
			int size = in.readInt();
			int number = in.readInt();
			int numberPlus4 = number +4;
			boolean partial = size % 4 !=0;
			//read in and print seqs 
			int count = 1;
			StringBuffer sb = new StringBuffer(number+4);
			while (count <= number){
				//load in bytes
				while (sb.length()< size){
					byte b = in.readByte();
					sb.append(byteSeq[b + 128]);
				}
				//print
				out.println(">"+count);
				out.println(sb.substring(0, size));
				//reset StringBuffer
				if (partial) sb = new StringBuffer(sb.subSequence(size, sb.length()));
				else sb = new StringBuffer(numberPlus4);
				count++;
			}
			out.close();
			in.close();
			return true;
		}
		catch (Exception ioe){
			System.err.println("\nProblem uncompressing "+compressed);
			ioe.printStackTrace();
			return false;   
		} 
	}

}
