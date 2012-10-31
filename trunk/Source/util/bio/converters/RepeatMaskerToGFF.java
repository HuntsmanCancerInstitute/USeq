package util.bio.converters;
import java.io.*;

import util.gen.*;

/**
Converts a RepeatMasker output to GFF annotation. Prints to screen, sorted by score
 * arg1 = fullpath text to RepeatMaskerFile, arg2(optional) = off set for start (ie 398 
 * indicate that base 1 is really 398)
 */
public class RepeatMaskerToGFF {
	public static void main(String[] args) {
		new RepeatMaskerToGFF(args);
	}
	public RepeatMaskerToGFF(String[] args){
	//read in line by line and print to this format
	//chromosome, source, feature, start, stop, score, strand,frame,attributes/comments
	String line;
	try {
		if (args.length==0) {
			System.out.println("Enter the text of a RepeatMasker.out file then optionally, a base off set to designate the real start nt\n\n e.g. java RepeatMaskerToGFF /home/nix/rmfile7.out 3958867\n\n");
			System.exit(0);
		} 
		
		BufferedReader in = new BufferedReader(new FileReader(args[0]));
		StringBuffer gff = new StringBuffer();
		int offSet=0;
		if (args.length>1) offSet = Integer.parseInt(args[1])-1; //subtract 1 from given offset
		while ((line = in.readLine()) !=null) {
			line = line.trim();                        //kill whitespace and test if exists
			if (line.length() == 0) continue;          //skip blank lines
			if (line.matches("^\\d.+")==false) continue; //make sure line begins with a number
			String[] items = line.split("\\s+");
			//build gff line and print
			StringBuffer sb = new StringBuffer("unknown\tRepeatMasker\t");
			//sb.append(items[10]+"\t");	//feature
			sb.append("repeat\t");
			sb.append((Integer.parseInt(items[5])+offSet)+"\t");	//start
			sb.append((Integer.parseInt(items[6])+offSet)+"\t");	//stop
			sb.append(items[0]+"\t");	//score
			if (items[8].equals("+")) sb.append(items[8]+"\t");  //strand
			else sb.append("-\t");	//strand
			sb.append(".\t");	//frame
			sb.append("text="+items[9]+"; class="+items[10]+"; percDiv="+items[1]+"; percDel="+items[2]+
				"; percIns="+items[3]+"; querySeq="+items[4]+"; beginInRepeat="+items[11]+
				"; endInRepeat="+items[12]);	//attributes
			System.out.println(sb);
			gff.append(sb);
			gff.append("\n");
			
		}
		in.close();
		IO.writeString(gff.toString(), args[0]+".gff");
		System.out.println("Wrote file: "+args[0]+".gff");
	}
	catch (IOException e) {
		e.printStackTrace();
	}
	
	}
	/**Prints out the text array on a line divided by the separator, good for debugging*/
	public static void printStringArray(String[] array, String separator){
		int len = array.length;
		for (int i=0; i<len; i++){
			System.out.print (array[i]+separator);
		}
		System.out.println();
	}
}
