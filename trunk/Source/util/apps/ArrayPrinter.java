package util.apps;
import java.io.*;
import java.util.*;

import util.gen.IO;



/**
 * Prints serialized arrays to screen.
 * @author nix
 *
 */
public class ArrayPrinter {
	public static void main(String[] args) {
		if (args.length!=3){
			System.out.println ("\nTo print a serialized array to screen, enter:" +
				"\n\t1) a number for the start index (enter 0 for the first index)," +
				"\n\t2) a number for the stop index (enter 0 for the last index)," +
				"\n\t3) the text of the serialized array file." +
				"\n(ie java nix/util/gen/ArrayReader 0 9 UVZeste.normTns)\n");
			System.exit(0);
		}

		Object vals = IO.fetchObject(new File(args[2]));
		Object valObj = vals.getClass();
		
		float[] floatArray = new float[0];
		short[] shortArray = new short[0];
		int[] intArray = new int[0];
		byte[] byteArray = new byte[0];
		double[] doubleArray = new double[0];
		ArrayList al = new ArrayList();
		int[][] intIntArray = new int[0][0];
		
		int start = Integer.parseInt(args[0]);
		int end = Integer.parseInt(args[1]);
		
		if (floatArray.getClass().equals(valObj)){
			floatArray =(float[])vals;
			if (end==0 || end>floatArray.length) end = floatArray.length;
			System.out.println("\n"+floatArray.length+" values.  Printing float[] from "+start+" to "+end);
			for (int i=start; i<end; i++)System.out.println(floatArray[i]);
		}
		else if(shortArray.getClass().equals(valObj)){
			shortArray =(short[])vals;
			if (end==0 || end>shortArray.length) end = shortArray.length;
			System.out.println("\n"+shortArray.length+" values.  Printing short[] from "+start+" to "+end);
			for (int i=start; i<end; i++)System.out.println(shortArray[i]);
		}
		else if (intArray.getClass().equals(valObj)) {
			intArray =(int[])vals;
			if (end==0 || end>intArray.length) end = intArray.length;
			System.out.println("\n"+intArray.length+" values.  Printing int[] from "+start+" to "+end);
			for (int i=start; i<end; i++)System.out.println(intArray[i]);
		}
		else if (byteArray.getClass().equals(valObj)) {
			byteArray =(byte[])vals;
			if (end==0 || end>byteArray.length) end = byteArray.length;
			System.out.println("\n"+byteArray.length+" values.  Printing byte[] from "+start+" to "+end);
			for (int i=start; i<end; i++)System.out.println(byteArray[i]);
		}
		else if (doubleArray.getClass().equals(valObj)) {
			doubleArray =(double[])vals;
			if (end==0 || end>doubleArray.length) end = doubleArray.length;
			System.out.println("\n"+doubleArray.length+" values.  Printing double[] from "+start+" to "+end);
			for (int i=start; i<end; i++)System.out.println(doubleArray[i]);
		}
		else if (al.getClass().equals(valObj)){
			al = (ArrayList)vals;
			if (end==0 || end>al.size()) end = al.size();
			System.out.println("\n"+al.size()+" values.  Printing ArrayList from "+start+" to "+end);
			for (int i=start; i<end; i++)System.out.println(al.get(i));			
		}
		else if (intIntArray.getClass().equals(valObj)){
			intIntArray = (int[][])vals;
			if (end==0 || end>intIntArray.length) end = intIntArray.length;
			System.out.println("\n"+intIntArray.length+" values.  Printing int[][] from "+start+" to "+end);
			int len;
			for (int i=start; i<end; i++){
				len = intIntArray[i].length;
				for (int j=0; j<len; j++){
					System.out.println("int["+i+"]"+"["+j+"]= "+intIntArray[i][j]);
				}
			}			
		}	
		//attempt to print object[]	ala toString()
		else {
			try{
				Object[] objs = (Object[])vals;
				if (end==0 || end>objs.length) end = objs.length;
				System.out.println("\n"+objs.length+" values.  Printing object[] from "+start+" to "+end);
				for (int i=start; i<end; i++)System.out.println(objs[i]);
			}catch (Exception e){
				System.out.println ("Sorry, can't seem to read this serialize object file?!");
			}
		}
		
		System.out.println();
	}

}
