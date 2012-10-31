package trans.qc;
import java.io.*;
import util.gen.*;
import java.util.*;

public class DestypeIntensityCalculator {

	public static void main(String[] args){
		/*
		 * ndb[0] = zeroA;
			ndb[1] = threeA;
			ndb[2] = negThreeA;
			ndb[3] = sixA;
			ndb[4] = negSixA;
			ndb[5] = pm;
		 */
		int[][][] ccs = (int[][][])IO.fetchObject(new File(args[0]));
		
		float[][] cela = (float[][])IO.fetchObject(new File(args[1]));
		
		for (int i=0; i< ccs.length; i++){
			System.out.println("\nNum CCs "+ccs[i].length);
			float[] ints = Num.fetchMatrixValues(cela, ccs[i]);
			Arrays.sort(ints);
			System.out.println("Mean "+Num.mean(ints));
			System.out.println("Median "+Num.median(ints));
		}
		
	}
	
}
