package trans.misc;
import java.io.*;

import util.gen.*;

/** Slides a 50 index window along an array averaging the contents. */
public class AverageDoubles {

	
	public static void main(String[] args){
		double[] dx = Num.loadDoubles(new File(args[0]));
		
		double[] aves = Num.windowAverageScores(dx, 50);
		
		Misc.printArray(aves);
	}

}
