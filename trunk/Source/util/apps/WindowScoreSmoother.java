package util.apps;
import java.io.*;
import util.gen.*;

/**Slides a window across an array of numbers averaging the scores.*/
public class WindowScoreSmoother {

	public static void main(String[] args) {
		if (args.length !=2)Misc.printExit("\nEnter a full path file text for a text file containing a column of numbers and the window size to scan.\n");
		
		//load double[]
		double[] scores = Num.loadDoubles(new File(args[0]));
		//parse window size
		int window = Integer.parseInt(args[1]);
		
		//run window across scores printing average
		int counter = 0;
		while (true){
			int end = window + counter;
			double total = 0;
			for (int i=counter; i< end; i++){
				if (i == scores.length)System.exit(0);
				total += scores[i];
			}
			counter++;
			System.out.println(total/window);
		}
		
	}

	

	
	

	}


