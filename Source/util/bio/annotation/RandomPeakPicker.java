package util.bio.annotation;

import java.io.File;
import java.util.Random;

import edu.utah.hci.misc.Gzipper;

public class RandomPeakPicker {

	public static void main( String[] args) {
		try {
			int numToCreate = 5000;
			
			Bed[] chrEnd = Bed.parseFile(new File(args[0]), 0, 0);
			int numChrs = chrEnd.length;

			Random rd = new Random(System.currentTimeMillis());
			Gzipper out = new Gzipper(new File(args[1]));
			
			//for each random region to create
			for (int x =0; x< numToCreate; x++) {
				//pick a chrom
				int chromIndex = rd.nextInt(numChrs);
				Bed region = chrEnd[chromIndex];
				
				//pick a point within the region
				int len = region.getLength();
				int pos = rd.nextInt(len)+ region.getStart();
				
				//write out bed line 
				boolean strand = rd.nextBoolean();
				String s = null;
				if (strand) s= "+";
				else s="-";
				out.println(region.getChromosome()+"\t"+pos+"\t"+(pos+1)+"\tR\t0\t"+s);
				
			}
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}


	}


}


