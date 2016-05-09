package edu.utah.seq.barcodes;

import java.util.Random;

public class BarcodeMatchSimulator {

	public static void main(String[] args) {
		int N = 7;
		double minFraction = 1;

		double total = 0;
		for (int x = 0; x< 100; x++){
			Random r = new Random();
			char[] first = fetchRandom(N, r);

			double numMatches = 0;
			for (int i=0; i< 1000000; i++){
				char[] test = fetchRandom(N, r);
				double frac = similariy (first, test);
				if (frac >= minFraction) numMatches++;
			}
			System.out.println(numMatches+"\t"+numMatches/(10000000.0));
			total+= numMatches;
		}
		double ave = (total/100);
		System.out.println(N+"\t"+minFraction+"\t1 in "+(1000000.0/ave));
	}

	private static double similariy(char[] first, char[] test) {
		double numMatches = 0;
		for (int i=0; i< first.length; i++){
			if (first[i] == test[i])numMatches++;
		}
		return numMatches/((double)first.length);
	}

	static char[] bases = {'G','A','T','C'};

	private static char[] fetchRandom(int num, Random r) {
		char[] build = new char[num];
		for (int i=0; i< num; i++){
			build[i] = bases[r.nextInt(4)];
		}
		return build;
	}

}
