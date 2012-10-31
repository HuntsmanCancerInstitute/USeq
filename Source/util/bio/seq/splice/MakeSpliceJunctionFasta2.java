package util.bio.seq.splice;

import java.util.Arrays;

public class MakeSpliceJunctionFasta2 {


	public MakeSpliceJunctionFasta2(){
		//start, end, index
		Exon[] exons = new Exon[4];
		exons[0] = new Exon(0,10,0);
		exons[1] = new Exon(20,30,1);
		exons[2] = new Exon(40,50,2);
		exons[3] = new Exon(60,70,3);
		
		int maxNumToSkip = exons.length-1;
		
		/*for (int i=0; i< maxNumToSkip; i++){
			System.out.println("NumToSkip "+i);
			//iterate through exons making indexes to skip
			for (int j=0; j< exons.length; j++){
				int[] indexesToSkip = new int[i];
				
			}
		}*/
		int[][] arrays;
		//skip zero
		/*System.out.println("Skip zero");
		Exon[] p = fetchExonConcatinate(exons, new int[]{});
		printExonIndexes(p);
		
		//skip one
		//System.out.println("\nSkip one");
		arrays = skipOne(exons);
		for (int i=0; i< arrays.length; i++){
			Exon[] parsed = fetchExonConcatinate(exons, arrays[i]);
			printExonIndexes(parsed);
		}
		
		//skip two
		System.out.println("\nSkip two");
		arrays = skipTwo(exons);
		for (int i=0; i< arrays.length; i++){
			Exon[] parsed = fetchExonConcatinate(exons, arrays[i]);
			printExonIndexes(parsed);
		}*/
		
		System.out.println("\nSkip three");
		arrays = skipThree(exons);
		for (int i=0; i< arrays.length; i++){
			Exon[] parsed = fetchExonConcatinate(exons, arrays[i]);
			printExonIndexes(parsed);
		}
		
		
	}
	public int[][] skipThree(Exon[] exons){
		int num = (exons.length-2) *3;
		int index = 0;
		int[][] toSkip = new int[num][3];
		for (int i=0; i< exons.length; i++){
			for (int j=i+1; j< exons.length; j++){
				for (int k=j+1; k< exons.length; k++){
				toSkip[index++] = new int[]{i,j,k};
				System.out.println("ToSkipThree "+i+" "+j+" "+k);
				}
			}
		}
		return toSkip;
	}
	
	public int[][] skipTwo(Exon[] exons){
		int num = (exons.length-1) *2;
		int index = 0;
		int[][] toSkip = new int[num][2];
		for (int i=0; i< exons.length; i++){
			for (int j=i+1; j< exons.length; j++){
				toSkip[index++] = new int[]{i,j};
				//System.out.println("ToSkipTwo "+i+" "+j);
			}
		}
		return toSkip;
	}
	
	public int[][] skipOne(Exon[] exons){
		int[][] toSkip = new int[exons.length][1];
		for (int i=0; i< exons.length; i++){
			toSkip[i] = new int[]{i};
		}
		return toSkip;
	}
	
	public Exon[] fetchExonConcatinate(Exon[] exons, int[] indexesToSkip){
		//System.out.println("FetchExons "+Arrays.toString(indexesToSkip));
		//make reference copy
		Exon[] clone = new Exon[exons.length];
		for (int i=0; i< exons.length; i++)  clone[i] = exons[i];
		//null references in clone to skip
		for (int i=0; i< indexesToSkip.length; i++) clone[indexesToSkip[i]] = null;
		//make new array
		Exon[] parsed = new Exon[exons.length - indexesToSkip.length];
		int index =0;
		for (int i=0; i< clone.length; i++){
			if (clone[i] != null) parsed[index++] = clone[i];
		}
		return parsed;
	}
	
	public void printExonIndexes(Exon[] exons){
		if (exons.length ==0) System.out.println("No exons!");
		System.out.print(exons[0].index);
		for (int i=1; i< exons.length; i++) System.out.print("_"+ exons[i].index);
		System.out.println();
	}
	
	
	
	
	public static void main(String[] args) {
		new MakeSpliceJunctionFasta2();

	}

}
