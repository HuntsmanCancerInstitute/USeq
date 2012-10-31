package trans.main;
import util.gen.*;
import java.util.*;

public class RandomizeLabels {

	//fields
	private HashSet groupings = new HashSet();
	private ArrayName[] ans;
	private int numberOfTs;
	private int numberOfCs;

	//constructor
	public RandomizeLabels(float[][] t, float[][] c){
		numberOfTs = t.length;
		numberOfCs = c.length;
		makeArrayNames(t, c);
		groupings.add(groupName(ans));
	}
	
	public static void main(String[] args){
		float[][] t = {
				{0,0,0,0,0},
				{1,1,1,1,1},		
		};
		float[][] c = {
				{2,2,2,2,2},
				{3,3,3,3,3},	
		};

		RandomizeLabels rl = new RandomizeLabels(t,c);
		System.out.println(rl.concatinate()+" starting ");
		
		for (int i=0; i< 100; i++){
			float[][][] tc = rl.fetchPermutatedArrays();
			if (tc != null){
				System.out.println("T "+tc[0][0][0]+" "+tc[0][1][0] + "  C "+tc[1][0][0]+" "+tc[1][1][0]);
			}
			else break;
		}
	}
	
	/**Attempts to fetch a unique, unseen permutation of the t[][] and control[][] arrays.
	 * Returns null if no more unique permutations are possible.
	 * Otherwise returns a float[2][][] of treatment then control arrays.*/
	public float[][][] fetchPermutatedArrays(){
		ans = randomize(ans);
		if (ans == null) {
			System.out.println("No more unique permutations.");
			return null;
		}
		else {
			int counter = 0;
			//System.out.println(groupings);
			float[][] t = new float[numberOfTs][];
			for (int i=0; i<numberOfTs; i++){
				t[i]= ans[counter++].intensities;
			}
			float[][] c = new float[numberOfCs][];
			for (int i=0; i<numberOfCs; i++){
				c[i]= ans[counter++].intensities;
			}
			return new float[][][]{t,c};
		}
	}
	
	//methods
	public ArrayName[] randomize(ArrayName[] ans){
		//try up to 10000 times to get a new permutation, return first.
		for (int i=0; i< 10000; i++){
			permutate();
			String gn = groupName(ans);
			if (groupings.contains(gn) == false){
				groupings.add(gn);
				return ans;
			}
		}
		return null;
	}
	
	/**Returns a sorted groupName where the sort is on the treatments and controls, not both.*/
	public String groupName(ArrayName[] ans){
		//fetch numbers
		int[] ints = new int[ans.length];
		for (int i=0;i< ans.length; i++){
			ints[i] = ans[i].name;
		}
		Arrays.sort(ints, 0, numberOfTs);
		Arrays.sort(ints, numberOfTs, ints.length);
		return Num.intArrayToString(ints, "_");
	}
	
	/**Randomly permutates the ArrayName[].*/
	public void permutate (){
		int len = ans.length;
		ArrayName current;
		ArrayName random;
		int index;
		Random generator = new Random();
		for (int i=0; i<len; i++){
			index = generator.nextInt(len);
			current = ans[i];
			random = ans[index];
			ans[i] = random;
			ans[index] = current;
		}
	}
	
	/**Concat text numbers to a String.*/
	public String concatinate(){
		StringBuffer sb = new StringBuffer(ans[0].name+"");
		for (int i=1; i<ans.length; i++){
			sb.append("_");
			sb.append(ans[i].name);
		}
		return sb.toString();
	}
	
	/**Convert intensity arrays into an ArrayName[].*/
	public void makeArrayNames(float[][] t, float[][] c){
		int num = t.length + c.length;
		ans = new ArrayName[num];
		int counter = 0;
		//make t's
		for (int i=0; i<t.length; i++){
			ans[counter] = new ArrayName(counter, t[i]);
			counter++;
		}
		//make c's
		for (int i=0; i<c.length; i++){
			ans[counter] = new ArrayName(counter, c[i]);
			counter++;
		}
	}
	
	private class ArrayName {
		int name;
		float[] intensities;
		public ArrayName (int name, float[] intensities){
			this.name = name;
			this.intensities = intensities;
		}
	}
	
}
