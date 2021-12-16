package util.gen;





public class Delme {

	public static void main(String[] args){
		//String[] x = new String[] {"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u"};
		//String[] x = new String[] {"a","b","c","d","e","f","g","h","i","j"};
		String[] x = new String[50];
		for (int i=0; i< x.length; i++) {
			x[i] = new Integer(i).toString();
		}
		IO.pl("Len "+x.length);
		Object[][] c = chunkWithMax(x, 10);
		
		for (int i=0; i< c.length; i++) {
			Object[] o = c[i];
			for (Object oi: o) {
				IO.p(oi.toString()+" ");
			}
			IO.pl();
		}
		 
	}

	/**Splits an object[] into chunks containing a max number in each. Note this is by reference, the array is not copied. */
	public static Object[][] chunkWithMax (Object[] s, int maxInEach){
		//watch out for cases where the min can't be met
		int numChunks = s.length/maxInEach;
		if (numChunks == 0) return new Object[][]{s};

		double numLeftOver = (double)s.length % (double)maxInEach;
		//IO.pl(numChunks + " num chunks init");
		//IO.pl(numLeftOver +" num leftOver init");

		if (numLeftOver != 0) numChunks++;
		//IO.pl(numChunks + " num chunks post");

		//load with max
		int[] numInEach = new int[numChunks];
		for (int i=0; i< numChunks; i++) numInEach[i] = maxInEach;
		int numToSubtract = (numChunks*maxInEach) - s.length;
		//IO.pl(numToSubtract + " num to subtract");

		if (numToSubtract !=0) {
			int numSubtracted = 0;
			boolean loop = true;
			while (loop) {
				for (int i=0; i< numInEach.length; i++) {
					numInEach[i]--;
					numSubtracted++;
					if (numSubtracted == numToSubtract) {
						loop = false;
						break;
					}
				}
			}
		}
		//int totalObjects = 0;
		//for (int i=0; i< numInEach.length; i++) totalObjects+= numInEach[i];
		//IO.pl(totalObjects + " total objects\n");

		//build chunk array
		Object[][] chunks = new Object[numChunks][];
		int index = 0;
		//for each chunk
		for (int i=0; i< numChunks; i++){
			//create container and fill it
			Object[] sub = new Object[numInEach[i]];
			for (int j=0; j< sub.length; j++) sub[j] = s[index++];
			chunks[i] = sub;
		}
		return chunks;
	}
	
}