package util.bio.seq;
import java.io.*;
import java.util.*;
import util.gen.*;

/**Fetches sequence given a directory of indexed sequences, see IndexFastas.
 * Start coordinates are always zero based; stop coordinate is excluded if interbase is true,
 * or included if false.*/
public class FetchIndexedSequence {

	//fields
	private File indexDirectory;
	private HashMap chromIndexedSequences;
	private IndexedSequence[] indexedSequence;
	private boolean interbase = true;

	//constructor
	public FetchIndexedSequence( File indexDirectory, boolean interbaseCoordinates){
		this.indexDirectory = indexDirectory;
		this.interbase = interbaseCoordinates;
		makeIndexedSequences();
	}

	//main methods
	/**Loads indexes for a particular chromosome. Nulls old sequences. Call before fecthSequences()*/
	public boolean setIndexedSequence(String chromosome){
		//null old chromosome sequences to same space
		if (indexedSequence != null) nullIndexedSequences();
		//fetch array of sequence objects
		if (chromIndexedSequences.containsKey(chromosome) == false) return false;
		indexedSequence = (IndexedSequence[]) chromIndexedSequences.get(chromosome);
		return true;
	}
	
	
	/**Assumes startStop[index][0 start, 1 stop] is sorted by starts.
	 * Assumes indexedSequence is set! 
	 * Must call setIndexedSequence(chromosome) prior to calling.*/
	public String[] fetchSequences (int[][] startStop){
		int i=0;
		String[] subSeqs = new String[startStop.length];
		for (int a =0; a< startStop.length; a++){
			int start = startStop[a][0];
			int stop = startStop[a][1];

			//find start and stop indexes
			boolean found =false;
			//find start
			for (; i< indexedSequence.length; i++){
				if (indexedSequence[i].contained(start)) {
					found = true;
					break;
				}
			}
			if (found == false) {
				System.err.println (indexedSequence[0].getChromosome()+" Start not found "+start+" "+stop);
				return null;
			}
			found = false;
			//find stop
			int j = i;
			for (; j< indexedSequence.length; j++){
				if (indexedSequence[j].contained(stop)) {
					found = true;
					break;
				}
			}
			if (found == false) {
				System.err.println (indexedSequence[0].getChromosome()+" End not found "+start+" "+stop);
				return null;
			}
			//get sequence
			String seq;
			//on same sequence?
			if (i == j) seq = indexedSequence[i].getSequence();
			//multiple sequences, build concat
			else {
				StringBuffer sb = new StringBuffer();
				for (int x=i; x<= j; x++) sb.append(indexedSequence[x].getSequence());
				seq = sb.toString();
			}
			//It's important to use new String()! 
			if (interbase) subSeqs[a] = new String (Seq.fetchSubSequenceInterbaseCoordinates(start, stop, indexedSequence[i].getStart(), seq));
			else subSeqs[a] = new String (Seq.fetchSubSequence(start, stop, indexedSequence[i].getStart(), seq));
		}
		return subSeqs;
	}

	/**Returns null if chromsome doesn't exist or the start or stop indexedSequence out of range.
	 * For multiple fetches use the fetchSequences() method.*/
	public String fetchSequence (String chromosome, int start, int stop){
		int[][] startStop = {
				{start, stop}
		};
		String[] seqs = fetchSequences (chromosome, startStop);
		if (seqs != null) return seqs[0];
		return null;
	}
	
	/**Returns null if chromsome doesn't exist or the start or stop indexedSequence out of range.
	 * Assumes the startStop int[index][0 start, 1 stop] is sorted by starts.*/
	public String[] fetchSequences (String chromosome, int[][] startStop){
		//fetch array of sequence objects
		if (chromIndexedSequences.containsKey(chromosome) == false) return null;
		indexedSequence = (IndexedSequence[]) chromIndexedSequences.get(chromosome);
		return fetchSequences (startStop);
	}

	public void makeIndexedSequences(){
		File[] files = IO.extractFiles(indexDirectory);

		//load hash with chromosome specific IndexedSequence objects
		chromIndexedSequences = new HashMap();
		for (int i=0; i< files.length; i++){
			String name = Seq.extractChromosomeName(files[i].getName());
			if (name == null) System.out.println("\tWARNING: cannot extract chromosome text from "+files[i]+", skipping.");
			else {
				ArrayList al;
				if (chromIndexedSequences.containsKey(name)) al = (ArrayList)chromIndexedSequences.get(name);
				else {
					al = new ArrayList();
					chromIndexedSequences.put(name, al);
				}
				al.add(new IndexedSequence(files[i]));
			}
		}

		//convert to sorted arrays
		Iterator it = chromIndexedSequences.keySet().iterator();
		while (it.hasNext()){
			Object key = it.next();
			ArrayList al = (ArrayList) chromIndexedSequences.get(key);
			IndexedSequence[] is = new IndexedSequence[al.size()];
			al.toArray(is);
			Arrays.sort(is);
			chromIndexedSequences.put(key, is);
			/*System.out.println("\n"+key);
			for (int i=0; i<indexedSequence.length; i++){
				System.out.println(i+" "+indexedSequence[i].getStart()+" "+indexedSequence[i].getSequence());
			}
			 */
		}

	}
	
	/**Good to call after each call to fetchSequences() unless you'll be fetching more from
	 * the same chromosome.*/
	public void nullIndexedSequences(){
		for (int i=0; i<indexedSequence.length; i++){
			indexedSequence[i].setSequence(null);
		}
	}
	
	public static void main (String[] args){
		FetchIndexedSequence fetch = new FetchIndexedSequence (new File ("/Users/nix/Desktop/TestIndexes"), false);
		int[][] ss = {
				{0, 10},
				{0,20},
				{10,30},
				{21,50}
		};
		String[] seqs = fetch.fetchSequences("chrX", ss);
		Misc.printArray(seqs);
		
	}

}
