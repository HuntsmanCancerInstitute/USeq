package trans.cel;
import util.gen.*;
import trans.tpmap.*;
import java.io.*;
import java.util.*;

/**Converts a tsv file to xxx.celp files (float[]s) linked to the tpmap.
 * Very memory intensive. Values from tsv file are antiLog2ed.  Missing values set to 1.*/
public class ConvertTVSToCelps {

	//fields
	private MapFeature[] mapFeatures;
	private File tvs;
	private HashMap tvsHash;
	private int numCelpFiles;
	private String[] headerLine;
	
	public ConvertTVSToCelps(String[] args) {
		//collect args
		tvs = new File(args[0]);
		File tpmapDir = new File(args[1]);
		File featureFile = new File(tpmapDir, "tpmap.fa");
		
		//load MapFeature[]
		mapFeatures = (MapFeature[])IO.fetchObject(featureFile);
		
		//load tvs hash
		loadTVSHash();
		
		//make float[]s
		float[][] celps = new float[numCelpFiles][mapFeatures.length];
		
		//for each MapFeature line find associated tvs value
		int numMatches = 0;
		int numMisses = 0;
		for (int x=0; x< mapFeatures.length; x++){
			String key = mapFeatures[x].chromosome+"_"+mapFeatures[x].start;
			if (tvsHash.containsKey(key)){
				numMatches++;
				float[] values = (float[])tvsHash.get(key);
				//load celps
				for (int y =0; y<numCelpFiles; y++){
					celps[y][x] = values[y];
				}
			}
			else{
				numMisses++;
				//load celps with ones
				for (int y =0; y<numCelpFiles; y++){
					celps[y][x] = 1;
				}
				
			}
		}
		System.out.println(numMatches+" num matches  "+numMisses+" num misses  "+mapFeatures.length+" total");
		
		//save celps
		int index = 0;
		for (int i=2; i<headerLine.length; i++){
			File celpFile = new File (tvs.getParentFile(), headerLine[i]+"p");
			IO.saveObject(celpFile, celps[index]);
			//make histogram
			System.out.println("Making histogram for "+celpFile);
			Histogram his = new Histogram(0, 100, 100);
			his.countAll(celps[index]);
			his.printScaledHistogram();
			
			
			index++;
			
		}
		
		System.out.println("\nDone\n");
	}
	
	public static void main(String[] args) {
		new ConvertTVSToCelps(args);
	}
	
	public void loadTVSHash(){
		tvsHash = new HashMap(3000000);
		try {
			BufferedReader in = new BufferedReader( new FileReader(tvs));
			//read in and save header
			String line = in.readLine();
			headerLine = line.split("\\s+");
			//run through adding to hash
			while ((line=in.readLine()) !=null){
				String[] tokens = line.split("\\s+");
				if (tokens.length < 3) continue;
				String key = tokens[0]+"_"+tokens[1];
				numCelpFiles = tokens.length - 2;
				float[] values = new float[numCelpFiles];
				for (int i=0; i< numCelpFiles; i++){
					double value = Double.parseDouble(tokens[i+2]);
					values[i] = new Double(Math.pow(2,value)).floatValue();
				}
				if (tvsHash.containsKey(key)) Misc.printExit("\nDuplicate key "+key);
				tvsHash.put(key, values);
			}
			in.close();	
		} catch (IOException e){
			e.printStackTrace();
		}
		
	}
	
	

}