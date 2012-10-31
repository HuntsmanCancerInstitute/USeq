package trans.bkgrnd;
import java.io.*;
import java.util.*;
import trans.misc.*;
import trans.graphics.*;
import util.gen.*;

/**
 * Loads a bpmap file with oligo intensity measurements from multiple '.cel' files after median scaling to 50.
 * */
public class MapControlIntensities {
	public static final float TARGET_MEDIAN = 50;

	public static void main(String[] args) {
		//process argws
		if (args.length==0){
			System.out.println("\nbpmap celDirectoryOrFile");
			System.exit(0);
		}
		File bpmapFile = new File (args[0]);
		File celFileDir = new File (args[1]);
		File[] celFiles = IO.extractFiles(celFileDir, ".cel");
		
		//parse bpmap
		System.out.println("Parsing bpmap...");
		BackGroundOligo[] oligos = makeInitialBGOs (bpmapFile);
		
		//load bkoligos with 	median normalized intensities
		loadIntensities(oligos, celFiles);
		
		//print out oligos (seq gc tm 23,36,19;45,46,41;34,38)
		for (int i=0; i< oligos.length; i++){
			System.out.println(oligos[i]);
		}		
		//save serialized oligos
	}
	
	public static void loadIntensities (BackGroundOligo[] oligos, File[] celFiles){
		for (int i=0; i< celFiles.length; i++){
			System.out.println("Loading "+celFiles[i]);
			//load virtual slide with .cel file intensities
			float[][] intensities = Util.createVirtualCel (celFiles[i]);
			//median normalize
			Num.medianNormalize(intensities, TARGET_MEDIAN);
			//run thru oligos adding intensities
			for (int j =0; j<oligos.length; j++){
				//get the coordinates
				int[][] coor = oligos[j].getCoordinates();
				float[] vals = new float[coor.length];
				//load values array
				for (int k=0; k< coor.length; k++){
					int[] xy = coor[k];
					vals[k] = intensities[xy[0]][xy[1]];
				}
				//set values array
				oligos[j].getIntensityArrayLists().add(vals);
			}
		}
	}
	
	public static BackGroundOligo[] makeInitialBGOs(File bpmap){
		BackGroundOligo[] bkOligos = null;
		try{
			//read in bpmap file to hashMap
			System.out.println("\tReading...");			
			BufferedReader in = new BufferedReader(new FileReader (bpmap));
			String[] tokens;
			String line;
			HashMap map = new HashMap(10000);
			while ( (line=in.readLine()) !=null){
				tokens = line.split("\\t");
				//check if it exists
				if (map.containsKey(tokens[0])){
					//get Set and add new String x,y
					HashSet al = (HashSet)map.get(tokens[0]);
					al.add(tokens[4]+","+tokens[5]);
				}
				else{
					HashSet al = new HashSet();
					al.add(tokens[4]+","+tokens[5]);
					map.put(tokens[0], al);
				}
			}
			in.close();
			
			//convert HashMap seq:HashSet[String(x,y)] to BackGroundOligo
System.out.println("\tConverting...");			
			bkOligos = new BackGroundOligo[map.size()];
			int oligoCounter = 0;
			Iterator it = map.keySet().iterator();
			while (it.hasNext()){
				String seq = (String)it.next();
				//make int[][] of x,y coordinates
				HashSet coors = (HashSet)map.get(seq);
				int[][] coordinates = new int[coors.size()][2];
				String[] xy;
				Iterator it2 = coors.iterator();
				int counter =0;
				while (it2.hasNext()){
					String s = ((String)it2.next());
					xy = s.split(",");
					coordinates[counter] = new int[] { Integer.parseInt(xy[0]), Integer.parseInt(xy[1]) };
					counter++;
				}
				//make new oligo
				bkOligos[oligoCounter] = new BackGroundOligo(seq, coordinates);
				oligoCounter++;
			}
		}catch (IOException e){
			e.printStackTrace();
		}
		return bkOligos;
	}
	
	
	
	
	
	
	
}


