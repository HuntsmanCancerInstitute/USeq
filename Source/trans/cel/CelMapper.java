package trans.cel;
import java.io.*;
import java.util.*;
import trans.tpmap.*;
import util.gen.*;
import trans.misc.*;

/**
 * Saves float[]s of intensity values for each cel file ordered by the tpmapFeature file.
 * If useMMData is true, will double array size and set mm value afterward. 
 * 
 */
public class CelMapper {
	
	//fields
	private int numberMapFeatures;
	private boolean useMMData = true;
	private boolean printStats = true;
	private MapFeature[] features;
	private int sizeMappedFeatures;
	private HashMap chromosomeNumber = null;
	
	public CelMapper(MapFeature[] features, boolean useMMData, boolean printStats){
		//fetch tpmap MapFeature[]
		this.features = features;
		numberMapFeatures = features.length;
		this.printStats= printStats;
		this.useMMData = useMMData; 
	}
	
	public File[] mapCelFiles(File[] files){	
		File[] mappedFiles = new File[files.length];
		try{
			for (int i=0; i< files.length; i++){	
				System.out.println("\tLoading data from "+files[i]+" ...");
				
				//fetch cel file intensities
				float[][] intensities = (float[][])IO.fetchObject(files[i]);
				
				//initialize array to hold mapped intensity values
				float[] values;
				if (useMMData) values = new float[numberMapFeatures*2];
				else values = new float[numberMapFeatures];
				sizeMappedFeatures = values.length;
				
				//load values 
				System.out.println("\tMatching cel array with tpmap MapFeature array...");
				int counter = 0;				
				for (int x=0; x<numberMapFeatures; x++){
					//set PM
					values[counter++] = intensities[features[x].pmX][features[x].pmY];
					//set MM? immediately after
					if (useMMData) {
						values[counter++] = intensities[features[x].mmX][features[x].mmY];
					}
				}
				intensities = null;
				
				//save serialized float[] of mapped intensities
				File toSave = new File(files[i].getCanonicalPath()+"Map");
				IO.saveObject(toSave, values);	
				mappedFiles[i] = toSave;
				
				//calc and print stats 
				if (printStats) {
					System.out.println("\nStatistics for raw mapped cel file intensities: "+files[i]);
					Num.statFloatArray (values, true);
				}
				
			}
			return mappedFiles;
		}catch (IOException e){
			e.printStackTrace();
			return null;
		}
		
	}
	
	/**Maps float[][] to MapFeature[][] saving under the index matched File in matchedSaveAsFiles as a MultiSetQuantile[].
	 * Will perform a max(PM-MM,1) transform if useMMData is true.*/
	public boolean mapCelFiles(File[] files, File[] matchedSaveAsFiles){
		if (chromosomeNumber == null) chromosomeNumber = Util.getChromosomeNumber();
		String chromosome = null;
		try{
			for (int i=0; i< files.length; i++){	
				System.out.println("\tProcessing "+files[i]);
				
				//fetch cel file intensities
				float[][] intensities = (float[][])IO.fetchObject(files[i]);
				
				//initialize array to hold mapped intensity values
				MultiSetQuantile[] msq = new MultiSetQuantile[numberMapFeatures];
				sizeMappedFeatures = msq.length;
				
				//load values and perform subtration, (byte chromosomeNumber, int position, float value)
				int counter = 0;
				for (int x=0; x<numberMapFeatures; x++){
					chromosome = features[x].chromosome;
					byte chromNum = ((Byte)chromosomeNumber.get(chromosome)).byteValue();
					
					int position = features[x].start;
					float value;
					if (useMMData){
						value = intensities[features[x].pmX][features[x].pmY] - intensities[features[x].mmX][features[x].mmY];
						if (value < 1) value = 1;
					}
					else value = intensities[features[x].pmX][features[x].pmY];
					//set new MSQ
					msq[counter++] =  new MultiSetQuantile(chromNum, position, value);
					
				}
				intensities = null;
				
				//save serialized MultiSetQuantile[] of mapped intensities to the matched file
				IO.saveObject(matchedSaveAsFiles[i], msq);	
			}
			
			
		} catch (Exception e){
			System.out.println("\nError: is this chromosome "+chromosome+" one of the following "+chromosomeNumber);
			e.printStackTrace();
			return false;
		}
		return true;
	}
	
	public int getSizeMappedFeatures() {
		return sizeMappedFeatures;
	}
}
