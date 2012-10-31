package util.bio.parsers;
import java.io.*;
import java.util.*;
import util.gen.*;

public class ADFDataBaseSplitter {



	public static void main (String[] args){
		File column = new File (args[0]);
		try{
			//how many different dbs are present
			LinkedHashMap dbs = new LinkedHashMap();
			BufferedReader in = new BufferedReader ( new FileReader (column));
			String line;
			int index = 0;
			while ((line = in.readLine()) != null){
				if (line.trim().length() == 0) continue;
				String[] tokens = line.split("\\|");
				for (int i=0; i< tokens.length; i+=2){
					if (dbs.containsKey(tokens[i]) == false){
						dbs.put(tokens[i], new Integer (index++));
					}
				}
			}
			in.close();
			//print header
			Iterator it = dbs.keySet().iterator();
			while (it.hasNext()){
				String c = (String)it.next();
				System.out.print("Reporter BioSequence DatabaseEntry ["+c+"]\t");
			}
			System.out.println();
			//print lines
			in = new BufferedReader ( new FileReader (column));
			int numColumns = dbs.size();
			while ((line = in.readLine()) != null){
				if (line.trim().length() == 0) {
					System.out.println();
				}
				else {
					String[] columns = new String[numColumns];
					String[] tokens = line.split("\\|");
					for (int i=0; i< tokens.length; i+=2){
						//System.out.println(line+"Tok"+tokens[i]);
						int columnIndex = ((Integer)dbs.get(tokens[i])).intValue();
						if (columns[columnIndex]==null) columns[columnIndex]= tokens[i+1];
						else columns[columnIndex]= columns[columnIndex] +";"+tokens[i+1];
					}
					//print it
					for (int i=0; i< numColumns; i++){
						if (columns[i] == null) columns[i]="";
					}
					System.out.println(Misc.stringArrayToString(columns, "\t"));
				}
			}
			in.close();

		}catch (Exception e){
			e.printStackTrace();
		}



	}





}
