package util.apps;
import java.util.*;
import java.util.regex.*;
import java.io.*;
import util.gen.*;

public class ParseMLSHTML {

	//fields
	private File htmlFile;
	private String[][] listings;
	private HouseBean[] houseBeans;
	
	public ParseMLSHTML(String[] args){
		htmlFile = new File (args[0]);
		System.out.print("\nParsing html ");
		parseListingsFile();
		System.out.print(listings.length+"\nMaking house beans ");
		makeHouseBeans();
		System.out.print(houseBeans.length+"\n");
		printForGeoCode();
	}
	
	public void printForGeoCode(){
		//print header
		System.out.println("Address\tCity\tState\tZipcode\tName\tURL\tColor\tImage\tGroup\tMisc");
		//print each bean
		for (int i=0; i< houseBeans.length; i++){
			if (houseBeans[i].getGarage() == 0 || houseBeans[i].getSize() < 2000) continue;
			StringBuffer sb = new StringBuffer();
			sb.append(houseBeans[i].getStreetAddress());
			sb.append("\t");
			sb.append(houseBeans[i].getCity());
			sb.append("\t");
			sb.append(houseBeans[i].getState());
			sb.append("\t");
			sb.append(houseBeans[i].getZipCode());
			sb.append("\t");
			//price
			int price = houseBeans[i].getPrice()/1000;
			sb.append("$"); sb.append(price);
			sb.append(" - ");
			sb.append(houseBeans[i].getStreetAddress());
			sb.append("\t");
			//url
			sb.append("http://www.utahrealestate.com/");
			sb.append(houseBeans[i].getMls());
			sb.append("\t");
			//color of spot
			sb.append("FF0000\t");
			//image
			sb.append("http://photo.wfrmls.com/280x210/");
			sb.append(houseBeans[i].getMls());
			sb.append(".jpg\t");
			//group
			sb.append(houseBeans[i].getZipCode());
			sb.append("\t");
			//misc
			sb.append(houseBeans[i].getStatus()); sb.append(", ");
			sb.append(houseBeans[i].getSize()); sb.append(" SqrFt, "); 
			sb.append(houseBeans[i].getBedrooms()); sb.append(" BdRms, "); 
			sb.append(houseBeans[i].getBathrooms()); sb.append(" BthRms, ");
			sb.append(houseBeans[i].getFamilyRooms()); sb.append(" FmlyRms, ");
			sb.append(houseBeans[i].getGarage()); sb.append(" Grg, ");
			sb.append(houseBeans[i].getPercentFinishedBasement()); sb.append("% Bsmt, ");
			sb.append(houseBeans[i].getYearBuilt()); sb.append(" YrBlt");
			
			System.out.println(sb);
		}
	}
	
	/**Uses a hash to toss duplicate records, same mls#*/
	public void makeHouseBeans(){
		HashMap map = new HashMap();
		for (int i=0; i< listings.length; i++) {
			HouseBean hb = new HouseBean(listings[i]);
			if (hb.isBroken()){
				System.out.println("\nBroken, not adding "+i);
				Misc.printArray(listings[i]);
			}
			else map.put(new Integer(hb.getMls()), hb);
		}
		//run through map
		Iterator it = map.keySet().iterator();
		houseBeans = new HouseBean[map.size()];
		int counter = 0;
		while (it.hasNext()){
			houseBeans[counter] = (HouseBean) map.get(it.next());
			//System.out.println(houseBeans[counter]+"\n");
			counter++;
		}
	}
	
	/**Parses listing lines on 'List Price'*/
	public void parseListingsFile(){
		ArrayList records = new ArrayList();
		try {
			Pattern priceRE = Pattern.compile(".+List Price.+");
			Matcher mat;
			BufferedReader in = new BufferedReader ( new FileReader(htmlFile));
			String line;
			ArrayList al = new ArrayList();
			//find first price line
			while ((line= in.readLine()) !=null){
				mat = priceRE.matcher(line);
				if (mat.matches()) break;
			}
			if (line == null) Misc.printExit("\nNo price line found in html dump file!\n");
			//make records
			while (true){
				al.add(line);
				while ((line= in.readLine()) !=null){
					mat = priceRE.matcher(line);
					if (mat.matches()) break;
					else al.add(line);
				}
				//add record
				String[] record = new String[al.size()];
				al.toArray(record);
				records.add(record);
				//check if it was last line
				if (line == null) break;
				//nope clear ArrayList and parse next
				al.clear();
			}
			//convert records to String[][]
			listings = new String[records.size()][];
			for (int i=0; i<listings.length; i++){
				listings[i] = (String[]) records.get(i);
			}
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	public static void main (String[] args){
		if (args.length ==0) Misc.printExit("\nEnter the full path file text for an HTML dump from the MLS.\n");
		new ParseMLSHTML(args);
	}
	
	
	
	
	
}
