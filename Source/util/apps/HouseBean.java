package util.apps;
import java.util.regex.*;
import util.gen.*;

public class HouseBean {
	
	//fields
	private int mls = -1;
	private int price = -1;
	private String status = null;
	private String streetAddress = null;
	private String city = null;
	private String state = null;
	private int zipCode = -1;
	private double latitude = -1;
	private double longitude = -1;
	private int size = -1;
	private double bathrooms = -1;
	private int bedrooms = -1;
	private int familyRooms = -1;
	private int garage = -1;
	private int yearBuilt = -1;
	private int percentFinishedBasement = -1;
	private String[] htmlLines;
	private int lineIndex = 0;
	private boolean broken = false;
	
	//regular expressions for parsing fields
	private static final Pattern priceRE = Pattern.compile(".+List Price.+>\\$(.+)\\s*</a.+"); //must remove the , in 450,000
	private static final Pattern mlsRE = Pattern.compile(".+MLS.+B>\\s*(\\d+)\\s*<.+");
	private static final Pattern statusRE = Pattern.compile(".+Stat:.+B>\\s*(.+)\\s*</TD.+");
	private static final Pattern streetAddressRE = Pattern.compile(".+Address.+false'>\\s*(.+)\\s*</a.+"); 
	private static final Pattern cityStateZipRE = Pattern.compile(".+City:</B>\\s*(.+),\\s*&nbsp;\\s*(\\w\\w)\\s+(\\d+)\\s*</TD.+");
	private static final Pattern sizeRE = Pattern.compile(".+Tot Sq Ft:</b>\\s*(\\d+)\\s*<BR.+"); 
	private static final Pattern bedroomsRE = Pattern.compile(".+Tot Beds:</b>\\s*(\\d+)\\s*<BR.+");
	private static final Pattern bathroomsRE = Pattern.compile(".+Tot Baths:</b>\\s*([\\d\\.]+)\\s*<BR.+");
	private static final Pattern familyRoomsRE= Pattern.compile(".+Family Rms:</b>\\s*(\\d+)\\s*<BR.+");
	private static final Pattern garageRE = Pattern.compile(".+Gar\\|Port:</B>\\s*(\\d+)\\s*\\|.+");
	private static final Pattern yearBuiltRE = Pattern.compile(".+Yr Built:</b>\\s*(\\d+)\\s*<BR.+"); 
	private static final Pattern percentFinishedBasementRE = Pattern.compile(".+Fin Bsmnt:</b>\\s*(\\d+)%*\\s*<br.+");
	private Matcher mat;
	
	//constructors
	public HouseBean (String[] htmlLines){
		this.htmlLines = htmlLines;
		replaceQuotations();
		parseFields();
		htmlLines = null;
	}
	
	//methods
	public String toString(){
		StringBuffer sb = new StringBuffer();
		sb.append(mls+"\tMLS#\n");
		sb.append(price+"\tPrice\n");
		sb.append(status+"\tStatus\n");
		sb.append(streetAddress+"\tStreet Address\n");
		sb.append(city+"\tCity\n");
		sb.append(state+"\tState\n");
		sb.append(zipCode+"\tZip Code\n");
		sb.append(latitude+"\tLatitude\n");
		sb.append(longitude+"\tLongitude\n");
		sb.append(size+"\tSize\n");
		sb.append(bathrooms+"\tBathrooms\n");
		sb.append(bedrooms+"\tBedrooms\n");
		sb.append(familyRooms+"\tFamily Rooms\n");
		sb.append(garage+"\tGarage\n");
		sb.append(yearBuilt+"\tYear Built\n");
		sb.append(percentFinishedBasement+"\tPercent Finished Basement\n");
		return sb.toString();
	}
	public void replaceQuotations(){
		//replace any " with '
		Pattern pat = Pattern.compile("\"");
		for (int i=0; i< htmlLines.length; i++){
			mat = pat.matcher(htmlLines[i]);
			htmlLines[i]= mat.replaceAll("'");
		}
	}
	
	public void parseFields (){
		//price
		if (matchIt(priceRE)) price = Integer.parseInt(mat.group(1).replaceAll(",",""));
		//mls
		if (matchIt(mlsRE)) mls = Integer.parseInt(mat.group(1));
		//status
		if (matchIt(statusRE)) status = mat.group(1);
		//address
		if (matchIt(streetAddressRE)) streetAddress = mat.group(1);
		//city state zip
		if (matchIt(cityStateZipRE)) {
			city = mat.group(1);
			state = mat.group(2);
			zipCode = Integer.parseInt(mat.group(3));
		}
		//size
		if (matchIt(sizeRE)) size = Integer.parseInt(mat.group(1));
		//bedrooms
		if (matchIt(bedroomsRE)) bedrooms = Integer.parseInt(mat.group(1));
		//bathrooms
		if (matchIt(bathroomsRE)) bathrooms = Double.parseDouble(mat.group(1));
		//familyRooms
		if (matchIt(familyRoomsRE)) familyRooms = Integer.parseInt(mat.group(1));
		//garages
		if (matchIt(garageRE)) garage = Integer.parseInt(mat.group(1));
		//year built
		if (matchIt(yearBuiltRE)) yearBuilt = Integer.parseInt(mat.group(1));
		//percent finished basement
		if (matchIt(percentFinishedBasementRE)) percentFinishedBasement = Integer.parseInt(mat.group(1));
	}
	
	public boolean matchIt(Pattern pat){
		int oldIndex = lineIndex+1;
		//attempt match starting from prior line index
		for (; lineIndex< htmlLines.length; lineIndex++){
			mat = pat.matcher(htmlLines[lineIndex]);
			if (mat.matches()) return true;
		}
		//no match look at beginning lines
		System.out.println("\tNo Match, looking at previous lines...");
		for (int i=0; i< oldIndex; i++){
			mat = pat.matcher(htmlLines[i]);
			if (mat.matches()) return true;
		}
		System.out.println("\tNo Match! Skipping \n\t"+pat);
		//reset lineIndex to original line
		lineIndex = oldIndex -1;
		broken = true;
		return false;
	}

	
	public static void main (String[] args){
		
	String x = "<TD colspan='2'><B>City:</B> Salt Lake City,  &nbsp; UT 84109</TD>";	
	Pattern pat = Pattern.compile(".+City:</B>\\s*(.+),\\s*&nbsp;\\s*(\\w\\w)\\s+(\\d+)\\s*</TD.+");
	Matcher mat = pat.matcher(x);
	if (mat.matches()) {
		System.out.println(mat.group(1));
		System.out.println(mat.group(2));
		System.out.println(mat.group(3));
	}
	else System.out.println("No match");
	}

	public double getBathrooms() {
		return bathrooms;
	}

	public void setBathrooms(double bathrooms) {
		this.bathrooms = bathrooms;
	}

	public int getBedrooms() {
		return bedrooms;
	}

	public void setBedrooms(int bedrooms) {
		this.bedrooms = bedrooms;
	}

	public String getCity() {
		return city;
	}

	public void setCity(String city) {
		this.city = city;
	}

	public int getFamilyRooms() {
		return familyRooms;
	}

	public void setFamilyRooms(int familyRooms) {
		this.familyRooms = familyRooms;
	}

	public int getGarage() {
		return garage;
	}

	public void setGarage(int garage) {
		this.garage = garage;
	}

	public double getLatitude() {
		return latitude;
	}

	public void setLatitude(double latitude) {
		this.latitude = latitude;
	}

	public double getLongitude() {
		return longitude;
	}

	public void setLongitude(double longitude) {
		this.longitude = longitude;
	}

	public int getMls() {
		return mls;
	}

	public void setMls(int mls) {
		this.mls = mls;
	}

	public int getPercentFinishedBasement() {
		return percentFinishedBasement;
	}

	public void setPercentFinishedBasement(int percentFinishedBasement) {
		this.percentFinishedBasement = percentFinishedBasement;
	}

	public int getPrice() {
		return price;
	}

	public void setPrice(int price) {
		this.price = price;
	}

	public int getSize() {
		return size;
	}

	public void setSize(int size) {
		this.size = size;
	}

	public String getState() {
		return state;
	}

	public void setState(String state) {
		this.state = state;
	}

	public String getStatus() {
		return status;
	}

	public void setStatus(String status) {
		this.status = status;
	}

	public String getStreetAddress() {
		return streetAddress;
	}

	public void setStreetAddress(String streetAddress) {
		this.streetAddress = streetAddress;
	}

	public int getYearBuilt() {
		return yearBuilt;
	}

	public void setYearBuilt(int yearBuilt) {
		this.yearBuilt = yearBuilt;
	}

	public int getZipCode() {
		return zipCode;
	}

	public void setZipCode(int zipCode) {
		this.zipCode = zipCode;
	}

	public boolean isBroken() {
		return broken;
	}
	
	
}
