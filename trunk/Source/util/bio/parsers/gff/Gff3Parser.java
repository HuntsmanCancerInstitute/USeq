package util.bio.parsers.gff;
import java.io.*;
import java.util.*;
import java.util.regex.*;

/**
 * For parsing a GFF3 text file into GFF3Feature objects, be sure to set the regExTypes String filter to match which types you want to save!
 */
public class Gff3Parser {
	//fields
	private File GFFFile;
	/**Used to describe what types of gff lines are to be made into GFF3Feature objects.*/
	private String regExTypes = "gene|ncRNA|snoRNA|tRNA|rRNA|transposable_element|pseudogene|CDS|mRNA|exon|CRM|snRNA|region";
	private Gff3Feature[] features;
	private Gff3Feature[][] chromSplitFeatures;
	private boolean printNotes = true;
	private boolean relax = false;	//relaxes attribute type requirements
	
	//constructor
	public Gff3Parser(File GFF3File){
		parseIt(GFF3File);
	}
	public Gff3Parser(){
	}
	
	public boolean parseIt(File GFF3File){
		GFFFile = GFF3File;
		Pattern p = Pattern.compile(regExTypes);
		Matcher m;
		
		try{
			BufferedReader in = new BufferedReader(new FileReader(GFFFile));
			//PrintWriter out = new PrintWriter(new FileWriter(new File(GFFFile.getCanonicalPath()+"Str")));
			String line;
			Gff3Feature feature;
			int counter = 0;
			if (printNotes) System.out.println ("\nParsing "+GFFFile.getName()+"...\n");
			HashSet types = new HashSet();
			HashSet disgardedTypes = new HashSet();
			ArrayList goodFeatures = new ArrayList();
			
			while ((line = in.readLine())!=null){
				counter++;
				//limit line reading?
				//if (counter == 25) break;
				line = line.trim();
				//print comments and spaces
				if (line.startsWith("#")) System.out.println("Comment: "+line);
				else if (line.length()==0) {}
				else {
					feature = new Gff3Feature(line);
					//is feature good and has the correct type				
					m = p.matcher(feature.getType());					
					if (feature.isValid() || relax){
						if (m.matches()) {
							goodFeatures.add(feature);
							//out.println(feature);
							types.add(feature.getType());
						}
						else disgardedTypes.add(feature.getType());
					}
					else {
						System.out.println("ERROR! Line "+counter+" is malformed -> "+line);
						break;
					}
				}
			}
			in.close();
			//out.close();
			
			if (printNotes) System.out.println("\nSaved Types: "+types);
			if (printNotes) System.out.println("\nDisgarded Types: "+disgardedTypes);
			
			//convert and sort
			if (printNotes) System.out.println("\nMaking correct type features...\n");
			features = new Gff3Feature[goodFeatures.size()];
			goodFeatures.toArray(features);
			return true;
		}catch (Exception e){
			e.printStackTrace();
			return false;
		}
	}
	
	/**Subract one from start and stop of each feature.  Use for converting to zero base coordinates.*/
	public void subtractOneFromFeatures(){
		for (int i=0; i< features.length; i++){
			features[i].setStart(features[i].getStart()-1);
			features[i].setEnd(features[i].getEnd()-1);
		}
		chromSplitFeatures = null;
	}
	
	/**Returns the features split by chromosome.*/
	public Gff3Feature[][] getChromSplitFeatures(){
		if (chromSplitFeatures != null) return chromSplitFeatures;
		//add to hash
		LinkedHashMap hash = new LinkedHashMap();
		ArrayList al;
		for (int i=0; i< features.length; i++){
			if (hash.containsKey(features[i].getSeqId())){
				al = (ArrayList) hash.get(features[i].getSeqId());
			}
			else {
				al = new ArrayList();
				hash.put(features[i].getSeqId(), al);
			}
			al.add(features[i]);
		}
		//convert to array
		chromSplitFeatures = new Gff3Feature[hash.size()][];
		Iterator it = hash.keySet().iterator();
		int counter = 0;
		while (it.hasNext()){
			al = (ArrayList)hash.get(it.next());
			Gff3Feature[] chr = new Gff3Feature[al.size()];
			al.toArray(chr);
			chromSplitFeatures[counter++] = chr;
		}
		return chromSplitFeatures;
	}
	
	public static void main(String[] args) {
		new Gff3Parser(new File(args[0]));
	}
	public Gff3Feature[] getFeatures() {
		return features;
	}
	public File getGFFFile() {
		return GFFFile;
	}
	public void setGFFFile(File file) {
		GFFFile = file;
	}
	public boolean printNotes() {
		return printNotes;
	}
	public void setPrintNotes(boolean printNotes) {
		this.printNotes = printNotes;
	}
	public String getRegExTypes() {
		return regExTypes;
	}
	/**Use something like -> gene|ncRNA|snoRNA|tRNA|rRNA|transposable_element|pseudogene|CDS|mRNA|exon 
	 * Only those gff lines with a type matching this regular expression will be made into GFF3Feature's.
	 * If you want everything set this to '.+' .
	 * If you don't know what types are present in the gff file, run this parser and it will print all of the types.*/
	public void setRegExTypes(String regExTypes) {
		this.regExTypes = regExTypes;
	}
	public void setRelax(boolean relax) {
		this.relax = relax;
	}
}
