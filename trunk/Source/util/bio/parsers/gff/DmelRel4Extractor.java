package util.bio.parsers.gff;

import java.util.regex.*;
import java.util.*;
import java.io.*;

import util.bio.annotation.*;
import util.gen.*;

/**
 * Class for building gene models from GFF3Features, specifically tailored to release 4.0 drosophila melanogaster GFF.
 * The start is always less than the stop.  
 * The actual length of an element is stop+1-start. 
 * Orientation is 1,-1, or 0 corresponding to +,-, or .
 */
public class DmelRel4Extractor {
	//fields
	//case sensitive descriptors for what type names the gff file uses to call a gene, an exon, a transcript, and a translation
	private String regExGeneTypes = "gene|transposable_element|CRM"; //add new gene type items here!
	private String exonName = "exon";
	private String transcriptName = "mRNA";
	private String translationName = "CDS";
	private ArrayList geneGroups= new ArrayList();
	private String chromosomeNameAppender = null;  //set to something to prepend the text of the chromosome with, ie "chr" since the gff just uses the 3r, 2L
	
	//generic Features
	private LinkedHashSet genericFeaturesHash = new LinkedHashSet(); //used to keep track of the generic features
	private ArrayList genericFeaturesAL = new ArrayList();
	private boolean genericFeaturesFound = false;
	
	//internal use
	private ArrayList exons = new ArrayList();
	private ArrayList translations = new ArrayList();
	private ArrayList transcripts = new ArrayList();
	private ArrayList transGroups = new ArrayList();
	
	
	public static void main(String[] args){
		DmelRel4Extractor x = new DmelRel4Extractor();
		x.extract( new File(args[0]),true);
		ArrayList geneGroups = x.getGeneGroupArrayList();
		//print CG names
		int num = geneGroups.size();
		for (int i=0; i<num; i++){
			String name = ((GeneGroup)geneGroups.get(i)).getName();
			//if (text.startsWith("CG")) System.out.println(text);
			System.out.println(name);
		}
	}
	
	public void extract(File gffFile, boolean stripIt){
		//parse gff file
		System.out.println("\nParsing GFF3 file...");
		Gff3Parser parser = new Gff3Parser();
		parser.setPrintNotes(false);
		parser.setRegExTypes(regExGeneTypes+"|"+exonName+"|"+transcriptName+"|"+translationName);
		parser.parseIt(gffFile);
		Gff3Feature[] features = parser.getFeatures();
		int numFeatures = features.length;
		
		if (stripIt){
			//sort and strip any screwie features
			//all I can say is that this gff file is really inconsistent!  !$#%!#$@!#@!@#$!!!!!
			//sort features by CG|TE|CR number or substitute in parent CG|TE|CR number or use ID
			System.out.println("Checking features for gene models...");
			Pattern patCG = Pattern.compile("CG\\d*");
			Pattern patTE = Pattern.compile("TE\\d*");
			Pattern patCR = Pattern.compile("CR\\d*");
			
			ArrayList filteredFeatures = new ArrayList(numFeatures);
			for (int i=0; i<numFeatures; i++) {
				String id = null;
				//try ID
				id = stripIt(patCG, patTE, patCR, features[i].getId());
				//try Parent
				if (id == null) id = stripIt(patCG, patTE, patCR, features[i].getParent());
				//if still null toss it
				if (id == null) {
					System.out.println("ERROR: Skipping odd ball feature ->"+ features[i]);
				}
				else{
					features[i].setSortBy(id);
					filteredFeatures.add(features[i]);
				}
				//check fix exons with more than one parent, this is a really stupid thing to have to do!
				if (features[i].getType().equals(exonName)){
					String[] parents = features[i].getParent().split(",");
					int counter =1;
					for (int x=0; x<parents.length; x++){
						if (parents[x].indexOf(id) == -1){
							//make new feature
							Gff3Feature bastard = (Gff3Feature)features[i].clone();
							String nakedParent = stripIt(patCG, patTE, patCR, parents[x]);
							bastard.setId(nakedParent+":b"+counter);
							counter++;
							bastard.setSortBy(nakedParent);
							filteredFeatures.add(bastard);
						}
					}
				}
			}
			features = new Gff3Feature[filteredFeatures.size()];
			filteredFeatures.toArray(features);
			System.out.println("Sorting features...");
			Arrays.sort(features);
		}
		//set sortby
		else {
			for (int i=0; i<numFeatures; i++){
				features[i].setSortBy(features[i].getSeqId()+"_"+features[i].getStart()+"_"+features[i].getEnd());
			}
		}
		
		if (features.length == 0){
			System.out.println("\nSorry, no GFF features?!\n");
			System.exit(0);
		}

		//build annotation objects, note everything is compared to lower case names 
		//loop thru all features
		System.out.println("Building gene models...");
		ArrayList group = new ArrayList();
		Pattern geneTypes = Pattern.compile(regExGeneTypes);
		//add first feature
		String strippedId = features[0].getSortBy();
		group.add(features[0]);
		int last = numFeatures-1;
		for (int i = 1; i < numFeatures; i++) {
			//build gene group?  or add feature to arrays?
			if (features[i].getSortBy().equals(strippedId) == false || i==last) {
				//make gene group
				//load exon, transcript, translation arrays and the closer
				int numSub = group.size();
				Gff3Feature closer = null;
				Gff3Feature subFeat = null;
				for (int j=0; j<numSub; j++){
					subFeat = (Gff3Feature)group.get(j);
					String n = subFeat.getType();
					if (n.equals(exonName)) exons.add(subFeat);
					else if (n.equals(translationName)) translations.add(subFeat);	
					else if (n.equals(transcriptName)) transcripts.add(subFeat);
					else if (geneTypes.matcher(n).matches()) closer = subFeat;
					else closer = subFeat;
				}
				//Was there a closing group?
				if (closer == null){
					System.out.println("ERROR: This feature has no closing group! ->"+subFeat+"\n\tAccepted types: "+regExGeneTypes+"\n");
					//clear containers
					exons.clear();
					transcripts.clear();
					translations.clear();
					transGroups.clear();
				}
				//build a gene group or a generic
				else if ( geneTypes.matcher( closer.getType() ) .matches() ){
					buildGeneGrp(closer);
				}
				else {
					buildGenericFeature(closer);
					genericFeaturesHash.add(closer.getType());  //add text to hash to get list
				}
				//reset
				strippedId = features[i].getSortBy();
				group.clear();
				group.add(features[i]);
			}
			//add
			else group.add(features[i]);
		}
		
		//set whether generic features were found
		if (genericFeaturesHash.size()>0) {
			genericFeaturesFound = true;
			System.out.println("Generic Features Found!");
		}
		System.out.println("Extraction complete! "+ (geneGroups.size()+1)+ " Gene Groups");
	}
	
	/*Use to return one of three pattern matches.**/
	public static String stripIt(Pattern patCG, Pattern patTE, Pattern patCR, String toStrip){
		if (toStrip == null) return null;
		Matcher mat = patCG.matcher(toStrip);
		if (mat.find()) return mat.group();
		else {
			mat = patTE.matcher(toStrip);
			if (mat.find()) return mat.group();
			else {
				mat = patCR.matcher(toStrip);
				if (mat.find()) return mat.group();
				else {
					//all failed 
					return null;
				}
			}
		}
	}
	
	/**Contains all the feature types that were not recognized by the Extractor.
	 * Should be user added items like CRMs, enhancers, etc.*/
	public LinkedHashSet getGenericFeatureHash(){
		return genericFeaturesHash;
	}
	public ArrayList getGenericFeatures(){
		return genericFeaturesAL;
	}
	public void buildGenericFeature(Gff3Feature f){
		genericFeaturesAL.add(new GenericFeature(f));
	}
	
	public void buildTransGroups(){
		Gff3Feature f = null;
		ExonIntron[] ex = null;
		Translation translation = null;
		Transcript transcript = null;
		int numTranslations = translations.size();
		int numTranscripts = transcripts.size();
		int numExons = exons.size();
		
		//make ExonIntron[]
		if (numExons !=0) {
			ex = new ExonIntron[numExons];
			for (int i=0; i<numExons; i++){
				f = (Gff3Feature) exons.get(i);
				ex[i]= new ExonIntron(f);
			}
		}
		
		//if only one transgroup
		if (numTranscripts <= 1) {
			//transcripts
			if (numTranscripts !=0) {
				f = (Gff3Feature) transcripts.get(0);
				transcript = new Transcript(f);
			}
			//translation
			translation = null;
			if (numTranslations !=0) {
				f = (Gff3Feature) translations.get(0);
				translation = new Translation(f);
			}
			//trans group
			transGroups.add(new TransGroup( transcript,  translation, ex));
		}
		//more than one trans groups!
		else{
			//convert arrays
			Gff3Feature[] transcriptFeatures = new Gff3Feature[numTranscripts]; 
			Gff3Feature[] translationFeatures = new Gff3Feature[numTranslations];
			transcripts.toArray(transcriptFeatures);
			translations.toArray(translationFeatures);
			
			//run thru each transcript
			for (int i=0; i<numTranscripts; i++){
				//get text of transcript
				String transcriptName = transcriptFeatures[i].getId();
				//run thru translations and find match
				translation = null;
				
				for (int j=0; j<numTranslations; j++) {
					if (translationFeatures[j].getParent().equals(transcriptName)) {
						translation = new Translation(translationFeatures[j]);
						break;
					}
				}
				//run thru exon features
				ArrayList exonGroupAL = new ArrayList();
				for (int k =0; k<numExons; k++){
					if (ex[k].getParent().indexOf(transcriptName) != -1) exonGroupAL.add(ex[k]);
				}
				ExonIntron[] exonGroup = null;
				if (exonGroupAL.size()!=0) {
					exonGroup = new ExonIntron[exonGroupAL.size()];
					exonGroupAL.toArray(exonGroup);
				}
				//make trans group, translation and exons might be null
				transGroups.add(new TransGroup( new Transcript(transcriptFeatures[i]), translation, exonGroup));
			}
		}
	}
	/**Creates a new gene group, call this last after having loaded the exons, transcripts, and translations.*/
	public void buildGeneGrp(Gff3Feature f){
		//build transcription groups
		buildTransGroups();
		//make GeneGroup
		TransGroup[] transGroupArray = new TransGroup[transGroups.size()];
		transGroups.toArray(transGroupArray);
		String chrom;
		if (chromosomeNameAppender!=null) chrom = chromosomeNameAppender+f.getSeqId();
		else chrom = f.getSeqId();
		geneGroups.add(new GeneGroup(f.getId(), chrom,
				f.getStart(),f.getEnd(),
				f.getType(), transGroupArray,
				GeneGroup.convertPlusToNumOrientation(f.getStrand()),
				f.getAttributes()));
		//clear containers
		exons.clear();
		transcripts.clear();
		translations.clear();
		transGroups.clear();
	}
	
	public ArrayList getGeneGroupArrayList(){
		return geneGroups;
	}
	/**Final product from the extractor. Use in downstream analysis.*/
	public GeneGroup[] getGeneGroups(){
		GeneGroup[] gg = new GeneGroup[geneGroups.size()];
		geneGroups.toArray(gg);
		return gg;
	}
	public boolean isGenericFeaturesFound() {
		return genericFeaturesFound;
	}
	public String getRegExGeneTypes() {
		return regExGeneTypes;
	}
	public void setRegExGeneTypes(String regExGeneTypes) {
		this.regExGeneTypes = regExGeneTypes;
	}
	public ArrayList getGenericFeaturesAL() {
		return genericFeaturesAL;
	}
	public LinkedHashSet getGenericFeaturesHash() {
		return genericFeaturesHash;
	}

	public void setChromosomeNameAppender(String chromosomeNameAppender) {
		this.chromosomeNameAppender = chromosomeNameAppender;
	}
}
