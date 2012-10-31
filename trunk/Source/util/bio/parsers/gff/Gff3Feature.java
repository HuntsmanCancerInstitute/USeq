package util.bio.parsers.gff;

import java.util.*;
import java.io.*;

import util.gen.*;

/**
 A container for each of the values in a GFF3 line. See http://flybase.net/annot/gff3.html.
 Does not look for multiple values in the attributes tag=value1,value2,value3 keeps value as a single String.
 */
public class Gff3Feature implements Comparable, Cloneable{
	
	//fields 
	private String seqId;     
	private String source;     
	private String type;
	private int start;  
	private int end;			
	private double score;          
	private String strand;				
	private int phase;	
	
	//number flags, needed to provide a "." if not set when printing a GFF3 line
	private boolean startSet = false;
	private boolean endSet = false;
	private boolean scoreSet = false;
	private boolean phaseSet = false;
	
	//attributes
	private String attributes;
	private String[] customAttributes;
	private String id;
	private String name;
	private String alias;
	private String parent;
	private String target;
	private String gap;
	private String note;
	private String dbxref;
	private String ontologyTerm;
	
	private boolean valid = true;
	private String sortBy="";
	
	/**Sorts by the sortBy field so set it before sorting.*/
	public int compareTo(Object obj){
		Gff3Feature other = (Gff3Feature) obj;
		return other.sortBy.compareTo(this.sortBy) * -1;
	}
	
	public Object clone(){
		try{
			return super.clone();
		}catch (CloneNotSupportedException e){
			e.printStackTrace();
			return null;
		}
	}
	
	/**Generates a GFF3 line.*/
	public String toString(){
		StringBuffer sb = new StringBuffer();
		//columns
		loadField(sb, seqId);
		sb.append("\t");
		loadField(sb, source);
		sb.append("\t");
		loadField(sb, type);
		sb.append("\t");
		if (startSet) sb.append(start);
		else sb.append(".");
		sb.append("\t");
		if (endSet) sb.append(end);
		else sb.append(".");
		sb.append("\t");
		if (scoreSet) sb.append(score);
		else sb.append(".");
		sb.append("\t");
		loadField(sb, strand);
		sb.append("\t");
		if (phaseSet) sb.append(phase);
		else sb.append(".");
		sb.append("\t");
		//attributes
		loadAttribute(sb, "ID",id);
		loadAttribute(sb, "Name",name);
		loadAttribute(sb, "Alias",alias);
		loadAttribute(sb, "Parent",parent);
		loadAttribute(sb, "CompositCNV",target);
		loadAttribute(sb, "Gap",gap);
		loadAttribute(sb, "Note",note);
		loadAttribute(sb, "Dbxref",dbxref);
		loadAttribute(sb, "Ontology_term",ontologyTerm);
		//custom attributes 
		int numCustom = 0;
		if (customAttributes!=null) numCustom= customAttributes.length;
		if (numCustom !=0 ){
			sb.append(customAttributes[0]);
			for (int i=1; i<numCustom; i++){
				sb.append(";");
				sb.append(customAttributes[i]);
			}
		}
		//nip final ; off the stop if no custom attributes
		String line = sb.toString();
		if (line.endsWith(";")) return line.substring(0,line.length()-1);
		return sb.toString();
	}

	/**Generates a GFF3 line.*/
	public String toStringNoAttributes(){
		StringBuffer sb = new StringBuffer();
		//columns
		loadField(sb, seqId);
		sb.append("\t");
		loadField(sb, source);
		sb.append("\t");
		loadField(sb, type);
		sb.append("\t");
		if (startSet) sb.append(start);
		else sb.append(".");
		sb.append("\t");
		if (endSet) sb.append(end);
		else sb.append(".");
		sb.append("\t");
		if (scoreSet) sb.append(score);
		else sb.append(".");
		sb.append("\t");
		loadField(sb, strand);
		sb.append("\t");
		if (phaseSet) sb.append(phase);
		else sb.append(".");
		return sb.toString();
	}

	/**Returns a collection of sgr lines representing the gff feature:
	 * 1) seqId startMin1 0
	 * 2) seqId start score
	 * 3) seqId stop score
	 * 4) seqId endPlus1 0
	 * Does not include a final return, thus println(gff.sgrBundle()).*/
	public String sgrBundle(){
		StringBuffer sb = new StringBuffer();
		int startMinOne = start -1;
		if (startMinOne < 0) startMinOne = 0;
		sb.append(seqId); sb.append("\t"); sb.append(startMinOne); sb.append("\t0\n");
		sb.append(seqId); sb.append("\t"); sb.append(start); sb.append("\t"); sb.append(score); sb.append("\n");
		sb.append(seqId); sb.append("\t"); sb.append(end); sb.append("\t"); sb.append(score); sb.append("\n");
		sb.append(seqId); sb.append("\t"); sb.append(end+1); sb.append("\t0");
		return sb.toString();
	}
	
	public static void loadAttribute(StringBuffer sb, String reserveWord, String value){
		if (Misc.isNotEmpty(value)){
			sb.append(reserveWord);
			sb.append("=");
			sb.append(value);
			sb.append(";");
		}
	}
	
	public static void loadField(StringBuffer sb, String field){
		if (Misc.isEmpty(field)) sb.append(".");
		else sb.append(field);
	}
	
	//constructors
	public Gff3Feature(){}
	
	/**Always check if the GFF3Feature is valid after instantiating with this constructor.
	 * Invalid cases arise from not having 9 tab delimited columns or having an attribute 
	 * that does not split in two on '='. Escaped characters are OK. Will also fail if 
	 * numbers could not be parsed correctly.*/
	public Gff3Feature(String unParsedGff3Line){
		String[] items = unParsedGff3Line.trim().split("\\t");
		
		//check length
		if (items.length>7){
			try{
				//set columns
				seqId = items[0];
				source = items[1];
				type = items[2];
				if (items[3].equals(".") == false) {
					start = Integer.parseInt(items[3]); 
					startSet = true;
				}
				if (items[4].equals(".") == false) {
					end = Integer.parseInt(items[4]); 
					endSet = true;
				}
				if (items[5].equals(".") == false) {
					score = Double.parseDouble(items[5]); 
					scoreSet = true;
				}
				strand = items[6];
				if (items[7].equals(".") == false) {
					phase = Integer.parseInt(items[7]); 
					phaseSet = true; 
				}
			} catch (NumberFormatException e){
				valid = false;
			}
			//any attributes?
			if (items.length>8){
			//split attributes on ; but not /;
			attributes = items[8];
			String[] attributes = Misc.splitString(items[8], ";");
			int numAtts = attributes.length;
			ArrayList custom = new ArrayList();
			String[] tagValue;
			for (int i=0; i<numAtts; i++){
				//split tag on = but not /=
				tagValue = Misc.splitString(attributes[i],"=");
				//check if ok if not break
				if (tagValue.length!=2) {
					valid = false;
					break;
				}
				//attempt to assign reserved tags, don't use a hash since there may be duplicates
				if (tagValue[0].equals("ID")) id = tagValue[1];
				else if (tagValue[0].equals("Name")) name = tagValue[1];
				else if (tagValue[0].equals("Alias")) alias = tagValue[1];
				else if (tagValue[0].equals("Parent")) parent = tagValue[1];
				else if (tagValue[0].equals("CompositCNV")) target = tagValue[1];
				else if (tagValue[0].equals("Gap")) gap = tagValue[1];
				else if (tagValue[0].equals("Note")) note = tagValue[1];
				else if (tagValue[0].equals("Dbxref")) dbxref = tagValue[1];
				else if (tagValue[0].equals("Ontology_term")) ontologyTerm = tagValue[1];
				else {
					custom.add(attributes[i]);
				}
			}
			//convert custom attributes to a String[]
			customAttributes = new String[custom.size()];
			custom.toArray(customAttributes);
			}
			
		}
		//incorrect number of columns
		else valid = false;
	}
	
	/**Prints an array of GFF3Features as a plain text file*/
	public static void saveGFF3TextFile(Gff3Feature[] g, File newGFF3File){
		int num = g.length;
		StringBuffer sb = new StringBuffer();
		for (int i=0; i<num; i++){
			sb.append(g[i].toString());
			sb.append("\n");
		}
		IO.writeString(sb.toString(), newGFF3File);
	}
	
	/**Reserved Attribute.*/
	public String getAlias() {
		return alias;
	}
	public void setAlias(String alias) {
		this.alias = alias;
	}
	/**Any Attribute that is not a reserve attribute is grouped here. 
	 * Ie {"cat= dog", "tag=value1, value2, value3"}.*/
	public String[] getCustomAttributes() {
		return customAttributes;
	}
	public void setCustomAttributes(String[] customAttributes) {
		this.customAttributes = customAttributes;
	}
	/**Reserved Attribute.*/
	public String getDbxref() {
		return dbxref;
	}
	public void setDbxref(String dbxref) {
		this.dbxref = dbxref;
	}
	public int getEnd() {
		return end;
	}
	public void setEnd(int end) {
		this.end = end;
		endSet=true;
	}
	/**Reserved Attribute.*/
	public String getGap() {
		return gap;
	}
	public void setGap(String gap) {
		this.gap = gap;
	}
	/**Reserved Attribute.*/
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	/**Reserved Attribute.*/
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	/**Reserved Attribute.*/
	public String getNote() {
		return note;
	}
	public void setNote(String note) {
		this.note = note;
	}
	/**Reserved Attribute.*/
	public String getOntologyTerm() {
		return ontologyTerm;
	}
	public void setOntologyTerm(String ontologyTerm) {
		this.ontologyTerm = ontologyTerm;
	}
	/**Reserved Attribute.*/
	public String getParent() {
		return parent;
	}
	public void setParent(String parent) {
		this.parent = parent;
	}
	public int getPhase() {
		return phase;
	}
	public void setPhase(int phase) {
		this.phase = phase;
		phaseSet = true;
	}
	public double getScore() {
		return score;
	}
	public void setScore(double score) {
		this.score = score;
		scoreSet = true;
	}
	public String getSeqId() {
		return seqId;
	}
	public void setSeqId(String seqId) {
		this.seqId = seqId;
	}
	public String getSource() {
		return source;
	}
	public void setSource(String source) {
		this.source = source;
	}
	public int getStart() {
		return start;
	}
	public void setStart(int start) {
		this.start = start;
		startSet=true;
	}
	public String getStrand() {
		return strand;
	}
	public void setStrand(String strand) {
		this.strand = strand;
	}
	/**Reserved Attribute.*/
	public String getTarget() {
		return target;
	}
	/**Reserved Attribute.*/
	public void setTarget(String target) {
		this.target = target;
	}
	public String getType() {
		return type;
	}
	public void setType(String type) {
		this.type = type;
	}
	/**Call to find out if the unparsed GFF line correctly parsed.*/
	public boolean isValid() {
		return valid;
	}
	public void setValid(boolean valid) {
		this.valid = valid;
	}
	public boolean isEndSet() {
		return endSet;
	}
	public void setEndSet(boolean endSet) {
		this.endSet = endSet;
	}
	public boolean isPhaseSet() {
		return phaseSet;
	}
	public void setPhaseSet(boolean phaseSet) {
		this.phaseSet = phaseSet;
	}
	public boolean isScoreSet() {
		return scoreSet;
	}
	public void setScoreSet(boolean scoreSet) {
		this.scoreSet = scoreSet;
	}
	public boolean isStartSet() {
		return startSet;
	}
	public void setStartSet(boolean startSet) {
		this.startSet = startSet;
	}
	public String getSortBy() {
		return sortBy;
	}
	public void setSortBy(String sortBy) {
		this.sortBy = sortBy;
	}
	public String getAttributes() {
		return attributes;
	}
}
