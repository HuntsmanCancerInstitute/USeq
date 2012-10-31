package util.bio.parsers;
import java.io.*;
import util.gen.*;
import java.util.*;
import java.util.regex.*;

/**Parser for pulling gene annotations into a UCSC table format for Sodalis.  Ugly! Better to go with XML!*/
public class RefSeqAnnotationParser {

	//fields 
	private File[] files = IO.extractFiles(new File ("/Users/davidnix/HCI/PIs/Mulvey/Travis/SeqAnno/GenBankAnno/"), "gb");
	private String chromosome;
	private Pattern divider = Pattern.compile("\\s{5}\\w+\\s+(.+)");
	private int index = 1;
	
	public RefSeqAnnotationParser(){
		for (int i=0; i< files.length; i++){
			chromosome = Misc.removeExtension(files[i].getName());
			parseIt(files[i]);
		}
	}
	
	public void printUCSCTableFormat(LinkedHashMap hash){
		Pattern coor = Pattern.compile("[^\\d]*(\\d+)\\.{2}(\\d+).*");
		Iterator it = hash.keySet().iterator();
		while (it.hasNext()){
			String key = (String) it.next();
			LinkedHashSet attributes = (LinkedHashSet) hash.get(key);
//System.out.println("XXXXXXX key "+key+"\n\t"+attributes);
			//pull coordinates
			Matcher mat = coor.matcher(key);
			if (mat.matches() == false) Misc.printExit("Cannot parse coordinates "+key+" "+attributes);
			int start = Integer.parseInt(mat.group(1)) -1;
			int end = Integer.parseInt(mat.group(2));
			//pull strand
			String strand = "+";
			if (key.startsWith("comp")) strand = "-";
			//find a locus tag display text
			String locusTag = find(attributes,"/locus_tag");
			if (locusTag == null) {
				locusTag = find(attributes,"/db_xref=\"GI:");
				if (locusTag == null) {
					locusTag = find(attributes,"/db_xref=\"GeneID:");
					if (locusTag == null) {
						locusTag = find(attributes,"/note=");
						if (locusTag == null){
//System.err.println("XXXXXXXXXXXXXXXXX Cannot assign locus tag "+key+" "+attributes);
						continue;
						}
						else {
							locusTag = locusTag.replaceAll(" ", "_")+"_"+index;
							index++;
						}
					}
				}
			}
			//Split it and polish
			locusTag = splitAndClean(locusTag);
			//make a real text
			/*
			ArrayList realNameAL = new ArrayList();
			String real = find(attributes,"/db_xref=\"GI:");
			if (real !=null) realNameAL.add( splitAndClean(real));
			real = find(attributes,"/db_xref=\"GeneID:");
			if (real !=null) realNameAL.add( splitAndClean(real));
			real = find(attributes,"/gene=");
			if (real !=null) realNameAL.add( splitAndClean(real));
			real = find(attributes,"/product=");
			if (real !=null) realNameAL.add( splitAndClean(real));
			real = find(attributes,"/note=");
			if (real !=null) realNameAL.add( splitAndClean(real));
			real = find(attributes,"/bound_moiety=");
			if (real !=null) realNameAL.add( splitAndClean(real));
			String realName = Misc.stringArrayListToString(realNameAL, "; ");
			*/
			String realName = null;
			realName = find(attributes,"/gene=");
			if (realName !=null) realName= splitAndClean(realName);
			else realName = locusTag;
			
			// print it
			//displayName realName chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds	
			System.out.println(locusTag+"\t"+realName+"\t"+chromosome+"\t"+strand+"\t"+start+"\t"+end+"\t"+start+"\t"+end+"\t"+1+"\t"+start+"\t"+end);
			
			
		}
	}
	
	public String splitAndClean(String x){
		String[] items = x.split("=");
		return items[1].replaceAll("\"", "");
	}
	
	public String find(HashSet att, String prefix){
		Iterator it = att.iterator();
		while (it.hasNext()){
			String item = (String)it.next();
			if (item.startsWith(prefix)) {
				//advance until hitting a /
//System.out.println("ITEM "+item);
				while (it.hasNext()){
					String next = (String)it.next();
//System.out.println("\tNEXT "+next);
					if (next.trim().startsWith("/")) {
//System.out.println("\t\tBreak");
						break;
					}
					else {
						item = item + " "+next;
//System.out.println("\t\tConcat");
						
					}
				}
//System.out.println("\tReturning "+item);
				return item;
			}
		}
		return null;
	}
	
	
	public boolean skipToFeatures(BufferedReader in){
		try {
			String line;
			//walk through header to first version then features
			while ((line = in.readLine()) !=null){
				if (line.startsWith("VERSION")) {
					String[] tags = line.split("\\s+");
					chromosome = tags[1];
					while ((line = in.readLine()) !=null){
						if (line.startsWith("FEATURES"))return true;
					}
				}
			}
		} catch (Exception e){
			e.printStackTrace();
		}
		return false;
	}
		
	public void parseIt(File file){
		try {
			BufferedReader in = new BufferedReader (new FileReader (file));
			
			//walk through header to FEATURES
			if (skipToFeatures(in) == false) Misc.printExit("FEATURES line not found.");
			
			//collect lines
			String line;
			LinkedHashMap features = new LinkedHashMap();
			LinkedHashSet attributes = null;
			while ((line = in.readLine()) !=null){ 
				if (line.startsWith("ORIGIN")) {
					printUCSCTableFormat(features);
					if (skipToFeatures(in) == false) System.exit(0);
				}
				Matcher mat = divider.matcher(line);
				if (mat.matches()) {
					//does it exist?
					if (features.containsKey(mat.group(1))) attributes = (LinkedHashSet) features.get(mat.group(1));
					else {
						attributes = new LinkedHashSet();
						features.put(mat.group(1), attributes);
					}
				}
				else {
					attributes.add(line.trim());
				}
			}

			in.close();
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
		
		
		
	public static void main(String[] args) {
		new RefSeqAnnotationParser();

	}

}
