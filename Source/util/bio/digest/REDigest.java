package util.bio.digest;
import java.util.regex.*;
import java.util.*;

import util.gen.*;

/**
 StandAlone class to perform a restriction digestion on a sequence, a hack, uses as a template.
 */

public class REDigest {
	
	
	public static void main (String[] args){
		
		String gtype2c = "/usr/local/jakarta-tomcat-5.0.14/webapps/ROOT/WEB-INF/classes/bioroot/NEBRebasegtype2c.txt";
		NEBParser ez = new NEBParser();
		Enzyme[] ezs = ez.makeEnzymeArray(gtype2c, 6);
		System.out.println(ez.getHeader());
		
		String DNA = "gggcccgggcccgatcgatcatgcatatggatcc";
		ArrayList cutsites = restrictionMapHTML (DNA, ezs);
		
		for (int i =0; i < cutsites.size() ; i++) {
			System.out.print(cutsites.get(i)+"; ");
		}
	}
	
	/**Returns an arraylist of enzymes found to cut the sequence*/
	public static ArrayList restrictionMapHTML (String DNA, Enzyme[] e) {
		ArrayList al  = new ArrayList();
		for (int i = 0; i < e.length; i++) {
			Pattern p = Pattern.compile(e[i].getRegEx(), Pattern.CASE_INSENSITIVE);
			Matcher m = p.matcher(DNA);
			
			if (m.find()){
				int cut = e[i].getCutTop() + m.start();
				StringBuffer sb = new StringBuffer(e[i].getName() + "<small>("+e[i].getRecogSeq()+")</small>: " + cut);
				
				while (m.find()) {
					cut = e[i].getCutTop() + m.start();
					sb.append(", " + cut);
				}    
				al.add(sb.toString()); 
			}
		}
		return al;
	}
	public static String getRestrictionMapString(String seq, Enzyme[] enzymes){
		ArrayList cutSites = restrictionMapHTML (seq, enzymes);
		int len = cutSites.size();
		StringBuffer sb = new StringBuffer();
		for (int i =0; i < len ; i++) {
			sb.append(cutSites.get(i)+"; ");
		}
		return sb.toString();
	}
	/**Runs thru an enzyme array attempting to cut a particular sequence*/
	public static Enzyme[] restrictionMap(String DNA, Enzyme[] e) {
		ArrayList al  = new ArrayList();
		ArrayList cutSites;
		for (int i = 0; i < e.length; i++) {
			Pattern p = Pattern.compile(e[i].getRegEx(), Pattern.CASE_INSENSITIVE);
			Matcher m = p.matcher(DNA);
			
			if (m.find()){
				cutSites = new ArrayList();
				int cut = e[i].getCutTop() + m.start();
				cutSites.add(new Integer(cut));
				while (m.find()) {
					cut = e[i].getCutTop() + m.start();
					cutSites.add(new Integer(cut));
				}    
				e[i].setTopCuts(Misc.integerArrayListToIntArray(cutSites));	 
			}
			else e[i].setTopCuts(null);
			al.add(e[i]);
		}
		Enzyme[] cutters = new Enzyme[al.size()];
		al.toArray(cutters);
		return cutters;
	}
	

	/**Restriction map a sequence. Makes an ordered list of enzymes. First by number of cuts, second alphabetically.*/
	public static String restrictionMapList(String DNA, Enzyme[] e){
		StringBuffer sb = new StringBuffer("");
		StringBuffer zeros = new StringBuffer();
		Enzyme[] cut = restrictionMap(DNA, e);
		Arrays.sort(cut);
		int num = cut.length;
		for (int i=0; i<num; i++){
			//zero?
			if (cut[i].getNumberOfCuts()==0){
				zeros.append(cut[i].getName());
				zeros.append(" (");
				zeros.append(cut[i].getRecogSeq());
				zeros.append(")\n");
			}
			else {
				sb.append(cut[i].getName());
				sb.append(" (");
				sb.append(cut[i].getRecogSeq());
				sb.append(") ");
				sb.append(Misc.intArrayToString(cut[i].getTopCuts(),", "));
				sb.append("\n");
			}
		}
		sb.append("\n**Non Cutters**\n");
		sb.append(zeros);
		return sb.toString();
	}
}

