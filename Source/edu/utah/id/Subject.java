package edu.utah.id;

import java.io.IOException;
import util.gen.Misc;

public class Subject implements Comparable<Subject> {

	private String lastName = "";
	private String firstName = "";
	private int dobMonth = -1;
	private int dobDay = -1;
	private int dobYear = -1;
	private String gender = "";
	private String mrn = "";
	private String[] hciSubjectIds = null;
	private String[] coreIds = null;

	private double score = 0;
	private String[] comparisonKeys = null;
	private Subject[] topMatches = null;
	private double[] topMatchScores = null;

	//constructor
	public Subject(String dataLine, boolean addCoreId, CoreId coreId) throws IOException {

		//required: lastName firstName dobMonth dobDay dobYear gender mrn 
		//             0        1          2       3      4      5     6
		//optional: coreIds hciPersonIds
		//              7        8
		String[] t = Misc.TAB.split(dataLine);
		if (t.length< 7) throw new IOException("ERROR: too few fields in subject dataline : "+dataLine);
		//remove leading or trailing whitespace
		for (int i=0; i<t.length; i++) t[i] = t[i].trim();

		lastName = t[0];
		firstName = t[1];

		if (t[2].length()!=0) {
			dobMonth = Integer.parseInt(t[2]);
			if (dobMonth< 1 || dobMonth > 12) throw new IOException("ERROR: dob month field '"+t[2]+"' is malformed, must be 1-12, in subject dataline : "+dataLine);
		}

		if (t[3].length()!=0) {
			dobDay = Integer.parseInt(t[3]);
			if (dobDay< 1 || dobDay > 31) throw new IOException("ERROR: dob day field '"+t[3]+"' is malformed, must be 1-31, in subject dataline : "+dataLine);
		}

		if (t[4].length()!=0) {
			dobYear = Integer.parseInt(t[4]);
			if (dobYear< 1900 || dobYear > 2050) throw new IOException("ERROR: dob year field '"+t[4]+"' is malformed, must be 1900-2050, in subject dataline : "+dataLine);
		}

		if (t[5].length()!=0) {
			gender = t[5];
			if (gender.equals("M") == false && gender.equals("F") == false) throw new IOException("ERROR: the gender field '"+t[5]+"' is malformed, must be M or F, in subject dataline : "+dataLine);
		}

		if (t[6].length()!=0) mrn= t[6];

		if (t.length > 7) {
			//LLDLLDLL but no 0 or O, all upper case
			coreIds = Misc.COMMA.split(t[7]);
			for (String id: coreIds) if (CoreId.isCoreId(id)==false) throw new IOException("ERROR: the coreId '"+id+"' is malformed, (must be "+CoreId.CORE_ID_DESCRIPTION+") in subject dataline : "+dataLine);
		}
		else if (addCoreId) coreIds = new String[] {coreId.createCoreId()};

		if (t.length > 8) {
			hciSubjectIds = Misc.COMMA.split(t[8]);
		}

		makeComparisonKeys();

	}

	/**Leave missing data as "", these will be skipped.*/
	private void makeComparisonKeys() {
		String dob = "";
		if (dobMonth!=-1 && dobDay!=-1 && dobYear!=-1) dob = dobMonth+"/"+dobDay+"/"+dobYear;
		comparisonKeys = new String[] {
				lastName+ firstName,
				dob,
				gender,
				mrn
		};
	}

	public synchronized void addTopCandidates(Subject[] topHits) {
		// yet instantiated?
		if (topMatches == null) {
			topMatches = topHits;
			topMatchScores = new double[topHits.length];
			for (int i=0; i< topHits.length; i++) topMatchScores[i] = topHits[i].getScore();
		}
		else {
			//for each topHit from the chunk, compare to what this subject has already seen
			for (int i=0; i< topHits.length; i++) {
				if (topHits[i].getScore() < topMatchScores[i]) {
					topMatches[i] = topHits[i];
					topMatchScores[i] = topHits[i].getScore();
				}
			}
		}
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(lastName); sb.append("\t");
		sb.append(firstName); sb.append("\t");
		if (dobMonth != -1) sb.append(dobMonth);
		sb.append("\t");
		if (dobDay != -1) sb.append(dobDay);
		sb.append("\t");
		if (dobYear != -1) sb.append(dobYear);
		sb.append("\t");
		sb.append(gender); sb.append("\t");
		sb.append(mrn); sb.append("\t");
		if (coreIds !=null) {
			sb.append(coreIds[0]);
			for (int i=1; i< coreIds.length; i++) {
				sb.append(",");
				sb.append(coreIds[i]);
			}
		}
		sb.append("\t");
		if (hciSubjectIds !=null) {
			sb.append(hciSubjectIds[0]);
			for (int i=1; i< hciSubjectIds.length; i++) {
				sb.append(",");
				sb.append(hciSubjectIds[i]);
			}
		}
		return sb.toString();
	}
	
	public String toStringPretty() {
		StringBuilder sb = new StringBuilder(comparisonKeys[0]);
		for (int i=1; i< comparisonKeys.length; i++) {
			sb.append("\t");
			sb.append(comparisonKeys[i]);
		}
		if (coreIds == null && hciSubjectIds==null) return sb.toString();
		
		sb.append("\t");
		if (coreIds !=null) {
			sb.append(coreIds[0]);
			for (int i=1; i< coreIds.length; i++) {
				sb.append(",");
				sb.append(coreIds[i]);
			}
		}
		sb.append("\t");
		if (hciSubjectIds !=null) {
			sb.append(hciSubjectIds[0]);
			for (int i=1; i< hciSubjectIds.length; i++) {
				sb.append(",");
				sb.append(hciSubjectIds[i]);
			}
		}
		return sb.toString();
		
	}

	/**Sorts by score, smallest to largest*/
	public int compareTo(Subject o) {
		if (o.score < this.score) return 1;
		if (o.score > this.score) return -1;
		return 0;
	}

	public String[] getComparisonKeys() {
		return comparisonKeys;
	}
	public double getScore() {
		return score;
	}
	public void setScore(double score) {
		this.score = score;
	}
	public Subject[] getTopMatches() {
		return topMatches;
	}
	public double[] getTopMatchScores() {
		return topMatchScores;
	}

	public String[] getCoreIds() {
		return coreIds;
	}








}
