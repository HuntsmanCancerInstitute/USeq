package edu.utah.id;

import java.util.Arrays;
import java.util.Random;

import org.apache.commons.text.similarity.LevenshteinDistance;

import util.gen.IO;
import util.gen.Misc;

public class MatcherEngine implements Runnable {

	private LevenshteinDistance ld = LevenshteinDistance.getDefaultInstance();
	private boolean failed = false;
	// These subjects are only only present in this thread
	private Subject[] subjectChunk = null;
	private Subject[] testSubjects = null;
	private double missingOneKeyPenalty = 0;
	private double missingAdditionalKeyPenalty = 0;
	private int numMatchesToReturn = 0;
	
	
	
	public MatcherEngine(Subject[] subjectChunk, SubjectIdMatchMaker pm) {
		this.subjectChunk = subjectChunk;
		testSubjects = pm.getTestSubjects();
		missingOneKeyPenalty = pm.getMissingOneKeyPenalty();
		missingAdditionalKeyPenalty = pm.getMissingAdditionalKeyPenalty();
		numMatchesToReturn = pm.getNumberTopMatchesToReturn();
		
	}

	
	public void run() {	
		try {
			//create a random set of index positions
			int[] indexes = new int[testSubjects.length];
			for (int i=0; i< indexes.length; i++) indexes[i] = i;
			Misc.randomize(indexes, new Random());
			
			//for each subject, find the top hits from the threads chunk of db subjects
			for (int i=0; i< indexes.length; i++) {
				Subject test = testSubjects[indexes[i]];
				findTopMatches(test);
			}
			
		} catch (Exception e) {
			failed = true;
			System.err.println("Error: problem matching subjects" );
			e.printStackTrace();
		}
	}


	/*Find top matches*/
	private void findTopMatches(Subject test) {
		
		//set the scores and sort smallest to largest
		String[] testKeys = test.getComparisonKeys();
		for (Subject c: subjectChunk) c.setScore(scoreKeysLD(testKeys, c.getComparisonKeys()));
		Arrays.sort(subjectChunk);
		
		//add top hits to the test subject in a thread safe manner
		Subject[] topHits = new Subject[numMatchesToReturn];
		for (int i=0; i<topHits.length; i++)topHits[i] = subjectChunk[i];
		test.addTopCandidates(topHits);
		
	}


	/**Score keys using Levenshtein Distance
	 * If more than one key is missing, a value of 1 is added to the return score for each.  If just one, then it is ignored.
	 * Thus it's ok to be missing one key, but afterward the penalty is severe. */
	private double scoreKeysLD(String[] test, String[] db) {

//IO.pl("\nT: "+Misc.stringArrayToString(test, ",")+"\nD: "+Misc.stringArrayToString(db, ","));
			//for each key
			double sum = 0;
			int numMissing = 0;
			for (int i=0; i< test.length; i++) {
//IO.p("\t"+test[i]+" vs "+db[i]+" ");
				//missing?
				if (test[i].length() == 0 || db[i].length() == 0) {
					numMissing++;
//IO.pl("missing");
				}
				else {
					double edits = ld.apply(test[i], db[i]);
					double length = test[i].length();
					double ws = edits/length;
//IO.pl(edits+"/"+length+"="+ws);
					sum+= ws;
				}
			}
			//just one missing? then return sum of edits, otherwise add 1 for each
			if (numMissing !=0) {
				if (numMissing == 1) sum+= missingOneKeyPenalty;
				else {
					sum = sum + missingOneKeyPenalty + ((numMissing-1)* missingAdditionalKeyPenalty);
				}
//IO.pl("\tadding missing penalties");
			}
//IO.pl("\t\tScore: "+sum);
			return sum;
	}


	public boolean isFailed() {
		return failed;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	/*
	LevenshteinDistance ld = LevenshteinDistance.getDefaultInstance();
	//						FNLN	Gen		DoB				MRN
	String[] subject = {"DavidNix", "M", "08/11/1968", "12345678"};
	String[][] entries = {
			{"PaolaNix", "F", "06/09/1972", "88534667"},
			{"LarryNix", "M", "08/11/1968", "12345678"},
			{"DavidNix", "M", "09/11/1968", "12345678"},
			{"", "M", "09/11/1968", "12345678"}
	};
	
	for (String[] dbEntry: entries) {
		IO.pl(Misc.stringArrayToString(dbEntry, " "));
		double numEdits = 0;
		double length = 0;
		double indiDist = 0;
		for (int i=0; i< subject.length; i++) {
			//missing?
			if (subject[i].length() == 0 || dbEntry[i].length() == 0) continue;
			double e = ld.apply(subject[i], dbEntry[i]);
			numEdits += e;
			double l = subject[i].length();
			length += l;
			indiDist+= e/l;
			
		}
		double distanceMetric = numEdits/length;
		IO.pl("\t"+Num.formatNumber(distanceMetric, 3)+"\t"+Num.formatNumber(indiDist, 3)+"\t("+numEdits+"/"+length+")");
	}
	*/
	
}
