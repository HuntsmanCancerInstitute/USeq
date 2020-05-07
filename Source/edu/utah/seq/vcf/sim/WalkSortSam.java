package edu.utah.seq.vcf.sim;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;

import edu.utah.seq.data.sam.PicardSortSam;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import util.gen.IO;
import util.gen.Misc;

public class WalkSortSam implements Runnable {

	private boolean failed = false;
	private File modifiedBam;
	private File unModifiedBam;
	private File finalBam;
	private File filteredBam;
	private int id;
	
	private SamReaderFactory readerFactory = null;
	private SamReader modifiedReader = null;
	private SamReader unModifiedReader = null;
	private SAMRecordIterator modifiedIterator = null;
	private SAMRecordIterator unModifiedIterator = null;
	
	private ArrayList<SAMRecord> modifiedSamRecords = new ArrayList<SAMRecord>();
	private SAMRecord lastModifiedSamRecord = null;
	
	private ArrayList<SAMRecord> unModifiedSamRecords = new ArrayList<SAMRecord>();
	private SAMRecord lastUnModifiedSamRecord = null;
	private SAMFileWriter unsortedWriter = null;
	
	private long numModifiedSaved = 0;
	private long numUnModifiedSaved = 0;
	
	private boolean verbose = true;
	
	public WalkSortSam (File modifiedBam, File unModifiedBam, File finalBam, File filteredBam, int id) {
		this.modifiedBam = modifiedBam;
		this.unModifiedBam = unModifiedBam;
		this.finalBam = finalBam;
		this.filteredBam = filteredBam;
		this.id = id;
	}
	
	public void run() {	
		try {
			IO.pl(id+"\tWalking modified bam, adding unmodified alignments for "+ Misc.removeExtension(modifiedBam.getName()));

			//create readers
			readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			modifiedReader = readerFactory.open(modifiedBam);
			unModifiedReader = readerFactory.open(unModifiedBam);
						
			//fetch iterators and load first record
			modifiedIterator = modifiedReader.iterator();
/*while (modifiedIterator.hasNext()) {
modifiedIterator.next();
numModifiedSaved++;
}*/
			lastModifiedSamRecord = modifiedIterator.next();
			unModifiedIterator = unModifiedReader.iterator();
/*while (unModifiedIterator.hasNext()) {
				unModifiedIterator.next();
				numUnModifiedSaved++;
}*/
			lastUnModifiedSamRecord = unModifiedIterator.next();
//IO.pl(id+ "\tTotals # in ModMix: "+ numModifiedSaved+"  UnMod: "+ numUnModifiedSaved);
//System.exit(0);

			//create a writer for the unsorted output
			File unsortedFinalBam = new File (finalBam.getParentFile(), "unsorted_"+finalBam.getName());
			unsortedFinalBam.deleteOnExit();
			SAMFileWriterFactory writerFactory = new SAMFileWriterFactory().setCreateIndex(true).setTempDirectory(finalBam.getParentFile());
			unsortedWriter = writerFactory.makeBAMWriter(unModifiedReader.getFileHeader(), false, unsortedFinalBam);
			
			walkBams();
			
//IO.pl(id+ "\t# Saved from ModMix: "+ numModifiedSaved+"  UnMod: "+ numUnModifiedSaved);

			//add in the remaining alignments not associated with any vcf
			addInRemainingAlignments();
			
			//close IO
			modifiedIterator.close();
			unModifiedIterator.close();
			modifiedReader.close();
			unModifiedReader.close();
			unsortedWriter.close();
			
			//sort and save final
			IO.pl(id+"\tCoordinate sorting final bam...");
			new PicardSortSam (unsortedFinalBam, finalBam, true);
			
		} catch (Exception e) {
			failed = true;
			e.printStackTrace();
			Misc.printErrAndExit(id+"\tProblem concatinating sams or sorting sams ");
		} 
	}

	private void addInRemainingAlignments() throws IOException {
		IO.pl(id+"\tWriting out alignments not associated with any vcf...");
		SamReader filteredReader = readerFactory.open(filteredBam);
		SAMRecordIterator it = filteredReader.iterator();
		while (it.hasNext()) {
			unsortedWriter.addAlignment(it.next());
			numUnModifiedSaved++;
		}
		it.close();
		filteredReader.close();
	}

	private void walkBams() throws IOException {
			
			//walk every unmodified record
//int counter = 0;
			while (true) {
//if (counter++ == 10) return;
				
				//fetch next unmodified set of Sams
				String unModifiedReadName = loadNextUnmodified();
//if(verbose) IO.pl("xxxxx NewRound" +unModifiedReadName+" with unMod "+unModifiedSamRecords.size());
				
				//any records?
				if (unModifiedReadName == null) {
					writeOutAllModified();
					return;
				}
				
				//attempt to load modified records
				loadModified(unModifiedReadName);
//IO.pl("\tNum Match Mod "+modifiedSamRecords.size());
				
				//any matches?
				if (modifiedSamRecords.size() !=0) writeSamRecords();
				else {
					//write out the unmodified 
					numUnModifiedSaved+=unModifiedSamRecords.size();
					for (SAMRecord unMod: unModifiedSamRecords) unsortedWriter.addAlignment(unMod);
				}
					
				//any modified left?
				if (lastModifiedSamRecord == null) {
					writeOutAllUnModified();
					return;
				}
			}
	}

	private void writeSamRecords() throws IOException {
		//check numbers, should have 1 or 2 records
		int numMod = modifiedSamRecords.size();

		//more than two modified? rarely the mod will have more provided both the single and paired alignments realigned the same read
		if (numMod > 2) {
			ArrayList<SAMRecord> pair = fetchPair(modifiedSamRecords);
			if (pair == null) {
				StringBuilder sb = new StringBuilder();
				for (SAMRecord sam : modifiedSamRecords) sb.append("Md\t"+sam.getSAMString()+"\n");
				throw new IOException ("Failed to select just two modified alignments "+numMod+"\n"+sb);
			}
			else {
				modifiedSamRecords = pair;
				numMod =2;
			}
		}
		
		int numUn = unModifiedSamRecords.size();
		if (numMod>2 || numUn>2 || numMod==0 || numUn==0) {
			StringBuilder sb = new StringBuilder();
			for (SAMRecord sam : modifiedSamRecords) sb.append("Md\t"+sam.getSAMString()+"\n");
			for (SAMRecord sam : unModifiedSamRecords) sb.append("Un\t"+sam.getSAMString()+"\n");
			throw new IOException ("Incorrect pairing "+numMod+" : "+numUn+"\n"+sb);
		}

//boolean speak = true;
		//if(speak) IO.pl("\nNew ReadSet "+modifiedSamRecords.get(0).getReadName());

//if(speak) IO.pl("\t\tNumMod: "+numMod+"   NumUnMod"+numUn);
		
		numModifiedSaved+= numMod;

		//two mods? then write both and return
		if (numMod == 2) {
//if(speak) IO.pl("\t\tTwo mod present, saving both and returning");
			for (SAMRecord mod: modifiedSamRecords) unsortedWriter.addAlignment(mod);
			return;
		}

		//ok only one mod and there is one or more unMod
		SAMRecord modOne = modifiedSamRecords.get(0);
		//save it
		unsortedWriter.addAlignment(modOne);
		
		SAMRecord unModOne = unModifiedSamRecords.get(0);
		
		//check to see if modOne and unModOne are different, if so then save unModOne and return
//if(speak) IO.pl("\t\tChecking single modOne against unModOne");
//if(speak) IO.pl("\t\tMo "+modOne.getSAMString());
		if (matchSAMs(modOne, unModOne) == false) {
//if(speak) IO.pl("\t\tDifferent so writing unModOne");	
//if(speak) IO.pl("\t\tUn "+unModOne.getSAMString());
			unsortedWriter.addAlignment(unModOne);
			numUnModifiedSaved++;
			return;
		}
		else {
//if(speak) IO.pl("\t\tSame so looking to see if there are 2 unMods");	
		}
		
		//print unModTwo if it exists
		if (numUn == 2) {
//if(speak) IO.pl("\t\tWriting 2nd unMod, should be different");
//if(speak) IO.pl("\t\tUn "+unModifiedSamRecords.get(1).getSAMString());
			unsortedWriter.addAlignment(unModifiedSamRecords.get(1));
			numUnModifiedSaved++;
		}
	}

	private ArrayList<SAMRecord> fetchPair(ArrayList<SAMRecord> mergedSinglePairSams) {
		//split by readGroup
		HashMap<String, ArrayList<SAMRecord>> rgs = new HashMap<String, ArrayList<SAMRecord>>();
		for (SAMRecord sam: mergedSinglePairSams) {
			String name = fetchReadGroup(sam);
			ArrayList<SAMRecord> al = rgs.get(name);
			if (al == null) {
				al = new ArrayList<SAMRecord>();
				rgs.put(name, al);
			}
			al.add(sam);
		}
		//return first readGroup with 2
		for (ArrayList<SAMRecord> al: rgs.values()) {
			if (al.size()==2) return al;
		}
		
		return null;
	}

	private String fetchReadGroup(SAMRecord sam) {
		//this is sometimes null from pulling directly due to header issue?
		String[] fields = Misc.TAB.split(sam.getSAMString());
		for (String f: fields) if (f.startsWith("RG:Z:")) return f;
		return null;
	}

	private boolean matchSAMs(SAMRecord one, SAMRecord two) {
		//assumes same read name
		//different strands?
		if (one.getReadNegativeStrandFlag() != two.getReadNegativeStrandFlag()) return false;
		//different unclipped starts
		if (one.getUnclippedStart() != two.getUnclippedStart()) return false;
		//different chroms
		if ((one.getReferenceName().equals(two.getReferenceName())) == false) return false;
		return true;
	}

	private void loadModified(String unModifiedReadName) {
		if (lastModifiedSamRecord == null) return;
		
		//clear old
		modifiedSamRecords.clear();
		
		//is the last one what is wanted?
		String lastModifiedReadName = lastModifiedSamRecord.getReadName();
		if (lastModifiedReadName.equals(unModifiedReadName) == false) return;
		
		//ok it is! add old one
		modifiedSamRecords.add(lastModifiedSamRecord);
		lastModifiedSamRecord = null;
		while (modifiedIterator.hasNext()) {
			SAMRecord nextMod = modifiedIterator.next();
			//same name?
			if (nextMod.getReadName().equals(lastModifiedReadName)) modifiedSamRecords.add(nextMod);
			//nope diff
			else {
				lastModifiedSamRecord= nextMod;
				return;
			}
		}
		
		return;
	}

	/**This will prob never get called.*/
	private void writeOutAllModified() {
		IO.pl(id+ "\tWriting out all remaining modified...");
		if (lastModifiedSamRecord!= null) {
			unsortedWriter.addAlignment(lastModifiedSamRecord);
			numModifiedSaved++;
		}
		while (modifiedIterator.hasNext()) {
			unsortedWriter.addAlignment(modifiedIterator.next());
			numModifiedSaved++;
		}
	}
	
	private void writeOutAllUnModified() {
		IO.pl(id+"\tWriting out all remaining unmodified...");
		if (lastUnModifiedSamRecord!= null) {
			unsortedWriter.addAlignment(lastUnModifiedSamRecord);
			numUnModifiedSaved++;
		}
		while (unModifiedIterator.hasNext()) {
			unsortedWriter.addAlignment(unModifiedIterator.next());
			numUnModifiedSaved++;
		}
	}

	private String loadNextUnmodified() {
		
		//clear old and add last if it isn't null;
		unModifiedSamRecords.clear();
		if (lastUnModifiedSamRecord == null) {
//IO.pl("LastUnModified is null!");
			return null;
		}
		unModifiedSamRecords.add(lastUnModifiedSamRecord);
		
		String currentName = lastUnModifiedSamRecord.getReadName();
		lastUnModifiedSamRecord = null;
		
		//fetch new
		while (unModifiedIterator.hasNext()) {
			SAMRecord nextUnMod = unModifiedIterator.next();
//IO.pl("\tFetching next "+nextUnMod.getReadName());			
			//same name?
			if (nextUnMod.getReadName().equals(currentName)) {
				unModifiedSamRecords.add(nextUnMod);
//IO.pl("\t\tSame as current adding");
			}
			//diff name
			else {
				lastUnModifiedSamRecord = nextUnMod;
//IO.pl("\t\tDifferent than current returning");
				return currentName;
			}
		}
		return currentName;
		

	}

	public boolean isFailed() {
		return failed;
	}

}
