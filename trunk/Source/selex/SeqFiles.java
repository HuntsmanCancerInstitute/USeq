package selex;
import java.util.ArrayList;

import util.bio.parsers.*;


/**Continer for raw sequence information, seqFile, qualFile, and the actual reads from each.*/
public class SeqFiles {
    private String seqFile;
    private String qualFile;
    private int numReads;
    private SeqRead[] seqReads;
    private SelexParams sp;
    

    public SeqFiles(String aSeqFile, String aQualFile, SelexParams SP) {
        
        //make field assigments
        seqFile = aSeqFile;
        qualFile = aQualFile;
        sp = SP;
        
        //extract sequences and comment lines from seqFile
        MultiFastaParser mfp = new MultiFastaParser(seqFile);
        numReads = mfp.getNumReads();
        String[] seqs = mfp.getSeqs();
        String[] sComs = mfp.getNames();
        
        //extract qual data and comment lines from qualFile
        QualityFileParser qfp = new QualityFileParser(qualFile);
        String[] quals = qfp.getQuals();
        String[] qComs = qfp.getNames();
        
        //check that number of seqs = number of quals
        if (seqs.length != quals.length || sComs.length!=qComs.length){
            System.out.println("");
            System.out.println("Fatal Problem! The number of sequence reads in "+seqFile);
            System.out.println("    does not match the number of quality reads in "+qualFile);
            System.out.println("    Is this the real quality file for this sequence file?");
            System.out.println("");
            System.exit(0);
        }
        printSeqFiles();
        
        //modify ends? For fixing a blunted ligation ends
        if (sp.modifyEnds()) {
        		modifySequenceEnds(seqs, quals, sp.getLeftSideMatch(), sp.getLeftSideReplace(), sp.getRightSideMatch(), sp.getRightSideReplace());
        		System.out.println("Attempting to modify sequence ends....");
        }
        
        //make seqRead objects
        seqReads = new SeqRead[numReads];
        for (int i=0; i<numReads; i++){
            seqReads[i] = new SeqRead(seqs[i], sComs[i], quals[i], qComs[i], sp, this);
        }
        sp.incNumSeqReads(numReads);
    }
    
    
    
    //getter methods
    public String getSeqFile(){return seqFile;}
    public String getQualFile(){return qualFile;}
    public int getNumReads(){return numReads;}
    public SeqRead[] getSeqReads(){return seqReads;}
    
    //primary methods
    public void printSeqFiles(){
        sp.printSave(
        "***************New Sequence File Set***************\n"+
        "Name of Sequence File: "+seqFile+
        "\nName of Quality File: "+qualFile+
        "\nNumber of Sequence Reads in File: "+numReads+"\n\n\n");
    }
    
	/**Modifies the left side and right side of Sequences also adds (100's) or subtracts quanilty scores to match.
	 * Will modify the first occurance on the left, and the last occurance on the right.
	 * Usefull for fixing blunted sequence ends for Selex parsing.
	 * Case insensitive.*/
	public static void modifySequenceEnds(String[] seqs, String[] quals, String leftSideMatch, String leftSideReplace, String rightSideMatch, String rightSideReplace){
		//case convert sequences
		String lMatch = leftSideMatch.toLowerCase();
		String lReplace = leftSideReplace.toLowerCase();
		String rMatch = rightSideMatch.toLowerCase();
		String rReplace = rightSideReplace.toLowerCase();
		
		//calc lengths
		int lengthRMatch = rMatch.length();
		int lengthLMatch = lMatch.length();
		int diffL = lReplace.length() - lengthLMatch;
		int diffR = rReplace.length() - lengthRMatch;
		
		//run thru each sequence
		int num = seqs.length;
		for (int i=0; i<num; i++){
			seqs[i] = seqs[i].toLowerCase();
			//modify 1st occurance on left side
			int index= seqs[i].indexOf(lMatch);
			if (index != -1){
				//change sequence
				seqs[i] = seqs[i].substring(0,index) + lReplace + seqs[i].substring(index+lengthLMatch);
				//change quality String
				quals[i] = changeQuality(quals[i], index+ lengthLMatch -1, diffL);
			}
			//modify last occurance on right side
			index = seqs[i].lastIndexOf(rMatch);
			if (index != -1){
				seqs[i] = seqs[i].substring(0,index) + rReplace + seqs[i].substring(index+lengthRMatch);
				//change quality String
				quals[i] = changeQuality(quals[i], index+ lengthRMatch -1, diffR);
			}
		}
	}
	/**Modifies a quality file, adding or subtracting scores.*/
	public static String changeQuality(String qual, int index, int diff){
		if (diff == 0) return qual;
		String[] qs = qual.split("\\s+");
		ArrayList al = stringArrayToArrayList(qs);
		if (diff > 0) {
			//add scores
			int num = index+diff;
			for (int i=index; i<num; i++){
				al.add(index,"100");
			}
		}
		else{
			//remove scores
			int num = index - diff;
			for (int i=index; i<num; i++){
				al.remove(index);
			}
		}
		return stringArrayListToString(al," ");
	}
	
	/**Returns an ArrayList of String given a String[]*/
	public static ArrayList stringArrayToArrayList(String[] s){
		int num = s.length;
		ArrayList al = new ArrayList(num);
		for (int i=0; i<num; i++){
			al.add(s[i]);
		}
		return al;
	}

	/**Returns a String separated by the separator given an ArrayList of String.*/
	public static String stringArrayListToString(ArrayList stringAL, String separator){
		int len = stringAL.size();
		if (len==0) return "";
		if (len==1) return (String)stringAL.get(0);
		StringBuffer sb = new StringBuffer((String)stringAL.get(0));
		for (int i=1; i<len; i++){
			sb.append(separator);
			sb.append((String)stringAL.get(i));
		}
		return sb.toString();
	}
}