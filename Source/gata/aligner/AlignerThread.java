package gata.aligner;

import gata.main.*;

import java.util.*;
import java.util.regex.*;
import javax.swing.*;

import util.bio.parsers.*;
import util.gen.*;

import java.io.*;

/**
 * @author nix
 * Threaded version of primary GATAligner methods to actually perform and save the alignment.
 * Updates the progress monitor. 
 * */
public class AlignerThread extends Thread {
	ProgressMonitor pm;
	AlignerInputForm aif;
	
	public AlignerThread(ProgressMonitor pm, AlignerInputForm aif) {
		this.pm = pm;
		this.aif = aif;
	}
	public void run() {
		AlignerPreferences alignerPrefs = aif.getPanel().getAlignerPrefReference();
		
		//fetch sequences
		MultiFastaParser ref = new MultiFastaParser(alignerPrefs.getRefSeq());
		MultiFastaParser comp = new MultiFastaParser(alignerPrefs.getCompSeq());
		
		//check if fasta files are OK
		if (ref.isFastaFound() == false || comp.isFastaFound()==false ){
			JOptionPane.showMessageDialog(null,
					"Looks like one of your FASTA files is misformatted!\nCheck both and resubmit.",
					null,JOptionPane.WARNING_MESSAGE);
			pm.close();
			return;
		}
		
		//for each comp seq, make a AlignParms and Seqs objects and fire blast...
		int numCompSeqs = comp.getNumReads();
		int progress = 0;
		for (int i=0; i<numCompSeqs; i++){
			//make new AlignPrams and seqs
			AlignParams ap =new AlignParams(alignerPrefs);
			Seqs seqs = new Seqs(ap, ref, comp, i);
			
			//attempt to extract start position from names?
			if (ap.isEXTRACT()){
				Pattern pat = Pattern.compile("\\d+"); //only returns ints, not floaties
				Matcher matRef = pat.matcher(ref.getNames()[0]);
				Matcher matComp = pat.matcher(comp.getNames()[i]);
				if (matRef.find()){
					ap.setSTART_INDEX_REFSEQ(Integer.parseInt(matRef.group()));
					aif.getPanel().getRefTexF().setText(matRef.group());
				}
				if (matComp.find()){
					ap.setSTART_INDEX_COMPSEQ(Integer.parseInt(matComp.group()));
					aif.getPanel().getCompTexF().setText(matComp.group());
				}
			}
			
			int pair = i+1;
			
			//Blast Seqs
			pm.setNote("# "+pair+" Launching bl2seq/BLASTN...");
			progress += 10/numCompSeqs; //10
			pm.setProgress(progress);
			
			//write comp seq to temp fileif there is more than one
			String compSeqFileName;
			File tempCompFile = null;
			if (numCompSeqs==1) compSeqFileName = ap.getCompSeqFile();
			else {
				String path = new File (ap.getCompSeqFile()).getParent();
				tempCompFile = new File (path, "GATA"+new Random().nextInt(10000)+".working.tmp");
				String contents = ">"+comp.getNames()[i]+"\n"+comp.getSeqs()[i];
				IO.writeString(contents,tempCompFile);
				compSeqFileName = tempCompFile.getPath();
			}
			
			Blast b = new Blast(ap);
			String[] blastResults = b.blastIt(compSeqFileName);
			//check to see if anything came back from BLAST
			if (blastResults.length<20){
				//blast problem!
				GATAUtil.throwWarning(aif, aif.getPanel().getBl2F(), "Hmm, nothing came back from the bl2seq program!\nThe problem may be due to a lack of homology between your sequences.\nTry using the testDemo sequences.\n\nAlternatively, there may be something wrong with the bl2seq BLAST executable.\nTry a different version or download a new copy from \nftp://ftp.ncbi.nih.gov/blast/executables/LATEST-BLAST/\ne.g. blast-2.2.6-powerpc-macosx.tar.gz if on a macOSX");
			}
			
			//delete temp file
			if (tempCompFile!=null) tempCompFile.delete();
			
			ArrayList parsedBlast = b.parseBlast(blastResults);
			
			//make local alignments
			progress += 20/numCompSeqs; //30
			if (setProgMon("# "+pair+" Processing Local Alignments...", progress, pm) == false) return;
			LocalAlignment[] localAligns = seqs.makeLocalAlignments(parsedBlast);
			
			//make sub alignments
			progress += 10/numCompSeqs; //40
			if (setProgMon("# "+pair+" Making Sub Alignments...", progress, pm) == false) return;
			Alignment[] subAlignments = GATAUtil.makeSubAligns(localAligns);
			
			//save AP and alignments
			progress += 50/numCompSeqs; //90
			if (setProgMon("# "+pair+" Saving Alignment objects...", progress, pm) == false) return;
			//make fileName
			String fileName = ap.getBaseName();
			if (numCompSeqs!=1){
				String compName = ap.getNAME_COMPSEQ();
				compName = compName.replaceAll("\\s+","_");
				if (compName.length()>10) compName = compName.substring(0,10);
				fileName = fileName+pair+"_"+compName;
			}
			File file = new File (ap.getPathToResults(), fileName+".gata");
			ArrayList objects = new ArrayList(2);
			objects.add(ap);
			objects.add(subAlignments);
			GATAUtil.saveArrayList(file, objects);
			
			if (i+1 == numCompSeqs){
				setProgMon("Complete!", 99, pm);
				try {
					sleep(500);
				} catch (InterruptedException e) {
					pm.close();
					return;
				}
			}
		}
		pm.close();
	}
	public boolean setProgMon(
			String newNote,
			int progress,
			ProgressMonitor pm) {
		if (pm.isCanceled()) {
			pm.close();
			return false;
		}
		pm.setNote(newNote);
		pm.setProgress(progress);
		return true;
	}
}
