package gata.plotter;

import gata.aligner.*;
import gata.main.*;

import java.util.*;
import java.io.*;
import javax.swing.*;

import util.bio.parsers.*;

/**
 * @author Nix
 Takes an array of Alignments and converts it into two conserved text sequences
 */
public class ConservedSeqs {
	private static JFileChooser chooser;
	private String[] refSeq; //fasta header and sequence
	private String[] compSeq; //ditto
	private AlignPanel panel;
	private Alignment[] aligns;
	private int refIndex;
	private int compIndex;
	private Console console;
	private AlignParams ap;
	private JTextField minBit;
	private JTextField maxBit;
	private JTextField minE;
	private JTextField maxE;
	private ToolsFrame tf;
	private String formattedRef;
	private String formattedComp;

	public ConservedSeqs(GATAParams params) {
		ap = params.getAlignParams();
		refIndex = ap.getSTART_INDEX_REFSEQ();
		compIndex = ap.getSTART_INDEX_COMPSEQ();
		panel = params.getAlignPanel();
		console = new Console (50, 50, 500, 500, "Conserved Sequences");
		JTextArea ta = console.getTextArea();
		ta.setEditable(true);
		ta.setLineWrap(true);
		tf = params.getToolsFrameRef();
		minBit = tf.getMinBit();
		maxBit = tf.getMaxBit();
		minE = tf.getMinE();
		maxE = tf.getMaxE();
	}
	
	public boolean fetchSequences(){	
		if (refSeq==null)refSeq = fetchSeq(ap.getRefSeqFile());
		if (compSeq== null)compSeq = fetchSeq(ap.getCompSeqFile());
		if (refSeq!=null && compSeq!=null) return true;
		return false;
	}
	
	public void fetchNewRegions(int leftBaseRef, int rightBaseRef, int leftBaseComp, int rightBaseComp){
		//reset title
		console.setTitle("Conserved Seqs: Ref("+(leftBaseRef+refIndex)+"-"+(rightBaseRef+refIndex-1)+
			") Comp("+(leftBaseComp+compIndex)+"-"+(rightBaseComp+compIndex-1)+")");
		generateConservedSeqs();
		printTruncatedSeqsToConsole(leftBaseRef, rightBaseRef, leftBaseComp, rightBaseComp);
	}
		
	public void generateConservedSeqs() {
		if (fetchSequences()==false) return;
		aligns = panel.getVisableAlignments();
		//convert Alignment[] into to int[][] representing the ref seq's start
		//	stops and the comp seq's start stops
		ArrayList ref = new ArrayList();
		ArrayList comp = new ArrayList();

		int len = aligns.length;
		for (int i = 0; i < len; i++) {
			int[] coor = aligns[i].getCoordinates();
			//check if +/-
			if (coor[0] < coor[1])
				ref.add(new StartStop(coor[0] - refIndex, coor[1] - refIndex));
			else
				ref.add(new StartStop(coor[1] - refIndex, coor[0] - refIndex));
			if (coor[2] < coor[3])
				comp.add(
					new StartStop(coor[2] - compIndex, coor[3] - compIndex));
			else
				comp.add(
					new StartStop(coor[3] - compIndex, coor[2] - compIndex));
		}
		StartStop[] refSS = new StartStop[ref.size()];
		StartStop[] compSS = new StartStop[comp.size()];
		ref.toArray(refSS);
		comp.toArray(compSS);

		//convert to blocks
		int[][] refBlocks = convertToBlocks(refSS);
		int[][] compBlocks = convertToBlocks(compSS);

		//request how the user would like the seq's reformatted
		String format =
			(String) JOptionPane.showInputDialog(
				null,
				"Choose a means to reformat the non-conserved sequence blocks.  Enter\n   'UPPERCASE', 'lowercase', or a single character (e.g. 'x' or 'n')",
				"Conserved Sequence Formatting",
				JOptionPane.PLAIN_MESSAGE,
				null,
				null,
				"lowercase");
		if (GATAUtil.isStringEmpty(format))
			format = "lowercase"; //if user hits cancel or an empty text

		formattedRef = reformatSeq(refBlocks, refSeq[1], format.trim());
		formattedComp = reformatSeq(compBlocks, compSeq[1], format.trim());
	}
	
	public void printSeqsToConsole(){
		String header = "(Conserved Seqs: Min Score "+minBit.getText()+"bits "+minE.getText()+" Expect; "+
			"Max Score "+maxBit.getText()+"bits "+maxE.getText()+"Expect Window: "+ap.getWIN_SIZE()+")\n";
		console.setTextArea(">"+refSeq[0]+ header + formattedRef);
		console.printToTextArea("\n\n>"+compSeq[0]+ header + formattedComp);
		console.show();

	}
	
	public void printTruncatedSeqsToConsole(int leftBaseRef, int rightBaseRef, int leftBaseComp, int rightBaseComp){
		//trim formatted seqs
		String subRef = formattedRef.substring(leftBaseRef, rightBaseRef);
		String subComp = formattedComp.substring(leftBaseComp, rightBaseComp);

		String head = "(Conserved Seqs: Min Score "+minBit.getText()+"bits "+minE.getText()+" Expect; "+
			"Max Score "+maxBit.getText()+"bits "+maxE.getText()+"Expect Window: "+ap.getWIN_SIZE()+")";

		console.setTextArea(">"+refSeq[0]+ head+
			" ("+(leftBaseRef+refIndex)+"-"+(rightBaseRef+refIndex-1)+")\n" + subRef);
		console.printToTextArea("\n\n>"+compSeq[0]+ head+
			" ("+(leftBaseComp+compIndex)+"-"+(rightBaseComp+compIndex-1)+")\n" + subComp);
		console.show();
	}	

	/**Reformats a sequence given and int[][] containing start stops*/
	public static String reformatSeq(
		int[][] blocks,
		String seq,
		String format) {
		String formattedSeq;
		int len = blocks.length;	

		if (format.equalsIgnoreCase("lowercase")) {
			StringBuffer sb = new StringBuffer(seq.toLowerCase());
			for (int i = 0; i < len; i++) {
				String sub = sb.substring(blocks[i][0], blocks[i][1] + 1);
				sub = sub.toUpperCase();
				sb.replace(blocks[i][0], blocks[i][1] + 1, sub);
			}
			return sb.toString();
		} else if (format.equalsIgnoreCase("uppercase")) {
			StringBuffer sb = new StringBuffer(seq.toUpperCase());
			for (int i = 0; i < len; i++) {
				String sub = sb.substring(blocks[i][0], blocks[i][1] + 1);
				sub = sub.toLowerCase();
				sb.replace(blocks[i][0], blocks[i][1] + 1, sub);
			}
			return sb.toString();
		} else { //replace with a char
			char c = format.charAt(0);
			StringBuffer sb = new StringBuffer(len);
			int seqLen = seq.length();
			for (int i = 0; i < seqLen; i++) {
				sb.append(c);
			}
			for (int i = 0; i < len; i++) {
				String sub = seq.substring(blocks[i][0], blocks[i][1] + 1);
				sb.replace(blocks[i][0], blocks[i][1] + 1, sub);
			}
			return sb.toString();
		}
	}

	public static String[] fetchSeq(String file) {
		boolean ok = false;
		//check to see if files exist, if so fetch first seq
		File refFile = new File(file);
		if (refFile.canRead()) {
			MultiFastaParser mfp = new MultiFastaParser(file);
			if (mfp.isFastaFound())
				return new String[] { mfp.getNames()[0], mfp.getSeqs()[0] };
		}
		//if this has failed to return, throw warning
		JOptionPane.showMessageDialog(
			null,
			"Sorry, "
				+ file
				+ " cannot be located or properly parsed!\n   Please find the FASTA formatted sequence file used in creating the original alignment.",
			null,
			JOptionPane.WARNING_MESSAGE);
		//get file
		if (chooser == null)
			chooser = new JFileChooser();
		chooser.setCurrentDirectory(refFile);
		int result = chooser.showDialog(null, "Choose");
		if (result == JFileChooser.APPROVE_OPTION) {
			File newfile = chooser.getSelectedFile();
			try {
				return fetchSeq(newfile.getCanonicalPath());
			} catch (IOException IOe) {
				JOptionPane.showMessageDialog(
					null,
					"Problem reading file! ",
					null,
					JOptionPane.WARNING_MESSAGE);
				IOe.printStackTrace();
				return null;
			}
		} else
			return null;
		//the user didn't use chooser and hit cancel			
	}

	public static int[][] convertToBlocks(StartStop[] ss) {
		//sort arrays based on start position
		Arrays.sort(ss);

		//extract int[][] arrays
		int len = ss.length;
		int[][] blocks = new int[len][2];
		for (int i = 0; i < len; i++) {
			int[] strStp = ss[i].getStartStop();
			blocks[i] = new int[] { strStp[0], strStp[1] };			
		}

		//combine start stops into blocks
		return growIntArrays(blocks);
	}

	public static int[][] growIntArrays(int[][] ints) {
		int len = ints.length;
		if (len <= 1)
			return ints; //if only one block!	
		ArrayList bigInts = new ArrayList();
		int[] seed = new int[] { ints[0][0], ints[0][1] };
		int lenMinus = len - 1;
		for (int i = 1; i < len; i++) {
			int[] second = new int[] { ints[i][0], ints[i][1] };
			//first check if they overlap
			if (seed[1] >= second[0]
				|| seed[1] >= second[0] - 1) { //overlap or adjacent
				//check if stop is bigger
				if (second[1] > seed[1]) { //second has bigger but
					seed[1] = second[1];
					//save last one
					if (i == lenMinus)
						bigInts.add(seed);
				}
				if (i == lenMinus)
					bigInts.add(seed);
				//if not then skip
			} else { //they dont overlap so save seed
				bigInts.add(seed);
				seed = second;
				//save last one
				if (i == lenMinus)
					bigInts.add(seed);
			}

		}
		int[][] goods = new int[bigInts.size()][2];
		bigInts.toArray(goods);
		return goods;
	}
	public Console getConsole() {
		return console;
	}

}