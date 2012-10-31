package gata.main;

import gata.aligner.*;

import java.util.*;
import java.util.regex.*;
import java.io.*;
import java.awt.geom.*;
import java.awt.*;
import java.awt.font.*;
import java.text.*;
import javax.swing.*;

import edu.utah.seq.useq.data.Region;

import util.bio.annotation.*;

import java.awt.event.*;


/**
 * @author Nix
 * Master utility class for GATAPlotter.  All static.
 */
public class GATAUtil {	
	
	/**Takes an array of boolean[] and sets all to the boolean value*/
	public static void setBooleanArray (boolean[] array, boolean value){
		for (int i= array.length-1; i>=0; i--) array[i] = value;					
	}
	
	/**Takes an array of float[] and sets all to the float value*/
	public static void setFloatArray (float[] array, float value){
		for (int i= array.length-1; i>=0; i--) array[i] = value;					
	}
	
	
	/**Sets the scroll bar move increments, vert and horiz, in a scroll pane*/
	public static void setScrollPaneScrollBarIncrement(JScrollPane scrollPane, int increment){
		scrollPane.getHorizontalScrollBar().setUnitIncrement(increment);
		scrollPane.getVerticalScrollBar().setUnitIncrement(increment);
	}
	
	public static double makeDoublePercentInputDialog(
			String title,
			String comment,
			double suggestion) {
		Object ob =
			JOptionPane.showInputDialog(null,comment,title,JOptionPane.PLAIN_MESSAGE,null,null,Double.toString(suggestion));
		if (ob == null) { // if user hits cancel or closes the window
			return suggestion;
		}
		double num;
		try{
			num = Double.parseDouble((String) ob);
		}
		catch (NumberFormatException e){
			num = makeDoublePercentInputDialog(title,"Decimals Only! "+comment,suggestion);
		}
		if (num >1 || num<0) return makeDoublePercentInputDialog(title,"Fractional Decimal Only (0-1)! "+comment,suggestion);
		else  return num;
	}
	
	public static double makeDoubleInputDialog(
			String title,
			String comment,
			double suggestion) {
		Object ob =
			JOptionPane.showInputDialog(null,comment,title,JOptionPane.PLAIN_MESSAGE,null,null,Double.toString(suggestion));
		if (ob == null) { // if user hits cancel or closes the window
			return suggestion;
		}
		double num;
		try{
			num = Double.parseDouble((String) ob);
		}
		catch (NumberFormatException e){
			num = makeDoubleInputDialog(title,"Problem parsing number! "+ob+" "+comment,suggestion);
		}
		return num;
	}
	
	public static int makeIntInputDialog(String title,String comment,int suggestion) {
		Object ob =JOptionPane.showInputDialog(null,comment,title,JOptionPane.PLAIN_MESSAGE,null,null,Integer.toString(suggestion));
		if (ob == null) { // if user hits cancel or closes the window
			return suggestion;
		}
		try{
			return Integer.parseInt((String) ob);
		}
		catch (NumberFormatException e){
			return makeIntInputDialog(title,"Integers only! "+comment,suggestion);
		}
		
	}
	
	
	public static void throwWarning(
			JFrame frame,
			JTextField field,
			String text) {
		JOptionPane.showMessageDialog(
				frame,
				text,
				null,
				JOptionPane.WARNING_MESSAGE);
		field.requestFocusInWindow();
	}
	
	public static int[] setButton(JButton but, int x, int y) {
		Dimension dim = but.getPreferredSize();
		int[] wh = {(int) dim.getWidth(), (int) dim.getHeight()};
		but.setBounds(x, y, wh[0], wh[1]);
		return wh;
	}
	
	public static JButton makeButton(
			String text,
			int size,
			ActionListener action,
			JPanel panel) {
		
		JButton but = new JButton(text);
		but.setFont(new Font("Dialog", 0, size));
		but.setBackground(Color.WHITE);
		panel.add(but);
		
		but.addActionListener(action);
		return but;
	}
	
	public static void showFileChooserDialogBox(
			JTextField fieldRef,
			JFrame frame,
			JFileChooser chooser) {
		//int result = chooser.showDialog(frame, "Choose"); //problems with selecting folders!
		int result = chooser.showOpenDialog(frame);
		if (result == JFileChooser.APPROVE_OPTION) {
			File file = chooser.getSelectedFile();
			String text = "Problem fetching this file/directory!";
			try {
				text = file.getCanonicalPath();
			} catch (IOException IOe) {
				System.out.println("Problem with getting path/file info!");
				IOe.printStackTrace();
			}
			fieldRef.setText(text);
			chooser.setCurrentDirectory(file);
		}
		fieldRef.requestFocusInWindow();
	}
	
	public static boolean checkFile(String par) {
		File f = new File(par);
		if (f.exists())
			return true;
		if (f.isDirectory())
			return true;
		JOptionPane.showMessageDialog(
				null,
				"Sorry, I cannot find ' "+ par+" '!",
				null, JOptionPane.WARNING_MESSAGE);
		return false;
	}
	public static JTextField makeField(
			JPanel panel,
			String text,
			JFrame frame,
			JFileChooser chooser,
			FileFocusListener focus) {
		JTextField tex = new JTextField(text);
		tex.addFocusListener(focus);
		focus.setRefs(tex, frame, chooser);
		panel.add(tex);
		return tex;
	}
	public static JTextField makeFieldSimple(
			JPanel panel,String text) {
		JTextField tex = new JTextField(text);
		panel.add(tex);
		return tex;
	}
	public static JLabel makeLabel(
			JPanel panel,
			String text,
			int size,
			int thickness) {
		JLabel lab = new JLabel(text);
		lab.setFont(new Font("Dialog", thickness, size));
		panel.add(lab);
		return lab;
	}
	public static int[] setChooser(
			JFileChooser choos,
			int x,
			int y,
			int width) {
		Dimension dim = choos.getPreferredSize();
		int[] wh = {(int) dim.getWidth(), (int) dim.getHeight()};
		choos.setBounds(x, y, wh[0], wh[1]);
		return wh;
	}
	public static int[] setField(JTextField field, int x, int y, int width) {
		Dimension dim = field.getPreferredSize();
		int[] wh = { width, (int) dim.getHeight()};
		field.setBounds(x, y, wh[0], wh[1]);
		return wh;
	}
	public static int[] setLabel(JLabel lab, int x, int y) {
		Dimension dim = lab.getPreferredSize();
		int[] wh = {(int) dim.getWidth(), (int) dim.getHeight()};
		lab.setBounds(x, y, wh[0], wh[1]);
		return wh;
	}
	public static void deleteFile(String fullPathFileName) {
		try {
			File file = new File(fullPathFileName);
			if (!file.delete())
				System.out.println(
						"Warning: could not delete the file " + fullPathFileName);
		} catch (SecurityException e) {
			System.out.println(
					"Warning: could not delete the file " + fullPathFileName);
			e.printStackTrace();
		}
	}
	
	/**Converts a raw S score to a bit score given lambda (in nats) and K*/
	public static double convertRawScoreToBit(
			double lambdaNats,
			double K,
			double rawS) {
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(1);
		double bs = (lambdaNats * rawS - Math.log(K)) / Math.log(2);
		String bsS = f.format(bs);
		bsS = bsS.replaceAll(",", "");
		return Double.parseDouble(bsS);
	}
	public static double convertBitScoreToRaw(
			double lambdaNats,
			double K,
			double bitS) {
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(1);
		double bs = Math.round((bitS * Math.log(2) + Math.log(K)) / lambdaNats);
		String bsS = f.format(bs);
		bsS = bsS.replaceAll(",", "");
		return Double.parseDouble(bsS);
	}
	public static String calculateEValue(
			double lambdaNats,
			double K,
			int rawS,
			double effectiveN,
			double effectiveM) {
		double num =
			K * effectiveN
			* effectiveM
			* Math.exp(-1 * lambdaNats * (double) rawS);
		DecimalFormat f = new DecimalFormat("0E0");
		String numS = f.format(num);
		if (numS.equals("0E0"))
			return "<1E-99";
		return numS;
		
	}
	/**Calculates effective m and n using BLAST parameters, returning int[m', n']*/
	public static int[] calculateEffectiveMandN (double m, double n, double K, double H){
		double l= (Math.log(K*m*n))/H;
		return new int[]{(int)Math.round(m-l), (int)Math.round(n-l)};
	}
	public static int DNAStat(
			String seq1,
			String seq2,
			int match,
			int misMatch,
			int gapCreate,
			int gapExtend) {
		/*scores a DNA alignment, returns the score */
		Character dash = new Character('-');
		int score = 0;
		boolean newGap = true;
		int matches = 0;
		int len = seq1.length();
		seq1 = seq1.toUpperCase();
		seq2 = seq2.toUpperCase();
		for (int i = 0; i < len; i++) {
			Character s1 = new Character(seq1.charAt(i));
			Character s2 = new Character(seq2.charAt(i));
			int comp = s1.compareTo(s2);
			//check to see if they are the same
			if (comp == 0) {
				score += match;
				newGap = true;
				matches++;
			}
			//if they are not the same check to see if it is a gap
			else {
				int comp1 = s1.compareTo(dash);
				int comp2 = s2.compareTo(dash);
				if (comp1 == 0 || comp2 == 0) { //there is a gap
					//check to see if this is a new gap
					if (newGap) {
						score += gapCreate;
						newGap = false;
					} else
						score += gapExtend;
				} else { //it's a mismatch
					score += misMatch;
					newGap = true;
				}
			}
		}
		return score;
	}
	
	public static int findHighestInt(int[] ints) {
		int len = ints.length;
		int index = 0;
		int score1 = ints[0];
		for (int i = 1; i < len; i++) {
			int score2 = ints[i];
			if (score2 > score1) {
				index = i;
				score1 = score2;
			}
		}
		return index;
	}
	
	public static String genDashes(String seq1, String seq2) {
		int len = seq1.length();
		StringBuffer dashes = new StringBuffer(len);
		for (int i = 0; i < len; i++) {
			Character s1 = new Character(seq1.charAt(i));
			Character s2 = new Character(seq2.charAt(i));
			int comp = s1.compareTo(s2);
			if (comp == 0) {
				dashes.append("|");
			} else {
				dashes.append(" ");
			}
		}
		return new String(dashes);
	}
	
	public static int[] getBaseNumArray(String s, int index) {
		//used to find out real numbers in seq, ignores gaps
		int len = s.length();
		int[] ref = new int[len];
		String gap = "-";
		int counter = index;
		for (int i = 0; i < len; i++) {
			if (s.regionMatches(i, gap, 0, 1) == false) {
				ref[i] = counter;
				counter++;
			}
		}
		return ref;
	}
	
	public static Alignment[] makeSubAligns(LocalAlignment[] locs) {
		int len = locs.length;
		ArrayList subAligns = new ArrayList(len * 5);
		
		for (int i = 0; i < len; i++) {
			//get an arraylist of alignment objects foreach local alignment and merge
			ArrayList subsTemp = locs[i].processLocalAlignment();
			if (subsTemp != null)
				subAligns.addAll(locs[i].processLocalAlignment());
		}
		//convert the ArrayList of Alignments to an array
		len = subAligns.size();
		Alignment[] sortedSubAlignments = new Alignment[len];
		subAligns.toArray(sortedSubAlignments);
		//sort by score, lowest to heighest
		Arrays.sort(sortedSubAlignments);
		return sortedSubAlignments;
	}
	
	public static LocalAlignment[] mergeLAArrays(
			LocalAlignment[] a,
			LocalAlignment[] b) {
		int lenA = a.length;
		int lenB = b.length;
		LocalAlignment[] c = new LocalAlignment[lenA + lenB];
		if (lenA > 0 && lenB > 0) { //bits in both merge!
			System.arraycopy(a, 0, c, 0, lenA);
			System.arraycopy(b, 0, c, lenA, lenB);
		} else if (lenA > 0 && lenB == 0)
			c = a; //no b
		else if (lenA == 0 && lenB > 0)
			c = b; // no a
		return c;
	}
	
	public static int[] optimizeAlignmentScore(
			String seq1,
			String seq2,
			int match,
			int misMatch,
			int gapCreate,
			int gapExtend,
			int scoreCutOff) {
		/**optimizes the score of an alignment by trimming back either stop, returns an ArrayList with the
		 int score, String
		 seq1, String seq2, int index 5', int index 3' where 0 is the first base in the window*/
		//create an array contining the score for each possible alignment resulting from trimming ends
		//trim 5'
		//System.out.println("UNopti :\n"+seq1+"\n"+seq2);
		int minBases = (scoreCutOff / match) - 1;
		//only need to trimback partial
		if (minBases < 0)
			minBases = 0;
		int len = seq1.length();
		int trim = len - minBases;
		if (trim <= 0)
			System.out.println(
			"\n\nYour score cut off is set too high for the given window size!\n\n");
		int[] scores5 = new int[trim];
		for (int i = 0; i < trim; i++) {
			String subSeq1 = seq1.substring(i);
			String subSeq2 = seq2.substring(i);
			scores5[i] =
				DNAStat(
						subSeq1,
						subSeq2,
						match,
						misMatch,
						gapCreate,
						gapExtend);
		}
		//trim 3'
		int[] scores3 = new int[len];
		for (int i = len - 1; i >= minBases; i--) {
			String subSeq1 = seq1.substring(0, i + 1);
			String subSeq2 = seq2.substring(0, i + 1);
			scores3[i] =
				DNAStat(
						subSeq1,
						subSeq2,
						match,
						misMatch,
						gapCreate,
						gapExtend);
		}
		//find highest index
		int index5 = findHighestInt(scores5);
		int index3 = findHighestInt(scores3);
		//find if trimback ran to stop or if score is negative, if so then don't trim
		if (index5 == trim || scores5[index5] <= 0) {
			index5 = 0;
		}
		if (index3 == (len - trim) || scores3[index3] <= 0) {
			index3 = len - 1;
		}
		if (index5 >= index3) { //if messed up coordinates
			index5 = 0;
			index3 = len - 1;
		}
		//make optimized alignment
		String subSeq1 = seq1.substring(index5, index3 + 1);
		String subSeq2 = seq2.substring(index5, index3 + 1);
		
		int optScore =
			DNAStat(subSeq1, subSeq2, match, misMatch, gapCreate, gapExtend);
		int[] x = new int[3];
		x[0] = index5; //start where 0 was pos in alignment
		x[1] = index3; //stop
		x[2] = optScore;
		return x;
	}
	
	public static String parseDNAFileName(String s) {
		// this should pull the actual file text from a full path text on dos and unix
		Pattern pat = Pattern.compile("[-~\\.\\w]+$");
		Matcher mat = pat.matcher(s);
		if (mat.find())
			return mat.group();
		System.out.println(
				"There is a problem with your DNA file text ("
				+ s
				+ ").  Are there any non standard characters associated with it?");
		System.exit(0);
		return "";
	}
	
	
	public static String spaces(String seq, int number) {
		int numSpaces = seq.length() - Integer.toString(number).length() - 1;
		if (numSpaces <= 0)
			return " ";
		StringBuffer sb = new StringBuffer(numSpaces);
		for (int i = 0; i < numSpaces; i++) {
			sb.append(' ');
		}
		return new String(sb);
	}
	
	public static void writeString(String data, String fullPathFileName) {
		try {
			PrintWriter out = new PrintWriter(new FileWriter(fullPathFileName));
			out.println(data);
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	/** Map of hints for high quality/ slow rendering*/
	public static RenderingHints fetchRenderingHints () {
		RenderingHints hints = new RenderingHints(null);
		hints.put(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		hints.put(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		hints.put(RenderingHints.KEY_DITHERING, RenderingHints.VALUE_DITHER_DISABLE);
		hints.put(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
		hints.put(RenderingHints.KEY_FRACTIONALMETRICS, RenderingHints.VALUE_FRACTIONALMETRICS_ON);
		hints.put(RenderingHints.KEY_ALPHA_INTERPOLATION, RenderingHints.VALUE_ALPHA_INTERPOLATION_QUALITY);
		hints.put(RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_QUALITY);
		return hints;
	}
	
	/**Color array for tracks, yup I'm lazy*/
	public static final Color[] COLORS =
		new Color[] {
			Color.MAGENTA,
			Color.GREEN,
			Color.CYAN,
			Color.PINK,
			Color.ORANGE,
			Color.LIGHT_GRAY,
			Color.MAGENTA,
			Color.CYAN,
			Color.GREEN,
			Color.PINK,
			Color.ORANGE,
			Color.LIGHT_GRAY,
			Color.MAGENTA,
			Color.CYAN,
			Color.GREEN,
			Color.PINK,
			Color.ORANGE,
			Color.LIGHT_GRAY,
			Color.MAGENTA,
			Color.CYAN,
			Color.GREEN,
			Color.PINK,
			Color.ORANGE,
			Color.LIGHT_GRAY };
	
	public static void drawLabels(ArrayList labels, Font font, Graphics2D g2) {
		FontRenderContext context = g2.getFontRenderContext();
		int len = labels.size() - 1;
		for (int i = len; i >= 0; i -= 2) {
			double[] nums = (double[]) labels.get(i); //x,y  center pt
			String name = (String) labels.get(i - 1);
			Rectangle2D bounds = font.getStringBounds(name, context);
			double halfWidth = (bounds.getWidth()) / 2;
			double x = nums[0] - halfWidth;
			double height = bounds.getHeight();
			double y = nums[1] + 8 - (height / 2); //fudge
			g2.drawString(name, (float) x, (float) y);
		}
	}
	
	public static ArrayList makeLabel(
			String label,
			int start,
			int end,
			double pixPerNt,
			double ntAtX,
			double yCoor) {
		ArrayList labels = new ArrayList(2);
		int mid = ((1+end - start) / 2) + start;
		double xMidPtPix = (mid - ntAtX) * pixPerNt;
		labels.add(label);
		labels.add(new double[] { xMidPtPix, yCoor });
		return labels;
	}
	
	public static Line2D.Double[] ALToLine2DArray(ArrayList al) {
		int len = al.size();
		Line2D.Double[] lines = new Line2D.Double[len];
		for (int i = len - 1; i >= 0; i--)
			lines[i] = (Line2D.Double) al.get(i);
		return lines;
	}
	
	/**draws Shapes*/
	public static void drawLines(Line2D.Double[] lines, Graphics2D g2) {
		for (int i = lines.length - 1; i >= 0; i--)
			g2.draw(lines[i]);
	}
	
	/**Converts an ArrayList of int[2]'s (startStop) in nucleotides to an ArrayList of Line2D.Double's*/
	public static ArrayList makeSegments(
			ArrayList segs,
			double yCoord,
			double pixPerNt,
			double ntAtX) {
		int len = segs.size();
		ArrayList lines = new ArrayList(len);
		for (int i = len - 1; i >= 0; i--) {
			int[] startStop = (int[]) segs.get(i);
			double s = ((double) startStop[0] - ntAtX) * pixPerNt;
			double e = ((double) startStop[1] - ntAtX) * pixPerNt;
			lines.add(new Line2D.Double(s, yCoord, e, yCoord));
			//System.out.println("Line Segment: "+s+" "+yCoord+" "+e+" "+yCoord);			
		}
		return lines;
	}
	/**Converts an int[x][2]'s (startStop) in nucleotides to an ArrayList of Line2D.Double's*/
	public static ArrayList makeSegmentsInt(
			int[][] segs,
			double yCoord,
			double pixPerNt,
			double ntAtX) {
		int len = segs.length;
		ArrayList lines = new ArrayList(len);
		for (int i = len - 1; i >= 0; i--) {
			double s = (segs[i][0] - ntAtX) * pixPerNt;
			double e = (segs[i][1] - ntAtX) * pixPerNt;
			lines.add(new Line2D.Double(s, yCoord, e, yCoord));
			//System.out.println("Line Segment: "+s+" "+yCoord+" "+e+" "+yCoord);			
		}
		return lines;
	}
	
	/**Converts an int[2] (startStop) in nucleotides to a Line2D.Double*/
	public static Line2D.Double makeLine(
			int start,
			int end,
			double yCoord,
			double pixPerNt,
			double ntAtX) {
		double s = (start - ntAtX) * pixPerNt;
		double e = (end+1 - ntAtX) * pixPerNt;
		return new Line2D.Double(s, yCoord, e, yCoord);
	}
	/**Returns a rectangle given start stop in nucleotides and y coords in pixels*/
	public static Rectangle2D.Double makeRectangle(
			int start, 
			int end,
			double startY,
			double stopY,
			double pixPerNt,
			double ntAtX) {
		double s = (start - ntAtX) * pixPerNt;
		double e = (end+1 - ntAtX) * pixPerNt;
		return new Rectangle2D.Double(s, stopY, e - s, startY - stopY);
	}
	
	/**Makes an arrow and butt but no base line.  Half arrow/butt up for -1 or down for 1.  Start
	 * and stop are really just 5 ' and 3' coordinates, start should always be less
	 * than stop.  Orientation is 1 for + (arrow pointing right), -1 for - (left).
	 * yCoord is for the y-coordinate to draw the main line.  Returns an 
	 * ArrayList of Line2D.doubles.*/
	public static ArrayList makeArrowAndButt(
			int start,
			int end,
			int upOrDown,
			double yCoord,
			int orientation,
			double pixPerNt,
			double ntAtX,
			double height) { //arrow height
		ArrayList lines = new ArrayList(2);
		
		//account for orientation
		double xS;
		double xE;
		if (orientation == 1) {
			xS = pixPerNt * ((double) start - ntAtX);
			xE = pixPerNt * ((double) end + 1 - ntAtX);
		} else {
			xS = pixPerNt * ((double) end +1 - ntAtX);
			xE = pixPerNt * ((double) start - ntAtX);
		}
		//make arrow head
		lines.add(
				new Line2D.Double(
						xE,
						yCoord,
						xE - (height * orientation),
						yCoord + (height * upOrDown)));
		//	  System.out.println("arrow line: "+xE+" "+yCoord+" "+(xE - (height * orientation))+" "+(yCoord + (height * upOrDown)));	
		//make but
		lines.add(
				new Line2D.Double(xS, yCoord, xS, yCoord + (height * upOrDown)));
		//	  System.out.println("butt line: "+xS+" "+yCoord+" "+xS+" "+(yCoord + (height * upOrDown)));				
		return lines;
	}
	
	/**Makes an arrow line with half arrow up for -1 or down for 1.  Start
	 * and stop are really just 5 ' and 3' coordinates, start should always be less
	 * than stop.  Orientation is 1 for + (arrow pointing right), -1 for - (left).
	 * yCoord is for the y-coordinate to draw the main line.  Returns an 
	 * ArrayList of Line2D.doubles.*/
	public static ArrayList makeArrowLine(
			int start,
			int end,
			int upOrDown,
			double yCoord,
			int orientation,
			double pixPerNt,
			double ntAtX,
			double height) {//arrow height
		ArrayList lines = new ArrayList(3);
		
		//account for orientation
		double xS;
		double xE;
		if (orientation == 1) {
			xS = pixPerNt * ((double) start - ntAtX);
			xE = pixPerNt * ((double) end +1 - ntAtX);
		} else {
			xS = pixPerNt * ((double) end +1 - ntAtX);
			xE = pixPerNt * ((double) start - ntAtX);
		}
		//make base line	
		lines.add(new Line2D.Double(xS, yCoord, xE, yCoord));
		//make arrow head
		lines.add(
				new Line2D.Double(
						xE,
						yCoord,
						xE - (height * orientation),
						yCoord + (height * upOrDown)));
		//make but
		lines.add(
				new Line2D.Double(xS, yCoord, xS, yCoord + (height * upOrDown)));				
		return lines;
	}
	
	/**Makes an arrow line with half arrow up for -1 or down for 1.  Start
	 * and stop are really just 5 ' and 3' coordinates, start should always be less
	 * than stop.  Orientation is 1 for + (arrow pointing right), -1 for - (left).
	 * yCoord is for the y-coordinate to draw the main line.  Line thickness is 
	 * needed to calculate whether to draw an arrow or not, for short segments.
	 * Will always return a line of at least 1 pixel.*/
	public static GeneralPath makeArrowLine(
			int start,
			int end,
			int upOrDown,
			float yCoord,
			int orientation,
			float pixPerNt,
			float ntAtX,
			float heightArrowButt,
			float lineThickness) {
		
		GeneralPath gp = new GeneralPath();
		//account for orientation
		float xS;
		float xE;
		if (orientation == 1) {
			xS = pixPerNt * ((float) start - ntAtX);
			xE = pixPerNt * ((float) end+1 - ntAtX);
		} else {
			xS = pixPerNt * ((float) end+1 - ntAtX);
			xE = pixPerNt * ((float) start - ntAtX);
		}	
		
		//account for line thickness
		//x= (0.5*thickness)/(tan22.5)
		float arrowExtend = lineThickness * 1.2071068F;
		float buttExtend = lineThickness * 0.5F;
		float check = xE - xS;
		if (check < 0)
			check *= -1;
		if ((check) > (arrowExtend + buttExtend + 2)) {
			if (orientation == 1) {
				xS += buttExtend;
				xE -= arrowExtend;
			} else {
				xS -= buttExtend;
				xE += arrowExtend;
			}
			//make but
			gp.moveTo(xS, yCoord + (heightArrowButt * upOrDown));
			gp.lineTo(xS, yCoord);
			//make base line
			gp.lineTo(xE, yCoord);
			//make arrow head
			gp.lineTo(
					xE - (heightArrowButt * orientation),
					yCoord + (heightArrowButt * upOrDown));
		} else {
			//check if >= 1 pixel wide
			if (check < 1)
				xE = xS + 1;
			//make single line
			gp.moveTo(xS, yCoord);
			gp.lineTo(xE, yCoord);
		}
		return gp;
	}
	
	/**Lambda, K, H, match, misMatch for DNA and BLASTN parameters*/
	final static String[][] blastParams =
	{
			{
				" 1.10 0.333 0.549 1 -1 ",
				" 0.264 0.0532 0.0722 2 -1 ",
				" ",
				" ",
				" ",
				" ",
				" ",
				" ",
				" ",
			" " },
			{
				" 1.33 0.621 1.12 1 -2 ",
				" 0.549 0.333 0.549 2 -2 ",
				" 0.271 0.130 0.222 3 -2 ",
				" 0.132 0.0532 0.0722 4 -2 ",
				" 0.0514 0.0250 0.0135 5 -2 ",
				" ",
				" ",
				" ",
				" ",
			" " },
			{
				" 1.37 0.711 1.31 1 -3 ",
				" 0.634 0.408 0.912 2 -3 ",
				" 0.366 0.333 0.549 3 -3 ",
				" 0.227 0.161 0.306 4 -3 ",
				" 0.144 0.0952 0.158 5 -3 ",
				" 0.0882 0.0532 0.0722 6 -3 ",
				" 0.0494 0.0291 0.0264 7 -3 ",
				" 0.0212 0.0268 0.00550 8 -3 ",
				" ",
			" " },
			{
				" 1.38 0.738 1.36 1 -4 ",
				" 0.666 0.621 1.12 2 -4 ",
				" 0.410 0.340 0.811 3 -4 ",
				" 0.275 0.333 0.549 4 -4 ",
				" 0.192 0.176 0.357 5 -4 ",
				" 0.136 1.63 0.222 6 -4 ",
				" 0.0957 0.0809 0.132 7 -4 ",
				" 0.0661 0.0532 0.0722 8 -4 ",
				" 0.0435 0.0327 0.0350 9 -4 ",
			" 0.0257 7.93 0.0135 10 -4 " },
			{
				" 1.39 0.747 1.38 1 -5 ",
				" 0.681 0.515 1.24 2 -5 ",
				" 0.432 0.396 0.997 3 -5 ",
				" 0.301 0.306 0.753 4 -5 ",
				" 0.220 0.333 0.549 5 -5 ",
				" 0.164 0.184 0.390 6 -5 ",
				" 0.124 0.140 0.270 7 -5 ",
				" 0.0945 0.103 0.182 8 -5 ",
				" 0.0713 0.0733 0.118 9 -5 ",
			" 0.0529 0.0532 0.0722 10 -5 " },
			{
				" 1.39 0.749 1.38 1 -6 ",
				" 0.687 0.711 1.31 2 -6 ",
				" 0.444 0.621 1.12 3 -6 ",
				" 0.317 1.17 0.912 4 -6 ",
				" 0.237 0.286 0.716 5 -6 ",
				" 0.183 0.333 0.549 6 -6 ",
				" 0.144 0.189 0.414 7 -6 ",
				" 0.114 1.41 0.306 8 -6 ",
				" 0.0904 1.63 0.222 9 -6 ",
			" 0.0718 1.79 0.158 10 -6 " },
			{
				" 1.39 0.750 1.39 1 -7 ",
				" 0.690 0.548 1.34 2 -7 ",
				" 0.451 0.455 1.21 3 -7 ",
				" 0.327 0.386 1.03 4 -7 ",
				" 0.249 0.327 0.854 5 -7 ",
				" 0.196 0.273 0.690 6 -7 ",
				" 0.157 0.333 0.549 7 -7 ",
				" 0.127 0.192 0.431 8 -7 ",
				" 0.104 0.161 0.334 9 -7 ",
			" 0.0855 0.132 0.256 10 -7 " },
			{
				" 1.39 0.750 1.39 1 -8 ",
				" 0.692 0.738 1.36 2 -8 ",
				" 0.455 0.472 1.27 3 -8 ",
				" 0.333 0.621 1.12 4 -8 ",
				" 0.257 0.357 0.965 5 -8 ",
				" 0.205 1.10 0.811 6 -8 ",
				" 0.167 0.263 0.671 7 -8 ",
				" 0.137 0.333 0.549 8 -8 ",
				" 0.114 0.194 0.445 9 -8 ",
			" 0.0958 1.31 0.357 10 -8 " },
			{
				" 1.39 0.750 1.39 1 -9 ",
				" 0.692 0.558 1.37 2 -9 ",
				" 0.458 0.711 1.31 3 -9 ",
				" 0.337 0.427 1.19 4 -9 ",
				" 0.263 0.378 1.05 5 -9 ",
				" 0.211 1.17 0.912 6 -9 ",
				" 0.174 0.295 0.779 7 -9 ",
				" 0.145 0.255 0.657 8 -9 ",
				" 0.122 0.333 0.549 9 -9 ",
			" 0.104 0.196 0.456 10 -9 " },
			{
				" 1.39 0.750 1.39 1 -10 ",
				" 0.693 0.747 1.38 2 -10 ",
				" 0.460 0.491 1.33 3 -10 ",
				" 0.340 1.06 1.24 4 -10 ",
				" 0.267 0.621 1.12 5 -10 ",
				" 0.216 1.03 0.997 6 -10 ",
				" 0.179 0.321 0.871 7 -10 ",
				" 0.151 1.07 0.753 8 -10 ",
				" 0.128 0.250 0.646 9 -10 ",
			" 0.110 0.333 0.549 10 -10 " }
	};
	
	/**Returns Lambda and K or null for a given match and misMatch from BLASTN.
	 Some combinations of match/misMatch are not computed by BLASTN!
	 */
	public static double[] fetchBLASTParams(int match, int misMatch) {
		String x = (blastParams[(misMatch * -1) - 1][match - 1]).trim();
		if (x.equals(""))
			return null;
		String[] nums = x.split("\\s");
		return new double[] {
				Double.parseDouble(nums[0]),
				Double.parseDouble(nums[1])};
	}
	
	public static void saveArrayList(File file, ArrayList al) {
		try {
			ObjectOutputStream out =
				new ObjectOutputStream(new FileOutputStream(file));
			out.writeObject(al);
			out.close();
		} catch (Exception e) {
			System.out.println(
					"\n\nThere appears to be a problem with saving this file: "
					+ file
					+ "\n\n");
			e.printStackTrace();
		}
	}
	
	public static void saveObject(File file, Object ob) {
		try {
			ObjectOutputStream out =
				new ObjectOutputStream(new FileOutputStream(file));
			out.writeObject(ob);
			out.close();
		} catch (Exception e) {
			System.out.println(
					"\n\nThere appears to be a problem with saving this file: "
					+ file
					+ "\n\n");
			e.printStackTrace();
		}
	}
	public static boolean isStringEmpty(String string){
		if (string==null) return true;
		if (string.trim().equals("")) return true;
		return false;
	}
	public static ArrayList fetchArrayList(File file) {
		ArrayList a;
		try {
			ObjectInputStream in =
				new ObjectInputStream(new FileInputStream(file));
			a = (ArrayList) in.readObject();
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
			return new ArrayList();
		}
		return a;
	}
	public static Object fetchObject(File file) {
		Object a;
		try {
			ObjectInputStream in =
				new ObjectInputStream(new FileInputStream(file));
			a = in.readObject();
			in.close();
		} catch (Exception e) {
			//e.printStackTrace();
			a=null;
		}
		return a;
	}	
	
	/**Method to check if a text is null or empty after trimming*/
	public static boolean isEmpty(String text) {
		if (text == null || text.trim().equals(""))	return true;
		return false;
	}
	/**Converts an Object[] to an ArrayList*/	
	public static ArrayList objectArraytoAL(Object[] objs){
		int len = objs.length;
		ArrayList al = new ArrayList(len);
		for (int i=0; i<len; i++) al.add(objs[i]);
		return al;
	}
	/**Deep copies start stops from ExonIntron objects to an ArrayList*/
	public static ArrayList extractStartEnds(ExonIntron[] exons){
		int len = exons.length;
		ArrayList al = new ArrayList(len);
		for (int i=0; i<len; i++){
			al.add(new int[]{exons[i].getStart(),exons[i].getEnd()});
		}
		return al;
	}
	public static int[][] extractStartEndsInts(ExonIntron[] exons){
		int len = exons.length;
		int[][] extra = new int[len][2];
		for (int i=0; i<len; i++){
			extra[i]= new int[]{exons[i].getStart(),exons[i].getEnd()};
		}
		return extra;
	}
	
	public static int convertPlusToNumOrientation(String orientation){
		if (orientation.equals("+")) return 1;
		if (orientation.equals("-")) return -1;
		return 0;
	}
	/**Takes an ArrayList of exon startStop int[2]'s, returns an ArrayList of int[2]'s representing
	 the start,stop for each intron*/
	public static ArrayList fetchIntrons(ArrayList exons) {
		int numIntrons = exons.size() - 1;
		ArrayList introns = new ArrayList();
		for (int i = 0; i < numIntrons; i++) {
			int[] ex1 = (int[]) exons.get(i);
			int[] ex2 = (int[]) exons.get(i + 1);
			introns.add(new int[] { ex1[1] + 1, ex2[0] - 1 });
		}
		return introns;
	}
	
	public static String ALToString(ArrayList al) {
		StringBuffer sb = new StringBuffer();
		int len = al.size();
		for (int i = 0; i < len; i++) {
			int[] ss = (int[]) al.get(i);
			sb.append("\n  " + ss[0] + "-" + ss[1]);
		}
		return sb.toString();
	}
	
	/**Sorts an ArrayList of int[2]'s representing start, stop*/
	public static ArrayList sortSegments(ArrayList segs){
		int num = segs.size();
		ArrayList sorted = new ArrayList (num);
		Region[] se = new Region[num];
		for (int i=0; i<num; i++){
			int[] stst = (int[])segs.get(i);
			se[i] = new Region(stst[0], stst[1]);
		}
		Arrays.sort(se);
		for (int i=0; i<num; i++){
			sorted.add(se[i].getStartStop());
		}
		return sorted;
	}
	
	/**Takes an ArrayList of int[2]'s representing start stop coords
	 * kills the dups, and those cotained within another seg, will
	 * merge overlapping segments.
	 */
	public static ArrayList joinSegments(ArrayList segs) {
		ArrayList mergedSegs = new ArrayList();
		int len = segs.size();
		int[] seg1 = new int[2];
		while (len > 0) {
			boolean test = true;
			seg1 = (int[]) segs.remove(0);
			len--;
			for (int i = 0; i < len; i++) {
				int[] seg2 = (int[]) segs.get(i);
				//if it entirely overlaps, contained within
				if (seg1[0] >= seg2[0] && seg1[1] <= seg2[1]) {
					test = false;
					break;
				} //don't save
				else if (seg1[0] <= seg2[0] && seg1[1] >= seg2[1]) {
					segs.remove(i);
					len--;
					i--;
					continue;
				}
				//if immediately adjacent fuse
				if (seg1[1] == (seg2[0] - 1)) {
					seg1 = new int[] { seg1[0], seg2[1] };
					segs.remove(i);
					len--;
					i = -1; //start over
					continue;
				} else if (seg1[0] == (seg2[1] + 1)) {
					seg1 = new int[] { seg2[0], seg1[1] };
					segs.remove(i);
					len--;
					i = -1; //start over
					continue;
				}
				//if it doesn't overlap then continue
				if (seg1[1] < seg2[0] || seg1[0] > seg2[1]) {
					continue;
				}
				//if it partially overlaps fuse together
				if (seg1[1] >= seg2[0] && seg1[0] <= seg2[0]) {
					seg1 = new int[] { seg1[0], seg2[1] };
					segs.remove(i);
					len--;
					i = -1; //start over
				} else if (seg1[0] <= seg2[1] && seg1[1] >= seg2[1]) {
					seg1 = new int[] { seg2[0], seg1[1] };
					segs.remove(i);
					len--;
					i = -1; //start over 
				}
			}
			//reached stop add exon
			if (test) {
				mergedSegs.add(seg1);
			}
		}
		return mergedSegs;
	}	
	//stop
}
