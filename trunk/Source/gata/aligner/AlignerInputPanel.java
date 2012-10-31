package gata.aligner;

import gata.main.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;

/**
 * @author nix
 * A swing panel for collecting user input needed to run GATAligner*/
public class AlignerInputPanel extends JPanel {
	//	fields
	private int A = 10;
	private int B = 20;
	private int C = 10;
	private int D = 25;
	private int E = 35;
	private int F = 0;
	private int G = 10;
	private int H = 40;
	private int I = 5;
	private int J = 30;
	private int K = 20;
	private int L = 50;
	private int sizeTextFields = 350;
	private int xRight;
	private int furthestRight;
	private int furthestDown;
	private String finalBaseName;
	private double[] stats;

	//object references
	private JTextField refF;
	private JButton refBrows;
	private JTextField compF;
	private JButton compBrows;
	private JTextField bl2F;
	private JButton bl2FBrows;
	private JTextField objF;
	private JButton objBrows;
	private JTextField baseNameF;
	private JButton baseNameBrows;

	private JSpinner matSpin;
	private JSpinner misSpin;
	private JSpinner createSpin;
	private JSpinner extSpin;
	private JRadioButton maskYes;
	private JRadioButton maskNo;
	private JRadioButton extractYes;
	private JRadioButton extractNo;
	
	private JTextField winTexF;
	private JTextField scoreTexF;
	private JTextField bitTexF;
	private JTextField refTexF;
	private JTextField compTexF;

	private JButton help;
	private JButton reset;
	private JButton go;

	private JFileChooser chooser;
	private JFileChooser chooserDirectories;
	private JFrame frame;
	private AlignerPreferences ap;
	private GATAligner gata;

	/**Complex class containing most everything to show and process the GATAligner inputs.*/
	public AlignerInputPanel(JFrame frame, GATAligner gata) {
		this.gata = gata;
		this.frame = frame;
		//modify panel 
		setBackground(Color.WHITE);
		setLayout(null);

		//look for AlignerPreferences object, will return a new one if none found on disk
		ap = getAlignerPreferences();
		
		//look to see if bl2seq has been specified
		if (ap.getBl2SeqProg().equals("")){
			//attempt to set to one of the included versions
			String os = System.getProperty("os.name").toLowerCase();
			File file;
			try{
				if (os.matches(".*windows.*")){
					file = new File("./Bl2Seqs", "Windows_bl2seq.exe");
					ap.setBl2SeqProg(file.getCanonicalPath());
				} 
				else if (os.matches(".*linux.*")){
					file = new File("./Bl2Seqs", "Linux_bl2seq");
					ap.setBl2SeqProg(file.getCanonicalPath());
				} 
				else if (os.matches(".*mac.*")){
					file = new File("./Bl2Seqs", "MacOSX_bl2seq");
					ap.setBl2SeqProg(file.getCanonicalPath());
				} 
			//else if (os.matches)  
			} catch (IOException e){
			e.printStackTrace();
			}   	
		}
		
		//make a JFileChooserObject for browsing, buggie, custom showDialog doesn't allow folder selection! 
		chooser = new JFileChooser(); 
		chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
		
		//make widgets
		//create caption
		JLabel cap =
			GATAUtil.makeLabel(this,
				"Graphic Alignment Tool For Comparative Sequence AnalysisBean - GATAligner (1.0)",
				16,
				0);

		int[] capSize = GATAUtil.setLabel(cap, A, A);
		//getting back width and height of label
		int yLeft = capSize[1] + B + A;
		int yRight = capSize[1] + B + A;

		//create Locations
		JLabel loc =
		GATAUtil.makeLabel(this,"Enter the path and text for the following:", 14, 0);
		loc.setForeground(Color.BLUE);
		int[] locSize = GATAUtil.setLabel(loc, A, yLeft);
		yLeft += locSize[1] + C;

		//Reference Sequence
		refF = GATAUtil.makeFieldSimple(this,ap.getRefSeq());

		int[] refFSize = GATAUtil.setField(refF, D, yLeft, sizeTextFields);

		int xLeft = refFSize[0] + D + I;
		refBrows = GATAUtil.makeButton("Browse", 10, new BrowseButtonAction(refF,frame,chooser),this);
		int[] refBrowsSize = GATAUtil.setButton(refBrows, xLeft, yLeft);
		yLeft += refFSize[1];
		//set xRight
		xRight = xLeft + refBrowsSize[0] + J;

		JLabel ref = GATAUtil.makeLabel(this,"Single-Fasta Reference Sequence", 10, 1);
		int[] refSize = GATAUtil.setLabel(ref, E, yLeft);
		yLeft += refSize[1] + G;

		//Comparative Sequence
		compF = GATAUtil.makeFieldSimple(this, ap.getCompSeq());
		int[] compFSize = GATAUtil.setField(compF, D, yLeft, sizeTextFields);
		compBrows = GATAUtil.makeButton("Browse", 10, new BrowseButtonAction(compF,frame,chooser),this);
		GATAUtil.setButton(compBrows, xLeft, yLeft);
		yLeft += compFSize[1];
		JLabel comp = GATAUtil.makeLabel(this,"Multi-Fasta Comparative Sequence(s)", 10, 1);
		int[] compSize = GATAUtil.setLabel(comp, E, yLeft);
		yLeft += compSize[1] + G;

		//Blast2Seq
		bl2F = GATAUtil.makeFieldSimple(this, ap.getBl2SeqProg());
		int[] bl2FSize = GATAUtil.setField(bl2F, D, yLeft, sizeTextFields);
		bl2FBrows = GATAUtil.makeButton("Browse", 10, new BrowseButtonAction(bl2F,frame,chooser),this);
		GATAUtil.setButton(bl2FBrows, xLeft, yLeft);
		yLeft += bl2FSize[1];
		JLabel bl2 = GATAUtil.makeLabel(this,"bl2seq BLAST Program", 10, 1);
		int[] bl2Size = GATAUtil.setLabel(bl2, E, yLeft);
		yLeft += bl2Size[1] + G;

		//Alignement Object Folder
		objF = GATAUtil.makeFieldSimple(this, ap.getStorageLoc());
		int[] objFSize = GATAUtil.setField(objF, D, yLeft, sizeTextFields);
		objBrows = GATAUtil.makeButton("Browse", 10, new BrowseButtonAction(objF,frame,chooser),this);
		GATAUtil.setButton(objBrows, xLeft, yLeft);
		yLeft += objFSize[1];
		JLabel obj =
		GATAUtil.makeLabel(this,
				"Select a folder for saving the GATA alignment",
				10,
				1);
		int[] objSize = GATAUtil.setLabel(obj, E, yLeft);
		yLeft += objSize[1] + G;

		//Alignement Object BaseName
		baseNameF = new JTextField(ap.getBaseName());
		add(baseNameF);
		int[] baseNameFSize = GATAUtil.setField(baseNameF, D, yLeft, sizeTextFields);
		yLeft += baseNameFSize[1];
		JLabel baseName =
		GATAUtil.makeLabel(this,"Enter a file base text for the GATA alignemnt", 10, 1);
		int[] baseNameSize = GATAUtil.setLabel(baseName, E, yLeft);
		yLeft += baseNameSize[1];

		//Column 2
		//BlastN Params

		JLabel bla = GATAUtil.makeLabel(this,"BLASTN Parameters", 14, 0);
		bla.setForeground(Color.BLUE);
		int[] blaSize = GATAUtil.setLabel(bla, xRight, yRight);
		yRight += blaSize[1] + C;

		//Match
		JLabel matL = GATAUtil.makeLabel(this,"Nt:  Match (1 to 8)", 12, 0);
		int xAdj = xRight + (D - A);
		int[] matLSize = GATAUtil.setLabel(matL, xAdj, yRight);
		xAdj += matLSize[0] + I;
		matSpin = new JSpinner(new SpinnerNumberModel(5, 1, 8, 1));
		matSpin.setValue(ap.getMatch());
		int[] matSize = setSpinner(matSpin, xAdj, yRight);
		matSpin.setBounds(xAdj, yRight, matSize[0] + 4, matSize[1]);
		//needed to avoid cramping
		xAdj += matSize[0] + K;
		//MisMatch
		JLabel misMatch = GATAUtil.makeLabel(this,"MisMatch (-1 to -10)", 12, 0);
		int[] misSize = GATAUtil.setLabel(misMatch, xAdj, yRight);
		xAdj += misSize[0] + I;
		misSpin = new JSpinner(new SpinnerNumberModel(-4, -10, -1, 1));
		misSpin.setValue(ap.getMisMatch());
		setSpinner(misSpin, xAdj, yRight);

		yRight += matSize[1] + G;
		xAdj = xRight + (D - A);

		//Gap Creation
		JLabel create = GATAUtil.makeLabel(this,"Gap:  Creation (0 to -50)", 12, 0);
		int[] createSize = GATAUtil.setLabel(create, xAdj, yRight);
		xAdj += createSize[0] + I;
		createSpin = new JSpinner(new SpinnerNumberModel(-10, -50, 0, 1));
		createSpin.setValue(ap.getCreate());
		int[] createSpinSize = setSpinner(createSpin, xAdj, yRight);
		xAdj += createSpinSize[0] + K;
		//Extension
		JLabel ext = GATAUtil.makeLabel(this,"Extension (0 to -50)", 12, 0);
		int[] extSize = GATAUtil.setLabel(ext, xAdj, yRight);
		xAdj += extSize[0] + I;
		extSpin = new JSpinner(new SpinnerNumberModel(-4, -50, 0, 1));
		extSpin.setValue(ap.getExtend());
		int[] extSpinSize = setSpinner(extSpin, xAdj, yRight);
		//set right most boundary
		furthestRight = xAdj + extSpinSize[0] + A;
		yRight += createSpinSize[1] + G;
		xAdj = xRight + (D - A);

		//Mask low complexity regions
		JLabel mask = GATAUtil.makeLabel(this,"Mask low complexity regions?", 12, 0);
		int[] maskSize = GATAUtil.setLabel(mask, xAdj, yRight);
		xAdj += maskSize[0] + I;
		ButtonGroup maskGroup = new ButtonGroup();
		maskYes = new JRadioButton("Yes", false);
		maskYes.setBackground(Color.WHITE);
		maskGroup.add(maskYes);
		maskNo = new JRadioButton("No", true);
		maskNo.setBackground(Color.WHITE);
		maskGroup.add(maskNo);
		Dimension dimYes = maskYes.getPreferredSize();
		Dimension dimNo = maskNo.getPreferredSize();
		maskYes.setBounds(
			xAdj,
			yRight,
			(int) dimYes.getWidth(),
			(int) dimYes.getHeight());
		add(maskYes);
		xAdj += dimYes.getWidth() + I;
		maskNo.setBounds(
			xAdj,
			yRight,
			(int) dimNo.getWidth(),
			(int) dimNo.getHeight());
		add(maskNo);
		if (ap.getMask().equals("Yes"))
			maskYes.setSelected(true);
		yRight += dimYes.getHeight() + B;

		//GataAligner Parameters
		JLabel gap = GATAUtil.makeLabel(this,"GATAligner Parameters", 14, 0);
		gap.setForeground(Color.BLUE);
		int[] gapSize = GATAUtil.setLabel(gap, xRight, yRight);
		yRight += gapSize[1] + C;
		xAdj = xRight + (D - A);
		//window size
		JLabel win = GATAUtil.makeLabel(this,"Window Size (7 to 100)", 12, 0);
		int[] winSize = GATAUtil.setLabel(win, xAdj, yRight);
		xAdj += winSize[0] + I;
		winTexF = makeField("24", new NumberFocusListener(24, 7, 100));
		winTexF.setText(ap.getWindow());
		int[] winTexFSize = GATAUtil.setField(winTexF, xAdj, yRight, 40);
		xAdj += winTexFSize[0] + K;
		yRight += G + winTexFSize[1];
		xAdj = xRight + (D - A);
		//score cut off
		JLabel score = GATAUtil.makeLabel(this,"Lower Cut Off Score", 12, 0);
		int[] scoreSize = GATAUtil.setLabel(score, xAdj, yRight);
		xAdj += scoreSize[0] + I;
		RawScoreFocusListener rawListener =
			new RawScoreFocusListener(70, 1, 10000);
		scoreTexF = makeField("70", rawListener);
		scoreTexF.setText(ap.getScore());
		int[] scoreTexFSize = GATAUtil.setField(scoreTexF, xAdj, yRight, 40);
		xAdj += scoreTexFSize[0];
		//raw label
		JLabel raw = GATAUtil.makeLabel(this,"(raw) ", 12, 0);
		int[] rawSize = GATAUtil.setLabel(raw, xAdj, yRight);
		xAdj += rawSize[0] + I + I;
		//bit field
		BitScoreFocusListener bitListener =
			new BitScoreFocusListener(21.9, 1, 1000);
		bitListener.setRawField(scoreTexF);
		bitTexF = makeField("21.9", bitListener);
		rawListener.setBitField(bitTexF);
		stats =
			GATAUtil.fetchBLASTParams(
				((Integer) matSpin.getValue()).intValue(),
				((Integer) misSpin.getValue()).intValue());
		bitTexF.setText(
			""
				+ GATAUtil.convertRawScoreToBit(
					stats[0],
					stats[1],
					Double.parseDouble(scoreTexF.getText())));
		int[] bitTexFSize = GATAUtil.setField(bitTexF, xAdj, yRight, 40);
		xAdj += bitTexFSize[0];
		JLabel bit = GATAUtil.makeLabel(this,"(bit) ", 12, 0);
		GATAUtil.setLabel(bit, xAdj, yRight);
		yRight += winTexFSize[1] + G;
		
		//start index reference sequence
		xAdj = xRight + (D - A);
		JLabel refS = GATAUtil.makeLabel(this,"Start Position Reference Sequence", 12, 0);
		int[] refSSize = GATAUtil.setLabel(refS, xAdj, yRight);
		xAdj += refSSize[0] + I;
		refTexF = makeField("1", new NumberFocusListener(1, 1, 1000000000));
		refTexF.setText(ap.getRefStart());
		int[] refTexFSize = GATAUtil.setField(refTexF, xAdj, yRight, 70);
		yRight += refTexFSize[1] + G;
		//start index comp seq
		xAdj = xRight + (D - A);
		JLabel compS = GATAUtil.makeLabel(this,"Start Position Comparative Sequence", 12, 0);
		int[] compSSize = GATAUtil.setLabel(compS, xAdj, yRight);
		xAdj += compSSize[0] + I;
		compTexF = makeField("1", new NumberFocusListener(1, 1, 1000000000));
		compTexF.setText(ap.getCompStart());
		int[] compTexFSize = GATAUtil.setField(compTexF, xAdj, yRight, 70);
		yRight += compTexFSize[1] + G;
		
		//Extract start index positions from fasta header option, defaults to yes
		xAdj = xRight + (D - A);
		JLabel extract = GATAUtil.makeLabel(this,"Extract Starts from FASTA '>' line?", 12, 0);
		int[] extractSize = GATAUtil.setLabel(extract, xAdj, yRight);
		xAdj += extractSize[0] + I;
		ButtonGroup extractGroup = new ButtonGroup();
		extractYes = new JRadioButton("Yes", false);
		extractYes.setBackground(Color.WHITE);
		extractGroup.add(extractYes);
		extractNo = new JRadioButton("No", true);
		extractNo.setBackground(Color.WHITE);
		extractGroup.add(extractNo);
		Dimension dimExtractYes = extractYes.getPreferredSize();
		Dimension dimExtractNo = extractNo.getPreferredSize();
		extractYes.setBounds(
			xAdj,
			yRight-5,
			(int) dimExtractYes.getWidth(),
			(int) dimExtractYes.getHeight());
		add(extractYes);
		xAdj += dimExtractYes.getWidth() + I;
		extractNo.setBounds(
			xAdj,
			yRight-5,
			(int) dimExtractNo.getWidth(),
			(int) dimExtractNo.getHeight());
		add(extractNo);
		if (ap.getExtract().equals("Yes"))
			extractYes.setSelected(true);
		yRight += dimExtractYes.getHeight() + B;
		
		//set furthest down
		furthestDown = yRight;
		
				
		

		//buttons
		yLeft += B;
		xAdj = L;
		help = GATAUtil.makeButton("Help", 10, new HelpButtonAction(),this);
		int[] helpSize = GATAUtil.setButton(help, xAdj, yLeft);
		xAdj += helpSize[0] + I;

		reset = GATAUtil.makeButton("Reset Form", 10, new ResetButtonAction(),this);
		int[] resetSize = GATAUtil.setButton(reset, xAdj, yLeft);
		xAdj += resetSize[0] + I + I + I + I;

		go = GATAUtil.makeButton("Align Sequences", 12, new GoButtonAction(),this);
		go.setForeground(new Color(99, 0, 0));
		int[] goSize = GATAUtil.setButton(go, xAdj, yLeft);
		
	}

	private class BitScoreFocusListener extends NumberFocusListener {
		JTextField rawField;
		public BitScoreFocusListener(double def, int min, int max) {
			super(def, min, max);
		}

		public void focusLost(FocusEvent fe) {
			//check number
			if (checkIfNumber() == false)
				return;
			//check if it's within bounds
			if (checkBounds() == false)
				return;
			//update raw score field
			stats =
				GATAUtil.fetchBLASTParams(
					((Integer) matSpin.getValue()).intValue(),
					((Integer) misSpin.getValue()).intValue());
			if (stats == null) {
				rawField.setText("0");
				field.setText("0");
				return;
			}
			rawField.setText(
				""
					+ GATAUtil.convertBitScoreToRaw(
						stats[0],
						stats[1],
						Double.parseDouble(field.getText())));

		}
		public void setRawField(JTextField x) {
			rawField = x;
		}

	}

	
	private class GoButtonAction implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			//check if each of the file fields are filled in
			if (refF.getText().equals(""))
				GATAUtil.throwWarning(frame,refF, "Sorry, GATAligner needs a reference sequence!");
			else if (compF.getText().equals(""))
			GATAUtil.throwWarning(frame,compF, "Sorry, GATAligner needs a comparative sequence!");
			else if (bl2F.getText().equals(""))
			GATAUtil.throwWarning(frame,bl2F, "No bl2seq program?!");
			else if (objF.getText().equals(""))
			GATAUtil.throwWarning(frame,objF, "Where do you want to save the GATAlignemnt?");
			else {
				if (checkFields()) {
					String maskOpt = "No";
					if (maskYes.isSelected())
						maskOpt = "Yes";
					String extractOpt = "No";
					if (extractYes.isSelected()){
						extractOpt = "Yes";
						compTexF.setText("1");
						refTexF.setText("1");
					}
						
					ap.setPreferences(
						refF.getText(),
						compF.getText(),
						bl2F.getText(),
						objF.getText(),
						finalBaseName,
						(Integer) matSpin.getValue(),
						(Integer) misSpin.getValue(),
						(Integer) createSpin.getValue(),
						(Integer) extSpin.getValue(),
						maskOpt,
						extractOpt,
						winTexF.getText(),
						scoreTexF.getText(),
						refTexF.getText(),
						compTexF.getText());
					gata.lauchAligner();	
					ap.saveAlignerPreferences();
				}
			}
		}
	}

	private class HelpButtonAction implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			TextWindow t =
				new TextWindow(
					"GATAliner Help",
					600,
					600,
					"file:./documents/AlignerHelp.html");
		}
	}

	private class NumberFocusListener implements FocusListener {
		double def;
		JTextField field;
		int max;
		int min;
		double num;
		public NumberFocusListener(double def, int min, int max) {
			this.def = def;
			this.min = min;
			this.max = max;
		}
		public boolean checkBounds() {
			if ((num >= min && num <= max) == false) {
				//throw warning
				JOptionPane.showMessageDialog(
					frame,
					"This number is out of range!",
					null,
					JOptionPane.WARNING_MESSAGE);
				field.setText(def + "");
				//set focus to this field
				field.requestFocusInWindow();
				return false;
			}
			return true;
		}

		public boolean checkIfNumber() {
			try {
				num = Double.parseDouble(field.getText());
			} catch (NumberFormatException e) {
				JOptionPane.showMessageDialog(
					frame,
					"That is not a number!",
					null,
					JOptionPane.WARNING_MESSAGE);
				field.requestFocusInWindow();
				return false;
			}
			return true;
		}
		public void focusGained(FocusEvent f) {
		}
		public void focusLost(FocusEvent fe) {
			//check number
			if (checkIfNumber() == false)
				return;
			//check if it's within bounds
			if (checkBounds() == false)
				return;

		}

		public void setRef(JTextField x) {
			field = x;
		}
	}
	private class RawScoreFocusListener extends NumberFocusListener {
		JTextField bitField;
		public RawScoreFocusListener(double def, int min, int max) {
			super(def, min, max);
		}

		public void focusLost(FocusEvent fe) {
			//check number
			if (checkIfNumber() == false)
				return;
			//check if it's within bounds
			if (checkBounds() == false)
				return;
			//update bit score field
			stats =
				GATAUtil.fetchBLASTParams(
					((Integer) matSpin.getValue()).intValue(),
					((Integer) misSpin.getValue()).intValue());
			if (stats == null) {
				bitField.setText("0");
				field.setText("0");
				return;
			}
			bitField.setText(
				""
					+ GATAUtil.convertRawScoreToBit(
						stats[0],
						stats[1],
						Double.parseDouble(field.getText())));

		}
		public void setBitField(JTextField x) {
			bitField = x;
		}

	}
	private class ResetButtonAction implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			resetForm();
		}
	}

	public boolean checkFields() {
		//check bl2seq
		if (bl2F.getText().matches(".*bl2seq.*") == false) {
			GATAUtil.throwWarning(frame,bl2F, "Where is the bl2seq program?");
			return false;
		}
		if (GATAUtil.checkFile(bl2F.getText())==false) return false;
		//check results directory
		File f = new File(objF.getText());
		if ((f.isDirectory() && f.canWrite()) == false) {
			GATAUtil.throwWarning(frame,objF, "Cannot write to/find your selected results folder!");
			return false;
		}
		//check seq files
		f = new File(refF.getText());
		if ((f.isFile() && f.canRead()) == false) {
			GATAUtil.throwWarning(frame,refF, "Cannot read/find your selected Reference sequence file!");
			return false;
		}
		f = new File(compF.getText());
		if ((f.isFile() && f.canRead()) == false) {
			GATAUtil.throwWarning(frame,
				compF,
				"Cannot read/find your selected Comparative sequence file!");
			return false;
		}
		//clean up and check base text
		finalBaseName = baseNameF.getText().trim();
		if (finalBaseName.equals("")) {
			GATAUtil.throwWarning(frame,baseNameF, "Choose an appropriate base text!");
			return false;
		}
		//check cut off score
		double topScore =
			Double.parseDouble(winTexF.getText())
				* (double) ((Integer) (matSpin.getValue())).intValue();
		if (topScore < Double.parseDouble(scoreTexF.getText())) {
			GATAUtil.throwWarning(frame,
				bitTexF,
				"Your Cut Off Score exceeds the maximum score obtainable for the given Window Size\n  Lower the Cut Off Score.");
			return false;
		}
		stats =
			GATAUtil.fetchBLASTParams(
				((Integer) matSpin.getValue()).intValue(),
				((Integer) misSpin.getValue()).intValue());
		if (stats == null) {
			GATAUtil.throwWarning(frame,
				bitTexF,
				"Your Match and MisMatch parameters are incompatible with the bl2seq/BLASTN program.\n  Please change them (try 5/-4 or 1/-3).");
			return false;
		}
		return true;
	}

	public AlignerPreferences getAlignerPreferences() {
		AlignerPreferences ap;
		try {
			ObjectInputStream in =
				new ObjectInputStream(
					new FileInputStream("GATAlignerPreferences"));
			ap = (AlignerPreferences) in.readObject();
			in.close();
		} catch (Exception e) {
			System.out.println(
				"Cannot find/open GATAlignerPreferences, making a new one....");
			return new AlignerPreferences();
		}
		return ap;
	}
	public int getFurthestDown() {
		return furthestDown;
	}
	public int getFurthestRight() {
		return furthestRight;
	}
	
	public JTextField makeField(String text, NumberFocusListener focus) {
		JTextField tex = new JTextField(text);
		tex.addFocusListener(focus);
		focus.setRef(tex);
		add(tex);
		return tex;
	}

	public void resetForm() {
		//JTextFields
		refF.setText("");
		compF.setText("");
		objF.setText("");
		baseNameF.setText("");
		winTexF.setText("24");
		scoreTexF.setText("70");
		bitTexF.setText("21.9");
		refTexF.setText("1");
		compTexF.setText("1");
		//spinners
		matSpin.setValue(new Integer(5));
		misSpin.setValue(new Integer(-4));
		createSpin.setValue(new Integer(-10));
		extSpin.setValue(new Integer(-4));
		//radio buttons
		maskNo.setSelected(true);
		extractNo.setSelected(true);
	}
	
	
	public int[] setSpinner(JSpinner s, int x, int y) {
		Dimension dim = s.getPreferredSize();
		int[] wh = {(int) dim.getWidth(), (int) dim.getHeight()};
		s.setBounds(x, y, wh[0], wh[1]);
		add(s);
		return wh;
	}

	
	public AlignerPreferences getAlignerPrefReference() {
		return ap;
	}
	public JTextField getBl2F() {
		return bl2F;
	}
	public JTextField getCompTexF() {
		return compTexF;
	}
	public JTextField getRefTexF() {
		return refTexF;
	}
}