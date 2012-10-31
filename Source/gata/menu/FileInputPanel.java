package gata.menu;

import gata.aligner.*;
import gata.main.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;

/**
 * @author nix
 * Panel used to fetch user inputs.
 */
public class FileInputPanel extends JPanel {
	//fields
	private int x = 10;
	private int y = 10;
	private JTextField alignerFile;
	private JTextField gffFileRef;
	private JTextField gffFileComp;
	private JFrame frame;
	private JFileChooser chooser;
	private int maxWidth;
	private int maxHeight;
	private GATAParams params;
	
	public FileInputPanel(GATAParams params, JFrame frame) {
		//modify panel 
		setBackground(Color.WHITE);
		setLayout(null);
		chooser = new JFileChooser();
		this.frame = frame;
		this.params = params;
		
		//get Aligner saved folder
		File objectsFile = params.getObjectsFile();
		String fileName = "";
		try{
			fileName = objectsFile.getCanonicalPath()+".gata";
		}
		catch(IOException e){
			e.printStackTrace();
		}
		chooser.setCurrentDirectory(objectsFile);
		//aligner file
		JLabel cap = GATAUtil.makeLabel(this,"Enter a GATAligner Alignment file (required).",14,0);
		cap.setForeground(Color.BLUE);
		int[] capSize = GATAUtil.setLabel(cap, x, y);
		y += capSize[1]+5;
		
		alignerFile = GATAUtil.makeFieldSimple(this,fileName);
		int[] alignerFileSize = GATAUtil.setField(alignerFile, x, y, 400);
		x+= alignerFileSize[0]+5;
		JButton alignerBrows = GATAUtil.makeButton("Browse", 12, new BrowseButtonAction(alignerFile,frame,chooser),this);
		int[] alignerBrowsSize = GATAUtil.setButton(alignerBrows, x,y);
		maxWidth = alignerBrowsSize[0]+10+x;
		y+=alignerBrowsSize[1]+20;
		x= 10;
		
		//gff file RefSeq
		JLabel gff = GATAUtil.makeLabel(this,"(Optional) Enter an annotation GFF file for the Reference Sequence.",14,0);
		gff.setForeground(Color.BLUE);
		int[] gffSize = GATAUtil.setLabel(gff, x, y);
		y += gffSize[1]+5;
		
		gffFileRef = GATAUtil.makeFieldSimple(this,"");
		int[] gffFileSize = GATAUtil.setField(gffFileRef, x, y, 400);
		x+= gffFileSize[0]+5;
		JButton gffBrows = GATAUtil.makeButton("Browse", 12, new BrowseButtonAction(gffFileRef,frame,chooser),this);
		int[] gffBrowsSize = GATAUtil.setButton(gffBrows, x,y);	
		y+=gffBrowsSize[1]+20;
		x=10;
		//gff file CompSeq
		JLabel gff2 = GATAUtil.makeLabel(this,"(Optional) Enter an annotation GFF file for the Comparative Sequence.",14,0);
		gff2.setForeground(Color.BLUE);
		int[] gffSize2 = GATAUtil.setLabel(gff2, x, y);
		y += gffSize2[1]+5;
		
		gffFileComp = GATAUtil.makeFieldSimple(this,"");
		int[] gffFileSizeComp = GATAUtil.setField(gffFileComp, x, y, 400);
		x+= gffFileSize[0]+5;
		JButton gffBrows2 = GATAUtil.makeButton("Browse", 12, new BrowseButtonAction(gffFileComp,frame,chooser),this);
		int[] gffBrowsSize2 = GATAUtil.setButton(gffBrows2, x,y);	
		y+=gffBrowsSize2[1]+20;
				
		
		//seperator
		JSeparator sep = new JSeparator();
		sep.setBounds(0,y,maxWidth,5);
		add(sep);
		y+=10;
		
		//Go
		JButton go = GATAUtil.makeButton(" Go ", 12, new GoButtonAction(),this);
		Dimension goDim = go.getPreferredSize();
		JButton cancel = GATAUtil.makeButton(" Cancel ", 12, new CancelButtonAction(),this);
		Dimension cancelDim = cancel.getPreferredSize();
		int wid = (int)(goDim.getWidth()+5+cancelDim.getWidth());
		x = maxWidth/2 - wid/2;
		int[] goSize = GATAUtil.setButton(go,x,y);
		x+= goSize[0]+5;
		GATAUtil.setButton(cancel,x,y);
		maxHeight = goSize[1]+y+35;

	}
	public int getMaxHeight() {
		return maxHeight;
	}
	public int getMaxWidth() {
		return maxWidth;
	}
	private class GoButtonAction implements ActionListener {
		
		public void actionPerformed(ActionEvent e) {
			//check for GATAlignment file
			String alignerFileString = alignerFile.getText();
			if (alignerFileString.equals("")){
				GATAUtil.throwWarning(frame, alignerFile,"No GATAligner alignment file?!");
				return;
			}
			if (GATAUtil.checkFile(alignerFileString)==false) return;
			File file = new File(alignerFileString);
			params.setObjectsFile(file);
			
			//check for gff files
			String gffFileString = gffFileRef.getText();
			if (gffFileString.equals("")==false){
				if (GATAUtil.checkFile(gffFileString)==false) return;
				file = new File(gffFileString);
				AnnoSpecParams asp = new AnnoSpecParams();
				asp.setGffFile(file);
				params.setGffRefPresent(true);
				params.setRefAnnoSpecParams(asp);
			}
			gffFileString = gffFileComp.getText();
			if(gffFileString.equals("")==false){
				if (GATAUtil.checkFile(gffFileString)==false) return;
				file = new File(gffFileString);
				AnnoSpecParams asp = new AnnoSpecParams();
				asp.setGffFile(file);
				params.setGffCompPresent(true);
				params.setCompAnnoSpecParams(asp);
			}
			frame.dispose();
			params.getGataPlotter().launchGATAPlotter();	
		}
	}
	private class CancelButtonAction implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			frame.dispose();
		}
	}
}
