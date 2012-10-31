package gata.aligner;

import gata.main.*;

import javax.swing.*;
import java.awt.event.*;


/**
 * @author nix
 * Browse button action for AlignerInputPanel.  Why isn't this an inner class?*/
public class BrowseButtonAction implements ActionListener {
		JTextField fieldRef;
		JFrame frame;
		JFileChooser chooser;
		public BrowseButtonAction(JTextField fieldRef,JFrame frame,JFileChooser chooser) {
			this.fieldRef = fieldRef;
			this.frame = frame;
			this.chooser = chooser;
		}
		public void actionPerformed(ActionEvent e) {
			GATAUtil.showFileChooserDialogBox(fieldRef,frame, chooser);
		}
	}