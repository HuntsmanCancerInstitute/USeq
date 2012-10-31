package gata.aligner;

import gata.main.*;

import javax.swing.*;
import java.awt.event.*;

/**
 * @author nix
 * Focus listener firing off upon loss of focus of a text field to check whether the entered file exists
 * */

public class FileFocusListener implements FocusListener {
		JTextField field;
		JFrame frame;
		JFileChooser chooser;
		public void focusGained(FocusEvent fe) {
		}
		public void focusLost(FocusEvent fe) {
			//check if anything was entered
			if (field.getText().equals("") == false) {
				//check if a real file or dir can be found
				if (GATAUtil.checkFile(field.getText()) == false) {
					//throw warning
					Object[] options = { "OK", "Browse" };
					int selected =
						JOptionPane.showOptionDialog(
							frame,
							"Sorry, I can't find your file, program, or requested directory!",
							null,
							JOptionPane.DEFAULT_OPTION,
							JOptionPane.WARNING_MESSAGE,
							null,
							options,
							options[0]);
					if (selected == 1)
						GATAUtil.showFileChooserDialogBox(field,frame,chooser);
					//set focus to this field
					else
						field.requestFocusInWindow();
				}
			}
		}

		public void setRefs(JTextField x, JFrame f, JFileChooser c) {
			field = x;
			frame = f;
			chooser = c;
		}
	}