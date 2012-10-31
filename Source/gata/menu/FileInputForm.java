package gata.menu;

import gata.main.*;

import javax.swing.*;
import java.awt.*;

/**
 * @author Nix
 * Container to hold the Plotter Input Panel, place for user inputs.
 * */

public class FileInputForm extends JFrame {
	//fields
	FileInputPanel panel;
	
	public FileInputForm(GATAParams params) {
		//modify frame
		setResizable(false);
		setTitle("GATAPlotter Input Files");
		
		//add panel
		panel = new FileInputPanel(params,this);
		Container contentPane = getContentPane();
		contentPane.add(panel);

		//set size of frame from panel info
		setBounds(100,100,panel.getMaxWidth(),panel.getMaxHeight());
		show();
	}
}
