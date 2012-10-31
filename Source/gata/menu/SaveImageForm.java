package gata.menu;

import gata.main.*;

import javax.swing.*;
import java.awt.*;

/**
 * @author Nix
 * Frame to hold panel for collecting how the user wants the image saved.
 * */
public class SaveImageForm extends JFrame {
	//fields
	SaveImagePanel panel;
	
	public SaveImageForm(GATAParams params) {
		//modify frame
		setTitle("GATAPlot Image Parameters");
		
		//add panel
		panel = new SaveImagePanel (params,this);
		Container contentPane = getContentPane();
		contentPane.add(panel);

		//set size of frame from panel info
		setBounds(100,100,panel.getMaxWidth(),panel.getMaxHeight());
		show();
	}
}
