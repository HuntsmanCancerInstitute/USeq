package gata.aligner;

import javax.swing.*;
import java.awt.*;

/**
 * @author nix
 * Frame for holding AlignerInputPanel, the place to fetch user inputs for GATAligner*/
public class AlignerInputForm extends JFrame {
	//fields
	AlignerInputPanel panel;
	
	public AlignerInputForm(GATAligner gata) {
		//modify frame
		setResizable(true);
		setTitle("GATAliner Input Parameters");

		//add panel
		panel = new AlignerInputPanel(this, gata);
		Container contentPane = getContentPane();
		contentPane.add(panel, BorderLayout.CENTER);

		//set size of frame from panel info
		setSize(panel.getFurthestRight(), panel.getFurthestDown() + 25);

	}
	public AlignerInputPanel getPanel() {
		return panel;
	}

}
