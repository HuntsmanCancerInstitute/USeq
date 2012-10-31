package gata.aligner;

import javax.swing.*;

/**
 * @author nix
 * This is the place to launch the GATAligner, multithreaded, fire at will.*/

public class GATAligner {
	//fields
	AlignerInputForm aif;

	public static void main(String[] args) {
		new GATAligner();
	}
	public GATAligner() {
		//make Aligner input panel to collect user inputs, the frame fires the panel
		aif = new AlignerInputForm(this);
		aif.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		aif.show();
		//make a Console to hold info
		int x = aif.getX();
		int y = aif.getY() + aif.getHeight() + 5;
		int width = aif.getWidth();
	}
	public void lauchAligner() {
		//make progress monitor
		String name = aif.getPanel().getAlignerPrefReference().getBaseName();
		ProgressMonitor pm =
			new ProgressMonitor(
				aif,
				"Running GATAligner for " + name,
				"Launching...",
				0,
				100);
		pm.setMillisToDecideToPopup(0);
		pm.setMillisToPopup(0);
		AlignerThread thread = new AlignerThread(pm, aif);
		thread.start();
	}
}

