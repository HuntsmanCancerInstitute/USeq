package util.apps;

import java.awt.Color;
import java.awt.Container;

import javax.swing.JFrame;

/**
 * Helper class for {@link ScatterPlot}, the frame.
 *
 */
public class ScatterDrawFrame extends JFrame{
	public ScatterDrawFrame(String[] args){
		setSize(700,700);
		
		ScatterDrawPanel panel = new ScatterDrawPanel(args);
		panel.setBackground(Color.WHITE);
		Container contentPane = getContentPane();
		contentPane.add(panel);
	}

}
