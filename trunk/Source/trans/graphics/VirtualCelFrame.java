package trans.graphics;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.io.*;

/**Frame for {@link VirtualCel}.*/
public class VirtualCelFrame extends JFrame{
	//fields
	private VirtualCelPanel panel;
	private JScrollPane scrollPane;
	private Container contentPane;
	private CelMasker celMasker;
	
	//constructors
	public VirtualCelFrame(){
	}
	
	public VirtualCelFrame(float[][] intensities, int maxIntensityValueToPlot, float maskValue){
		makePanel(intensities, maxIntensityValueToPlot, maskValue);
	}
	
	//methods
	public void makePanel(float[][] intensities, int maxIntensityValueToPlot, float maskValue){
		panel = new VirtualCelPanel(intensities, maxIntensityValueToPlot, maskValue, this);
		panel.setBackground(Color.BLACK);
		panel.setPreferredSize(new Dimension(intensities.length, intensities.length));
		setSize(500,500);
		
		scrollPane = new JScrollPane(panel);
		//kill scroll bar movements with arrow keys
		InputMap im = scrollPane.getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
	    im.put(KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, 0), "none");
	    im.put(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, 0), "none");
	    im.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, 0), "none");
	    im.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, 0), "none");
		contentPane = getContentPane();
		contentPane.add(scrollPane);
	}
	
	public void makePNGPanel(float[][] intensities, int maxIntensityValueToPlot, File pngFile){
		panel = new VirtualCelPanel(intensities, maxIntensityValueToPlot, pngFile);
	}

	public JScrollPane getScrollPane() {
		return scrollPane;
	}

	public CelMasker getCelMasker() {
		return celMasker;
	}

	public void setCelMasker(CelMasker celMasker) {
		this.celMasker = celMasker;
	}	
}