package gata.plotter;

import gata.aligner.*;
import gata.main.*;

import java.awt.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;

/**
 * 
 * @author nix
 * Class representing the sliders used to show or hide box-line-box shapes in the DNA alignment
 */
public class Sliders {
	
	//fields
	private double minValue;
	private double maxValue;
	private double range;
	private JSlider slider1;
	private JSlider slider2;
	private AlignParams AP;
	private GATAParams params;
	private AlignPanel aAlignPanel;
	private ToolsFrame toolsFrame;
	//constructor
	
	public Sliders(AlignParams ap, GATAParams params) {
		AP = ap;
		this.params = params;
		toolsFrame = params.getToolsFrameRef();
		minValue = AP.getMIN_SCORE();
		maxValue = AP.getMATCH() * AP.getWIN_SIZE();
		range = maxValue - minValue;
		//Make objects and add them to panel
		slider1 = new JSlider(0, 100, 0);
		slider2 = new JSlider(0, 100, 100);
		makeSlider(slider1, 10, 2);
		makeSlider(slider2, 10, 2);
		SliderHandler sh = new SliderHandler();
		slider1.addChangeListener(sh);
		slider2.addChangeListener(sh);
		aAlignPanel = params.getAlignPanel();
	}
	public JSlider[] getSliders() {
		return new JSlider[] { slider1, slider2 };
	}
	//helper method for sliders
	public void makeSlider(JSlider s, int major, int minor) {
		s.setOrientation(1);
		s.setBackground(Color.WHITE);
		s.setBorder(new EmptyBorder(2, 2, 2, 2));
		s.setPaintTicks(true);
		s.setPaintLabels(true);
		s.setMajorTickSpacing(major);
	}
	private class SliderHandler implements ChangeListener {
		int valSlider1;
		int valSlider2;
		int min = 0;
		int max = 100;
		public void stateChanged(javax.swing.event.ChangeEvent changeEvent) {
			//check to see if slider is still being moved
			JSlider source = (JSlider) changeEvent.getSource();
			//if (source.getValueIsAdjusting() == false){  //reinstate to get fast movement
			//check to see if sliders have been moved if so then print value
			valSlider1 = slider1.getValue();
			valSlider2 = slider2.getValue();
			if (valSlider1 != min) {
				double S = ((((double) valSlider1) / 100) * range) + minValue;
				aAlignPanel.setMinScore(S);
				min = valSlider1;
				toolsFrame.setScores(S, true);
			}
			else if (valSlider2 != max) {
				double S = ((((double) valSlider2) / 100) * range) + minValue;
				aAlignPanel.setMaxScore(S);
				max = valSlider2;
				toolsFrame.setScores(S, false);
			}
		}
	}
}