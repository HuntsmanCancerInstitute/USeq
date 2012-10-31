package edu.utah.bass;

import java.awt.*;

import javax.swing.*;


public class SequenceGraphFrame extends JFrame{

	
	//fields
	private SequenceGraph sequenceGraph;
	private static final long serialVersionUID = 1L;
	private int heightBuffer = 22;
	
	//constructors
	public SequenceGraphFrame(SequenceGraph sequenceGraph){
		this.sequenceGraph = sequenceGraph;
		makePanel();
	}
	
	//methods
	public void makePanel(){
		SequenceGraphPanel panel = new SequenceGraphPanel(sequenceGraph);
		getContentPane().add(panel);
		setSize((int)panel.getPanelWidth(),(int)panel.getPanelHeight() + heightBuffer);
	}
	
	
	
	
	
	
	
}
