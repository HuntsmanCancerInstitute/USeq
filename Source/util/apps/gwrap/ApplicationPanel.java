package util.apps.gwrap;

import info.clearthought.layout.TableLayout;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

@SuppressWarnings("serial")
public class ApplicationPanel  extends JPanel implements ComponentListener{

	//fields
	public static final Font FONT = new Font("Dialog",Font.PLAIN,12);
	public static final Font FONT_BOLD = new Font("Dialog",Font.BOLD,12);
	public static final Font FONT_ITALICS = new Font("Dialog", Font.ITALIC,12);
	public double panelHeight = 0;
	
	public ApplicationPanel(GWrap_GUI_ClickMe bWrapper) {
		double b = 10;
		double p = TableLayout.PREFERRED;
		//double size[][] = {{5, 250, 24, 5, 20, 5, 435, 20}, {b, p, b}};  // Layout: {{Columns}, {Rows}};
		double size[][] = {{5, TableLayout.FILL,5, 24, 5, 20, 5, 500, 20}, {b, p, b}};  // Layout: {{Columns}, {Rows}};
		TableLayout layout = new TableLayout(size);
		setLayout (layout);
		setBackground(Color.WHITE);
		
		//add results 

		// Add no text instruction
		layout.insertRow (2, TableLayout.PREFERRED);
		layout.insertRow (2, 6);
		JTextArea inst = new JTextArea("For parameters requiring no text, enter a 'Y'.");
		inst.setEditable(false);
		inst.setFont(ApplicationPanel.FONT_ITALICS);
		inst.setBackground(Color.WHITE);
		add(inst, "1, 3, 6, 3");

		// Add app options
		Option[] options = bWrapper.getOptions();      
		for (int i=options.length-1; i>=0; i--){
			// Add a row for spacing and a row for the color boxes
			layout.insertRow (2, TableLayout.PREFERRED);
			layout.insertRow (2, 6);
			add (options[i].getUserInputS(), "1, 3");
			if(options[i].getBrowseBtn() != null) {
				add(options[i].getBrowseBtn(), "3, 3, c, c");
			}
			add (options[i].getFlagS(), "5, 3");
			add (options[i].getDescriptionS(), "7, 3");            
		}
		// Add spacer
		layout.insertRow (2, TableLayout.PREFERRED);
		layout.insertRow (2, 5);
		add(new JTextArea(" "), "1, 3, 7, 3");
		
		// Add description of app
		layout.insertRow (2, TableLayout.PREFERRED);
		layout.insertRow (2, 5);
		add(bWrapper.getAppDescriptionS(), "1, 3, 7, 3");   

		// Add app title
		layout.insertRow (2, TableLayout.PREFERRED);
		layout.insertRow (2, 5);
		add(bWrapper.getAppName(), "1, 3, 7, 3");
		
		
	}

	public void componentHidden(ComponentEvent e) {}
	public void componentMoved(ComponentEvent e) {}
	public void componentShown(ComponentEvent e) {}

	public void componentResized(ComponentEvent e) {
		/*
	        if(frame != null) {
	          Dimension frameSize = frame.getSize();
	          Dimension panelSize = getPreferredSize();
	          panelSize.width = frameSize.width-20;
	          this.setPreferredSize(panelSize);	          
	        }
		 */
	}
}