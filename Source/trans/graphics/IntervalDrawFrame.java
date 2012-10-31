package trans.graphics;

import java.awt.*;
import javax.swing.*;
import trans.main.Interval;
import util.gen.TextFrame;


/**Frame for {@link IntervalPlotter}.*/
public class IntervalDrawFrame extends JFrame{
	private IntervalDrawPanel panel;
	private JScrollPane scrollPane;
	private Container contentPane;
	
	public IntervalDrawFrame(){
	}
	
	public IntervalDrawFrame(Interval interval, int intervalNumber, TextFrame textFrame, boolean showPSPM){
		makePanel(interval, intervalNumber, textFrame, showPSPM);
	}
	public void makePanel(Interval interval, int intervalNumber, TextFrame textFrame, boolean showPSPM){
		panel = new IntervalDrawPanel(interval, intervalNumber, textFrame, showPSPM);
		panel.setBackground(Color.BLACK);
		panel.setPreferredSize(new Dimension((int)panel.getRunningX()+10, (int)panel.getRunningY()+10));
		setSize((int)panel.getRunningX()+50,(int)panel.getRunningY()+60);
		setTitle(panel.getTitle());
		scrollPane = new JScrollPane(panel);
		contentPane = getContentPane();
		contentPane.add(scrollPane);
		
	}
	
	public IntervalDrawPanel getPanel() {
		return panel;
	}
}