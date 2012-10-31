package trans.qc;

import java.awt.*;
import javax.swing.*;
import org.jfree.chart.*;


/**Frame for JFree Panels.*/
public class ChartFrame extends JFrame{
	
	public ChartFrame (JFreeChart[] charts, String title, int numRows, int numColumns){
		int num = charts.length;
		setSize(1400,1000);
		setTitle(title);
		
		//set grid layout
		Container contentPane = getContentPane();
		contentPane.setLayout( new GridLayout( numRows, numColumns ) );
		//make and load panels into content pane
		for (int i=0; i< num; i++){
			contentPane.add(new ChartPanel(charts[i]));
		}
		setVisible(true);
		setDefaultCloseOperation( JFrame.EXIT_ON_CLOSE );
		
	}
}