package trans.graphics;
import java.awt.*;
import org.jfree.chart.*;

import javax.swing.*;
import util.gen.*;
import java.io.*;



/**Frame for {@link RankedSetAnalysis}.*/
public class RankedSetDrawFrame extends JFrame{
	private RankedSetDrawPanel panel;
	private JScrollPane scrollPane;
	private Container contentPane;
	private int heightAdjuster = 50;
	
	public RankedSetDrawFrame(RankedSetAnalysis main){
		panel = new RankedSetDrawPanel(main);
		panel.setBackground(Color.WHITE);
		panel.setPreferredSize(new Dimension((int)panel.getPanelWidth(), (int)panel.getPanelHeight()));
		int height = (int)panel.getPanelHeight();
		if (height > 1000) height = 1000;
		setSize((int)panel.getPanelWidth()*2,height + heightAdjuster);
		setTitle(main.getSetOneFile().getName()+" vs "+main.getSetTwoFile().getName());
		scrollPane = new JScrollPane(panel);
		ChartPanel cp = main.getChartPanel();
		JSplitPane sp = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,scrollPane, cp);
		contentPane = getContentPane();
		contentPane.add(sp);
		setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
		setVisible(true);
		sp.setDividerLocation(0.5);
		//save pngs?
		if (main.isSavePngs()) {
			String name = IO.getFullPathName(main.getSetOneFile())+"_vs_"+main.getSetTwoFile().getName();
			panel.saveBufferedImage(1,new File ( name+"_BoxLine.png"), false);
			try {
				ChartUtilities.saveChartAsPNG(new File(name+ "_IntGraph.png"), main.getChart(), 1000, 700);
			} catch (IOException e) {
				System.out.println("Problem occurred creating chart.");
			}
		}
	}
}