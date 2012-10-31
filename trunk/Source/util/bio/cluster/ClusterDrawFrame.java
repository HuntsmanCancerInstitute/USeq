package util.bio.cluster;

import java.awt.*;
import javax.swing.*;

/**
 * Helper class for {@link HierarchicalClustering}, the frame. Creates panel and scroll pane.
 *
 */
public class ClusterDrawFrame extends JFrame{
	
	public ClusterDrawFrame(HierarchicalClustering vhc){
		ClusterDrawPanel panel = new ClusterDrawPanel(vhc);
		int panelWidth = (int)panel.getPanelWidth();
		int panelHeight = (int)panel.getPanelHeight();		
		panel.setPreferredSize(new Dimension(panelWidth, panelHeight));
		JScrollPane scrollPane = new JScrollPane(panel);
		Container contentPane = getContentPane();
		contentPane.add(scrollPane);
		if (panelWidth> 1000) panelWidth = 1000;
		if (panelHeight> 750) panelHeight = 750;
		setSize (panelWidth+30, panelHeight+30);
		this.setTitle(vhc.getTitle());
		setVisible(true);
		setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
	}

}
