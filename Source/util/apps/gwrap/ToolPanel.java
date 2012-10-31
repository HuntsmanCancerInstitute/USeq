package util.apps.gwrap;

import java.awt.Color;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JPanel;

public class ToolPanel  extends JPanel implements ActionListener{

	private static final long serialVersionUID = 1L;
	//fields
	public static final Font SMALLFONT = new Font("Dialog",Font.PLAIN,10);
	private GWrap_GUI_ClickMe mainFrame;
	private JButton startButton;
	private JButton resetButton;

	public ToolPanel() {

	}

	public ToolPanel(GWrap_GUI_ClickMe aFrame) {
		mainFrame = aFrame;
		setBackground(Color.WHITE);
		startButton = makeButton("Start Job");
		resetButton = makeButton("Reset Form");
	}

	public void enableActionButtons() {
		startButton.setEnabled(true);
		resetButton.setEnabled(true);
	}

	private JButton makeButton(String label){
		JButton button = new JButton(label);
		button.addActionListener(this);
		button.setFont(SMALLFONT);
		button.setEnabled(false);
		add(button);
		return button;
	}

	public void actionPerformed(ActionEvent e) {
		String command = e.getActionCommand();
		if ("Start Job".equals(command)) {
			mainFrame.runCommand();
		} 
		else if ("Reset Form".equals(command)) {
			mainFrame.clearParams();
		}
		else {
			System.out.println("Kill Job");
		}
	}

}
