package util.apps.gwrap;
import java.awt.Container;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JScrollPane;

class AboutWindow extends JFrame {

	private static final long serialVersionUID = 1L;

	private Container pane;

	private String closeActionCommand = "cancel";


	/** Creates the reusable dialog. */
	public AboutWindow() {

		setTitle("About GWrap GUI Command Line Wrapper");

		pane = getContentPane();
		pane.setLayout(null);
		pane.setBounds(0, 0, 800, 600);

		InstuctionsPanel appPanel = new InstuctionsPanel();

		JScrollPane sp = new JScrollPane(appPanel);
		sp.setBounds(0, 0, 792, 530);

		add(sp);

		HelpButtonListener prefsButtonListener = new HelpButtonListener();



		JButton closeBtn = new JButton("Close");
		closeBtn.setActionCommand(closeActionCommand);
		closeBtn.setBounds(360, 540, 80, 20);
		closeBtn.addActionListener(prefsButtonListener);
		add(closeBtn);
	}


	private class HelpButtonListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			if(e.getActionCommand().compareTo(closeActionCommand)==0) {
				clearAndHide(); 
				return;
			}
		}
	}    


	/** This method clears the dialog and hides it. */
	public void clearAndHide() {
		//textField.setText(null);
		setVisible(false);
	}

}

