package util.gen;
import javax.swing.*;

import java.awt.*;

/**
 * Creates a free standing plain editable text window. Don't forget to call show()! */
public class TextFrame extends JFrame {
	//fields
	private JTextArea textArea;
	
	public TextFrame (int locX, int locY, int width, int height, String title){
		textArea = new JTextArea();
		getContentPane().add(new JScrollPane(textArea), BorderLayout.CENTER);
		//for JFrame
		setBounds(locX, locY, width, height);
		setTitle(title);
		setResizable(true);
		setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);	
	}
	
	public void append(String x){
		textArea.append(x);
		textArea.setCaretPosition(textArea.getDocument().getLength());
	}
	public JTextArea getTextArea() {
		return textArea;
	}
}