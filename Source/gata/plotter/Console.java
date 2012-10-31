package gata.plotter;

import java.awt.*;
import javax.swing.*;

/**
 * @author Nix
 	 Console creates a free standing window to which you can write strings
	 The window can be resized.  The text scrolls off the top and is unlimited
	 in length.
	 initialize with #rows, #columns, x coordinate, y coordinate
	 */
public class Console extends JFrame{
    private JTextArea textArea;
    public Console (int locX, int locY, int width, int height, String title){
		textArea = new JTextArea();
    	
        //for text area
        textArea.setEditable(false);
        getContentPane().add(new JScrollPane(textArea), BorderLayout.CENTER);
        textArea.setFont(new Font("Monospaced", Font.BOLD, 12));
        //for JFrame
		setBounds(locX, locY, width, height);
        setTitle(title);
        setResizable(true);
		setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
        
    }
    
    public void printToTextArea(String x){
        textArea.append(x);
        textArea.setCaretPosition(textArea.getDocument().getLength());
    }
    public void setTextArea(String text){
    	textArea.setText("");
    	textArea.setText(text);
    }
    public void updateConsole(){
    	textArea.repaint();
    	textArea.revalidate();
    	repaint();
    }
	public JTextArea getTextArea() {
		return textArea;
	}

}
