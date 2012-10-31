package gata.main;
import javax.swing.*;
import java.io.*;
import javax.swing.event.*;
import java.awt.*;

/**
 * @author Nix
 * Creates a free standing html text window with active links using a url like 
 * file:./docs/gataLinerHelp.html */
public class TextWindow extends JFrame {
	private JEditorPane editorPane;
	
	public TextWindow(String title, int height, int width, String url) {
		setTitle(title);
		setSize(height, width);
		setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
		editorPane = new JEditorPane();
		editorPane.setEditable(false);
		editorPane.addHyperlinkListener(new HyperlinkListener() {
			public void hyperlinkUpdate(HyperlinkEvent event) {
				if (event.getEventType()
					== HyperlinkEvent.EventType.ACTIVATED) {
					try {
						editorPane.setPage(event.getURL());
					} catch (IOException e) {
						editorPane.setText("Problem! Are you off line? " + e);
						e.printStackTrace();
					}
				}
			}
		});
		Container contentPane = getContentPane();
		contentPane.add(new JScrollPane(editorPane), BorderLayout.CENTER);
		loadNewURL(url);
	}
	public void loadNewURL(String url){
		try {
			editorPane.setPage(url);
		} 
		catch (IOException e) {
			editorPane.setText("Problem! Are you off line? " + e);
			e.printStackTrace();
		}
		show();
	}
}