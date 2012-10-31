package gata.main;
import javax.swing.*;
import java.io.*;
import javax.swing.event.*;
import java.awt.*;

/**
 * @author Nix
 * Creates a free standing html text window with active links using a text containing html*/
public class SimpleTextWindow extends JFrame {
	public SimpleTextWindow(String title, int height, int width, String text) {
		setTitle(title);
		setBounds(100,100,height, width);
		//setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		final JEditorPane editorPane = new JEditorPane("text/html",text);
		editorPane.setEditable(false);
		editorPane.addHyperlinkListener(new HyperlinkListener() {
			public void hyperlinkUpdate(HyperlinkEvent event) {
				if (event.getEventType()
					== HyperlinkEvent.EventType.ACTIVATED) {
					try {
						editorPane.setPage(event.getURL());
					} catch (IOException e) {
						editorPane.setText("Problem! " + e);
						e.printStackTrace();
					}
				}
			}
		});

		Container contentPane = getContentPane();
		contentPane.add(new JScrollPane(editorPane), BorderLayout.CENTER);
		show();
	}
}