package util.gen;

import java.awt.*;

import javax.swing.*;

import java.awt.event.*;
import java.awt.font.FontRenderContext;
import java.awt.font.TextLayout;
import java.io.*;


/**
 * Methods related to Swing apps.
 */
public class Swing {

	/** Map of hints for high quality/ slow rendering*/
	public static RenderingHints fetchRenderingHints () {
		RenderingHints hints = new RenderingHints(null);
		hints.put(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		hints.put(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		hints.put(RenderingHints.KEY_DITHERING, RenderingHints.VALUE_DITHER_DISABLE);
		hints.put(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
		hints.put(RenderingHints.KEY_FRACTIONALMETRICS, RenderingHints.VALUE_FRACTIONALMETRICS_ON);
		hints.put(RenderingHints.KEY_ALPHA_INTERPOLATION, RenderingHints.VALUE_ALPHA_INTERPOLATION_QUALITY);
		hints.put(RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_QUALITY);
		return hints;
	}
	
	/**Draws a JLabel on a Graphics2D. Set the tuningFraction to 1.05 for standard drawing or adjust to modify the y displacement.
	 * For some blasted reason, JLabels don't draw for buffered images when the layout manager is null. Thus this work around.*/
	public static void drawJLabel(Graphics2D g2, JLabel label, float tuningFraction){
		//fields
		String text = label.getText();
		Rectangle bounds = label.getBounds();
		Font font = label.getFont();
		float centerY = new Double (bounds.getCenterY()).floatValue();
		g2.setFont(font);
		g2.setColor(label.getForeground());
		FontRenderContext context = g2.getFontRenderContext();
		TextLayout layout = new TextLayout(text, font, context);
		float offset = layout.getAscent()/2.0f - layout.getDescent();
		float adjustedCenterY = centerY + offset*tuningFraction;
		g2.drawString(text, new Double(bounds.getX()).floatValue(), adjustedCenterY);
	}

	public static double makeDoublePercentInputDialog(
		String title,
		String comment,
		double suggestion) {
		Object ob =
			JOptionPane.showInputDialog(null,comment,title,JOptionPane.PLAIN_MESSAGE,null,null,Double.toString(suggestion));
		if (ob == null) { // if user hits cancel or closes the window
			return suggestion;
		}
		double num;
		try{
			num = Double.parseDouble((String) ob);
		}
		catch (NumberFormatException e){
			num = makeDoublePercentInputDialog(title,"Decimals Only! "+comment,suggestion);
		}
		if (num >1 || num<0) return makeDoublePercentInputDialog(title,"Fractional Decimal Only (0-1)! "+comment,suggestion);
		else  return num;
	}

	public static double makeDoubleInputDialog(
		String title,
		String comment,
		double suggestion) {
		Object ob =
			JOptionPane.showInputDialog(null,comment,title,JOptionPane.PLAIN_MESSAGE,null,null,Double.toString(suggestion));
		if (ob == null) { // if user hits cancel or closes the window
			return suggestion;
		}
		double num;
		try{
			num = Double.parseDouble((String) ob);
		}
		catch (NumberFormatException e){
			num = makeDoubleInputDialog(title,"Problem parsing number! "+ob+" "+comment,suggestion);
		}
		return num;
	}
		
	public static int makeIntInputDialog(String title,String comment,int suggestion) {
		//(Component parentComponent, Object message, String title, int messageType,          Icon icon, Object[] selectionValues, Object initialSelectionValue)
		Object ob =JOptionPane.showInputDialog(null,comment,title,  JOptionPane.PLAIN_MESSAGE,null,      null,                     Integer.toString(suggestion));
		if (ob == null) { // if user hits cancel or closes the window
			return suggestion;
		}
		try{
			return Integer.parseInt((String) ob);
		}
		catch (NumberFormatException e){
			return makeIntInputDialog(title,"Integers only! "+comment,suggestion);
		}
		
	}


	public static void throwWarning(
		JFrame frame,
		JTextField field,
		String text) {
		JOptionPane.showMessageDialog(
			frame,
			text,
			null,
			JOptionPane.WARNING_MESSAGE);
		field.requestFocusInWindow();
	}
	
	public static void throwError(JFrame frame, String text) {
			JOptionPane.showMessageDialog(
				frame,
				text,
				null,
				JOptionPane.ERROR_MESSAGE);
		}
	
	public static int[] setButton(JButton but, int x, int y) {
		Dimension dim = but.getPreferredSize();
		int[] wh = {(int) dim.getWidth(), (int) dim.getHeight()};
		but.setBounds(x, y, wh[0], wh[1]);
		return wh;
	}

	public static JButton makeButton(
		String text,
		int size,
		ActionListener action,
		JPanel panel) {

		JButton but = new JButton(text);
		but.setFont(new Font("Dialog", 0, size));
		but.setBackground(Color.WHITE);
		panel.add(but);

		but.addActionListener(action);
		return but;
	}

	public static void showFileChooserDialogBox(
		JTextField fieldRef,
		JFrame frame,
		JFileChooser chooser) {
		//int result = chooser.showDialog(frame, "Choose"); //problems with selecting folders!
		int result = chooser.showOpenDialog(frame);
		if (result == JFileChooser.APPROVE_OPTION) {
			File file = chooser.getSelectedFile();
			String text = "Problem fetching this file/directory!";
			try {
				text = file.getCanonicalPath();
			} catch (IOException IOe) {
				System.out.println("Problem with getting path/file info!");
				IOe.printStackTrace();
			}
			fieldRef.setText(text);
			chooser.setCurrentDirectory(file);
		}
		fieldRef.requestFocusInWindow();
	}

	public static JTextField makeFieldSimple(
		JPanel panel,String text) {
		JTextField tex = new JTextField(text);
		panel.add(tex);
		return tex;
	}
	public static JLabel makeLabel(
		JPanel panel,
		String text,
		int size,
		int thickness) {
		JLabel lab = new JLabel(text);
		lab.setFont(new Font("Dialog", thickness, size));
		panel.add(lab);
		return lab;
	}

	public static int[] setField(JTextField field, int x, int y, int width) {
		Dimension dim = field.getPreferredSize();
		int[] wh = { width, (int) dim.getHeight()};
		field.setBounds(x, y, wh[0], wh[1]);
		return wh;
	}
	public static int[] setLabel(JLabel lab, int x, int y) {
		Dimension dim = lab.getPreferredSize();
		int[] wh = {(int) dim.getWidth(), (int) dim.getHeight()};
		lab.setBounds(x, y, wh[0], wh[1]);
		return wh;
	}




}
