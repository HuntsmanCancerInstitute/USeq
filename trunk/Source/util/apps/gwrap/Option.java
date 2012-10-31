package util.apps.gwrap;
import java.awt.*;
import java.awt.dnd.DropTargetDragEvent;
import java.awt.dnd.DropTargetDropEvent;
import java.awt.dnd.DropTargetEvent;
import java.awt.dnd.DropTargetListener;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.regex.*;
import javax.swing.*;
import javax.swing.border.*;

import util.gen.FileDrop;
import util.gen.Swing;

public class Option implements DropTargetListener {
	private String flag;
	private String description;
	private String userInput;
	private JTextArea userInputS;
	private JTextArea flagS;
	private JTextArea descriptionS;
	private JButton browseBtn = null;
	private GWrap_GUI_ClickMe mainApp;
	private boolean hasFile = false;
	private boolean hasDirectory = false;
	public static final Pattern space = Pattern.compile("\\s+");
	public static final Border border = BorderFactory.createLineBorder(Color.black);

	public Option (String line, GWrap_GUI_ClickMe aFrame){
		mainApp = aFrame;
		Matcher mat = space.matcher(line);
		if (mat.find()){
			flag = line.substring(0,2);
			description = line.substring(mat.end());
			String descLower = description.toLowerCase();
			if(descLower.indexOf("file") >= 0) {
				hasFile = true;
			}
			if((descLower.indexOf("directory") >= 0) || (descLower.indexOf("directories") >= 0)) {
				hasDirectory = true;
			}
		}
	}

	public void makeSwingComponents(int optIndex)
	{
		userInputS = new JTextArea();
		userInputS.setFont(ApplicationPanel.FONT);
		userInputS.setColumns(10);
		userInputS.setBackground(Color.white);
		userInputS.setWrapStyleWord(true);
		userInputS.setLineWrap(true);
		userInputS.setBorder(border);

		BrowseButtonListener browseButtonListener = new BrowseButtonListener();

		if(hasFile || hasDirectory) {
			
			browseBtn = new JButton("...");
			browseBtn.setPreferredSize(new Dimension(24, 18));
			browseBtn.setActionCommand("" + optIndex);
			browseBtn.addActionListener(browseButtonListener);   

			new FileDrop(userInputS, BorderFactory.createLineBorder(Color.red), false, new FileDrop.Listener() {   
				public void filesDropped( java.io.File[] files ) {  
					String pathStr = "";
					// handle file drop
					for (int i = 0; i < files.length; i++) {
						if (i > 0) pathStr = pathStr + ",";                
						pathStr = pathStr + files[i].getAbsolutePath();
					}
					//look for spaces in path
					if (pathStr.contains(" ")) Swing.throwError(null, "No spaces in file paths!");
					else userInputS.setText(pathStr);
				}   
			}); 
		}

		flagS = new JTextArea(flag);
		flagS.setFont(ApplicationPanel.FONT);
		flagS.setEditable(false);
		flagS.setBackground(Color.WHITE);

		// Remove carriage returns
		description = description.replace("\n", "");

		String[] vals = description.split(" ");

		StringBuilder sb = new StringBuilder();
		for(String s : vals) { 
			if(s.trim().length() > 0){
				sb.append(s).append(" ");
			}
		}
		description = sb.toString().trim();
		descriptionS = new JTextArea(description);
		descriptionS.setFont(ApplicationPanel.FONT);
		descriptionS.setEditable(false);
		descriptionS.setBackground(Color.WHITE);
		descriptionS.setWrapStyleWord(true);
		descriptionS.setLineWrap(true);

	}

	public JButton getBrowseBtn() {
		return browseBtn;
	}

	public void setBrowseBtn(JButton browseBtn) {
		this.browseBtn = browseBtn;
	}

	private class BrowseButtonListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			//System.out.println(e.getActionCommand());
			javax.swing.JFileChooser fc = new javax.swing.JFileChooser();
			File file = mainApp.getLastFile();
			if(file != null) {
				fc.setCurrentDirectory(file);
			}
			if(hasDirectory) {
				fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);  
			}
			if(hasFile) {
				fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);  
			}

			if (fc.showOpenDialog(mainApp) == JFileChooser.APPROVE_OPTION) {
				file = fc.getSelectedFile();
				mainApp.setLastFile(file);
				String path = file.getPath();
				userInput = path;
				userInputS.setText(path);    
			}
		}
	} 	

	public String toString(){
		if (userInput == null) return flag+" "+description;
		else return flag+" "+description+"\n\t"+userInput;
	}

	public String getFlag() {
		return flag;
	}

	public void setFlag(String flag) {
		this.flag = flag;
	}

	public String getNotes() {
		return description;
	}

	public void setNotes(String notes) {
		this.description = notes;
	}

	public String getUserInput() {
		return userInput;
	}

	public void setUserInput(String userInput) {
		this.userInput = userInput;
	}

	public static Pattern getSpace() {
		return space;
	}

	public String getDescription() {
		return description;
	}

	public void setDescription(String description) {
		this.description = description;
	}

	public JTextArea getUserInputS() {
		return userInputS;
	}

	public void setUserInputS(JTextArea userInputS) {
		this.userInputS = userInputS;
	}

	public JTextArea getFlagS() {
		return flagS;
	}

	public void setFlagS(JTextArea flagS) {
		this.flagS = flagS;
	}

	public JTextArea getDescriptionS() {
		return descriptionS;
	}

	public void setDescriptionS(JTextArea descriptionS) {
		this.descriptionS = descriptionS;
	}

	public void dragEnter(DropTargetDragEvent arg0) {
		// TODO Auto-generated method stub

	}

	public void dragExit(DropTargetEvent arg0) {
		// TODO Auto-generated method stub

	}

	public void dragOver(DropTargetDragEvent arg0) {
		// TODO Auto-generated method stub

	}

	public void drop(DropTargetDropEvent arg0) {
		// TODO Auto-generated method stub

	}

	public void dropActionChanged(DropTargetDragEvent arg0) {
		// TODO Auto-generated method stub

	}

}
