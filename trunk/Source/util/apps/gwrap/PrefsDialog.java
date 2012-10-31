/*
 * Copyright (c) 1995 - 2008 Sun Microsystems, Inc.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   - Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *   - Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *
 *   - Neither the text of Sun Microsystems nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */ 

package util.apps.gwrap;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JTextField;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.Properties;
import java.util.regex.Pattern;

import util.gen.FileDrop;
import util.gen.IO;

class PrefsDialog extends JDialog {

	private static final long serialVersionUID = 1L;
	private JTextField javaPathField;
	private JTextField memoryParamField;
	private Container pane;
	private String jpActionCommand = "jpBrowse";
	private String saveActionCommand = "save";
	private String cancelActionCommand = "cancel";
	private GWrap_GUI_ClickMe mainApp;
	private static final Pattern memoryPattern = Pattern.compile("-Xmx\\d+[MGmg]");




	/** Creates the reusable dialog. */
	public PrefsDialog(GWrap_GUI_ClickMe aFrame) {
		super(aFrame, true);
		mainApp = aFrame;
System.out.println("Making prefs dialog");
		setTitle("Preferences");

		pane = getContentPane();
		pane.setLayout(null);
		pane.setPreferredSize(new Dimension(200,200));

		PrefsButtonListener prefsButtonListener = new PrefsButtonListener();

		// Java Path Inro
		JLabel jpLabel = new JLabel("Full Path to Java (e.g. /usr/bin/java)");        
		jpLabel.setBounds(10, 10, 500, 20);
		add(jpLabel);
		javaPathField = new JTextField(60);
		javaPathField.setBounds(10, 30, 490, 20);
		add(javaPathField);

		new FileDrop( javaPathField, new FileDrop.Listener() {   
			public void filesDropped( java.io.File[] files ) {  
				//add first file
				if (files.length !=0) javaPathField.setText(files[0].getAbsolutePath());
			}  
		});         

		JButton jpBrowseBtn = new JButton("Browse");
		jpBrowseBtn.setActionCommand(jpActionCommand);
		jpBrowseBtn.setBounds(500, 30, 80, 18);
		jpBrowseBtn.addActionListener(prefsButtonListener);
		add(jpBrowseBtn);

		// Java Memory Param Info
		JLabel mpLabel = new JLabel("Java Memory Parameter (e.g. -Xmx2G)");        
		mpLabel.setBounds(10, 60, 500, 20);
		add(mpLabel);
		memoryParamField = new JTextField(60);
		memoryParamField.setBounds(10, 80, 490, 20);
		add(memoryParamField);        

		JButton saveBtn = new JButton("Save");
		saveBtn.setActionCommand(saveActionCommand);
		saveBtn.setBounds(195, 120, 80, 20);
		saveBtn.addActionListener(prefsButtonListener);
		add(saveBtn);

		JButton cancelBtn = new JButton("Cancel");
		cancelBtn.setActionCommand(saveActionCommand);
		cancelBtn.setBounds(310, 120, 80, 20);
		cancelBtn.addActionListener(prefsButtonListener);
		add(cancelBtn);
	}

	public void checkPreferences() {	
System.out.println("CheckPreferences() called");		
		//java
		String javaPath = javaPathField.getText().trim();
		if (checkJava(javaPath) == false) return;
		//memory
		String memoryParam = memoryParamField.getText().trim();
		if (checkMemoryParameter(memoryParam) == false) return;
	}


	/**Checks to see java can be executed and if the current java is 1.6, 1.7, or 1.8 by executing and parsing 'java -version'*/
	public boolean checkJava(String java){
			String[] cmd = {java,"-version"};
			String[] results = IO.executeCommandLineNoError(cmd);
			if (results == null || results.length ==0) {
				//cannot find or execute java
				mainApp.throwError("Cannot find or execute java? -> "+java);
				return false;
			}
			else {
				for (int i=0; i< results.length; i++){
					String line = results[i].toLowerCase().trim();
					if (line.startsWith("java version")){
						if (line.contains("1.6.") || line.contains("1.7.")|| line.contains("1.8.")) return true;
						mainApp.throwError("Your java application is not >= 1.6 (type 'java -version' on the cmd line). Install the most recent java from http://www.java.com/en/download/ .\n");
						return false;
					}
				}
			}
			return false;
	}

	public boolean checkMemoryParameter(String p){
		boolean good =  memoryPattern.matcher(p).matches();
		if (good == false) {
			mainApp.throwError("Your Java Memory Parameter is malformed -> "+p);
		}
		return good;
	}

	private class PrefsButtonListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			//save
			if(e.getActionCommand().compareTo(saveActionCommand)==0) {
				//java
				String javaPath = javaPathField.getText().trim();
				if (checkJava(javaPath) == false) return;
				//memory
				String memoryParam = memoryParamField.getText().trim();
				if (checkMemoryParameter(memoryParam) == false) return;

				//save props
				Properties props = new Properties();
				try{
					props.setProperty("javaPath", javaPath);
					props.setProperty("memoryParams", memoryParam);
					FileOutputStream fos = new FileOutputStream(new File(mainApp.getSaveDirectory(),"config.properties"));
					props.store(fos, "GWrap Properties File.");
					fos.close();
				} catch(Exception e1) {
					e1.printStackTrace();           
				}         
				clearAndHide(); 
				return;
			}

			//cancel
			else if(e.getActionCommand().compareTo(cancelActionCommand)==0) {
				clearAndHide(); 
				return;
			}

			//file/ directory chooser
			else if(e.getActionCommand().compareTo(jpActionCommand)==0) {
				// Show open dialog; this method does not return until the dialog is closed
				//ThreadedFileChooser ft = new ThreadedFileChooser(this, (e.getActionCommand().compareTo(jpActionCommand)==0));
				javax.swing.JFileChooser fc = new javax.swing.JFileChooser();
				if(e.getActionCommand().compareTo(jpActionCommand)==0) {
					fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
				} else {
					fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				}
				if (fc.showOpenDialog(pane) == JFileChooser.APPROVE_OPTION) {
					String path = fc.getSelectedFile().getPath();
					if(e.getActionCommand().compareTo(jpActionCommand)==0) {
						javaPathField.setText(path);               
					}      
				}
				return;
			}
		}
	}    

	public void initTextFields() {
System.out.println("initTextFields() called");		
		String javaPathText = null;
		String memoryParamText = null;
		javaPathText = getProperty("javaPath");
		javaPathField.setText(javaPathText);
		memoryParamText = getProperty("memoryParams");
		memoryParamField.setText(memoryParamText);
	}

	/** This method clears the dialog and hides it. */
	public void clearAndHide() {
		setVisible(false);
	}

	/**
	 * Returns a property from the "config.properties" file.
	 * If the property is not found returns an empty text.  No Exceptions will be thrown.
	 * @param propertyName
	 * @return
	 */
	public String getProperty(String propertyName) {
		//set default
		String defaultValue = "";
		if(propertyName.compareTo("javaPath")==0) {
			defaultValue = "java";
		}
		else if(propertyName.compareTo("memoryParams")==0) {
			defaultValue = "-Xmx2G";
		}      
		//attempt to load from file
		Properties props = new Properties();
		try{
			props.load(new FileInputStream(new File(mainApp.getSaveDirectory(),"config.properties")));
		} catch(Exception e) {
			return defaultValue;
		}
		String toReturn = props.getProperty(propertyName);
		//if not found then return default
		if(isNullOrBlank(toReturn)){
			return defaultValue;
		}
		return toReturn;
	}


	public File getHistoryDirectory() {      
		File h = new File(mainApp.getSaveDirectory(),"History");
		if (h.exists() == false) h.mkdir();
		return h;
	}

	public static boolean isNullOrBlank(String str){
		if(str == null || str.trim().equals("")) return true;
		return false;
	}



}

