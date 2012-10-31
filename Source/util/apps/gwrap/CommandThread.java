package util.apps.gwrap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

import javax.swing.JFrame;
import javax.swing.JOptionPane;
import util.gen.Misc;

public class CommandThread extends Thread {

	private String[] command;
	private String[] commandDisplay;
	private String[] saveParams;
	private String fileName;
	private String[] helpCommand;
	private String[] result;
	private GWrap_GUI_ClickMe mainApp;
	private String runStartTime;
	private String commandThreadID;
	private String status;
	private Process p = null;



	public void killProcess() {
		if(p != null) {
			p.destroy();
			p = null;
		}
	}

	public CommandThread(GWrap_GUI_ClickMe aFrame) {
		mainApp = aFrame;
		status = "Running";
	}

	public void run() {   
		try {
			result = executeCommandLineReturnAll(command);
		} catch (Exception e1) {
			result = new String[1];
			result[0] = "Error running command: " + e1.getMessage();
		}
		if(result == null) {
			result = new String[1];
			result[0] = "Error: unable to run command.";   
		}

		File historyDirectory = mainApp.getPrefsDialog().getHistoryDirectory();
		if(!historyDirectory.exists()) {
			boolean fCreated = historyDirectory.mkdir();
			if(!fCreated) {
				showError("Error: Unable to create history directory");
				return;        
			}      
		}

		FileWriter fw = null;
		BufferedWriter writer = null;
		FileReader fr = null;
		BufferedReader reader = null;

		// Save params to history file
		try {
			// Create a temporary file to write the parameters to
			String tempFileName = getCommandThreadID()+".tmp";
			File outFile = new File(historyDirectory, tempFileName);        

			fw = new FileWriter(outFile);
			writer = new BufferedWriter(fw);

			for (int i = 0; i < saveParams.length; i++) {
				writer.write(saveParams[i]);
				writer.newLine();
			}

			writer.write("*********** " + this.getRunStartTime()+" ***********");
			writer.newLine();
			for (int i = 0; i < commandDisplay.length; i++) {
				writer.write(commandDisplay[i] + " ");	
			} 
			writer.newLine();
			writer.newLine();
			for (int i = 0; i < result.length; i++) {
				writer.write(result[i]);
				writer.newLine();
			}
			writer.newLine();
			File inFile = new File(historyDirectory, fileName+".txt");

			if(inFile.exists()) {
				// If file already exists then append its content to the new temp file
				fr = new FileReader(inFile);
				reader = new BufferedReader(fr);
				String thisLine = "";
				while ((thisLine = reader.readLine()) != null) {
					writer.write(thisLine);
					writer.newLine();            
				}
				reader.close();
				reader = null;
				fr.close();
				fr = null;
			}
			writer.flush();
			writer.close();
			writer = null;
			fw.close();
			fw = null;
			inFile.delete();
			outFile.renameTo(inFile);  

		} catch (FileNotFoundException e) {
			showError("Error: FileNotFoundException");       
		} catch (IOException e) {
			showError("Error: IOException");       
		} finally {
			try {
				if (writer != null) {
					writer.close();
				}
				if (fw != null) {
					fw.close();
				}
				if (reader != null) {
					reader.close();
				}
				if (fr != null) {
					fr.close();
				}                  
			} catch (IOException e) {
				showError("Error: IOException");       
			}
		}
		status = "Completed";
		JOptionPane.showMessageDialog(mainApp, fileName + " launched "+runStartTime+" is complete.");
		mainApp.getJobsDialog().refreshJobsList();
		mainApp.getResultsDialog().refreshResultsWindow(fileName);
		
	}

	/**Executes tokenized params on command line, use full paths.
	 * Put each param in its own String.  
	 * Returns null if a problem is encountered.
	 * Returns both error and data from execution.
	 */
	public String[] executeCommandLineReturnAll(String[] command){
		ArrayList<String> al = new ArrayList<String>();
		try {
			Runtime rt = Runtime.getRuntime();
			p = rt.exec(command);
			BufferedReader data = new BufferedReader(new InputStreamReader(p.getInputStream()));
			BufferedReader error = new BufferedReader(new InputStreamReader(p.getErrorStream())); //for debugging
			String line;
			while ((line = data.readLine()) != null){
				al.add(line);
			}
			while ((line = error.readLine()) != null){
				al.add(line);
			}
			data.close();
			error.close();

		} catch (Exception e) {
			System.out.println("Problem executing -> "+Misc.stringArrayToString(command," ")+" "+e.getLocalizedMessage());
			e.printStackTrace();
			return null;
		}
		String[] res = new String[al.size()];
		al.toArray(res);
		return res;
	}   

	private void showError(String errorMsg) {
		JOptionPane.showMessageDialog(new JFrame(), errorMsg, "Dialog",
				JOptionPane.ERROR_MESSAGE); 

	}
	
	public String[] getCommandDisplay() {
		return commandDisplay;
	}

	public void setCommandDisplay(String[] commandDisplay) {
		this.commandDisplay = commandDisplay;
	}
	
	public String getStatus() {
		return status;
	}

	public void setStatus(String status) {
		this.status = status;
	}

	public String getCommandThreadID() {
		return commandThreadID;
	}

	public void setCommandThreadID(String commandThreadID) {
		this.commandThreadID = commandThreadID;
	}

	public String getRunStartTime() {
		return runStartTime;
	}

	public void setRunStartTime(String runStartTime) {
		this.runStartTime = runStartTime;
	}
	public String getFileName() {
		return fileName;
	}

	public void setFileName(String fileName) {
		this.fileName = fileName;
	}

	public String[] getSaveParams() {
		return saveParams;
	}

	public void setSaveParams(String[] saveParams) {
		this.saveParams = saveParams;
	}

	public String[] getCommand() {
		return command;
	}

	public void setCommand(String[] command) {
		this.command = command;
	}

	public String[] getHelpCommand() {
		return helpCommand;
	}

	public void setHelpCommand(String[] helpCommand) {
		this.helpCommand = helpCommand;
	}

}
