package util.apps.gwrap;



import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.net.*;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.security.CodeSource;
import java.util.*;
import java.util.regex.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import javax.swing.*;
import javax.swing.event.*;

import edu.utah.seq.useq.SliceInfo;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Swing;
import info.clearthought.layout.*;

public class GWrap_GUI_ClickMe extends JFrame{

	private GWrap_GUI_ClickMe bw;
	private File jarFile;
	private JTextArea appName;
	private PrefsDialog prefsDialog;
	private JobsDialog jobsDialog;
	private ResultsDialog resultsDialog;
	private String[] menu;
	private JMenu appsMenuAL;
	private JMenu appsMenuMZ;
	private JMenu prefsMenu;
	private JMenu resultsMenu;
	private JMenu jobsMenu;
	private JMenu helpMenu;

	private String appDescription;
	private JTextArea appDescriptionS;
	private Container pane;
	private Option[] options;
	private String java = "java";
	private File lastFile = null;
	public static final Pattern startSpace = Pattern.compile("^\\s+");
	public static final Pattern flag = Pattern.compile("^-[a-z]");
	private boolean debug = false; 
	private ToolPanel toolPanel;
	protected String currentCommand = "";
	protected TreeMap<String, CommandThread> jobs;
	protected HashMap<String, String> helpMenuFiles;
	private static final long serialVersionUID = 1L;
	private File saveDirectory;
	private File appsDirectory;
	private File docsDirectory;

	public static void main (String args[]){
		new GWrap_GUI_ClickMe(args);
	}

	public GWrap_GUI_ClickMe (String args[]){
		//find/make the save directory for preferences and history
		saveDirectory = new File(new File (System.getProperty("user.home")),".GWrap");
		saveDirectory.mkdir();
		
		//locate Apps and Documentation directories
		findDirectories();

		bw = this;
		jobs = new TreeMap<String, CommandThread>(new JobsComparator());
		helpMenuFiles = new HashMap<String, String>();
		setTitle("GWrap Command Line GUI Wrapper");
		pane = getContentPane();

		double b = 1;
		double f = TableLayout.FILL;
		double p = TableLayout.PREFERRED;
		double size[][] = {{f}, {p, b, p, b, p, b}};
		TableLayout layout = new TableLayout(size);
		pane.setLayout (layout);

		//create prefsDialog and load preferences
		prefsDialog = new PrefsDialog(this);

		//add menus
		addMenu(pane);

		//add buttons
		toolPanel = new ToolPanel(this);
		toolPanel.setBackground(new Color(238, 238, 238));
		pane.add (toolPanel, "0, 4");      

		//throw up splash screen with instructions
		showAppInstructions();

		//kill app upon closing
		allowClosing();
		setVisible(true);

		//create other pages
		jobsDialog = new JobsDialog(this);
		resultsDialog = new ResultsDialog(this);     
	}
	
	public void findDirectories(){
		System.out.println("Here findDirectories()");
		//look in current working directory (will be found here if using standard SourceForge release)
		File dir = new File(System.getProperty("user.dir"));
System.out.println("Looking for dirs in "+dir);		
		if (dir != null){
			docsDirectory = new File (dir, "Documentation");
			appsDirectory = new File (dir, "Apps");
			if (checkDocAppDirs()) {
System.out.println("Docs and Apps exist, assigning pointers.");						
				return;
			}
			//attempt to get it from the jar file (e.g. webstart release)
			else {
System.out.println("Attempting to get them from webstart jar file.");				
				copyDocsApps();
				if (checkDocAppDirs() == false){
	System.out.println("appsDir "+appsDirectory);
	System.out.println("docsDir "+docsDirectory);
					throwError("Couldn't find your Apps or Documents directory? Did you move the GWrap_GUI_ClickMe.jar file out of the USeq_xxx dir?  Contact the USeq admin for help.");
				}
			}
		}
		
	}
	
	/**Looks for good Documentation and Apps directories.
	 * @return true if good, false if bad.*/
	public boolean checkDocAppDirs(){
		if (appsDirectory == null || docsDirectory == null) return false;
		File[] f = IO.extractOnlyFiles(appsDirectory);
		if (f == null || f.length == 0) return false;
		f = IO.extractOnlyFiles(docsDirectory);
		if (f == null || f.length == 0) return false;
		return true;
	}
	
	/**Copies the files from the jar file to the .GWrap dir in their home dir.  Needed for webstart since nested jar files not allowed.*/
	public void copyDocsApps(){
System.out.println("Attempting to pull files using copyDocsApps()");
		
		try {
			CodeSource src = GWrap_GUI_ClickMe.class.getProtectionDomain().getCodeSource();
			if( src != null ) {
				//make directories to hold files
				docsDirectory = new File (saveDirectory, "Documentation");
				docsDirectory.mkdir();
				appsDirectory = new File (saveDirectory, "Apps");
				appsDirectory.mkdir();
				
				//look for files inside the GWrap jar
				URL jar = src.getLocation();
System.out.println("src jar "+jar.getFile());				
				ZipInputStream zip = new ZipInputStream( jar.openStream());
				BufferedInputStream bis = new BufferedInputStream (zip);
				ZipEntry ze = null;
				FileOutputStream fos = null;
				DataOutputStream dos = null;
				int count;
				byte data[] = new byte[2048];
				
				//for each zip entry in the jar file
				while( ( ze = zip.getNextEntry() ) != null ) {
					//is it a doc or app?
					String entryName = ze.getName();				
					if (entryName.startsWith("Documentation/")){
						String trimmedName = entryName.replace("Documentation/", "");
						if (trimmedName.length() == 0) continue;
						fos = new FileOutputStream(new File (docsDirectory, trimmedName));
					}
					else if (entryName.startsWith("Apps/")){
						String trimmedName = entryName.replace("Apps/", "");
						if (trimmedName.length() == 0) continue;
						fos = new FileOutputStream(new File (appsDirectory, trimmedName));
					}
					else fos = null;
					//write it out?
					if (fos != null){
						dos = new DataOutputStream( new BufferedOutputStream (fos));
						while ((count = bis.read(data, 0, 2048))!= -1)  dos.write(data, 0, count);
						dos.close();
						fos.close();
					}
				}
				//close zip stream
				bis.close();
				zip.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public File getLastFile() {
		return lastFile;
	}

	public void setLastFile(File lastFile) {
		this.lastFile = lastFile;
	}   

	public void throwError(String text) {
		JOptionPane.showMessageDialog(this, text, null, JOptionPane.ERROR_MESSAGE);
	}

	/**Launches a command for a given application.*/
	protected void runCommand ()
	{     
		ArrayList<String> cmdList=new ArrayList<String>();      // Command to be run
		ArrayList<String> cmdDisplay=new ArrayList<String>();   // Command to be displayed
		ArrayList<String> paramsList=new ArrayList<String>();   // Active parameter list
		String javaPath = prefsDialog.getProperty("javaPath");
		cmdList.add(javaPath);
		cmdDisplay.add(javaPath);
		String memoryParam = prefsDialog.getProperty("memoryParams");
		if(memoryParam.length() > 0) {
			cmdList.add(memoryParam);
			cmdDisplay.add(memoryParam);

		}
		cmdList.add("-jar");
		cmdList.add(jarFile.toString());
		cmdDisplay.add("-jar " + jarFile.toString());

		paramsList.add("$<$");

		// Add parameters that have data in text field
		for(int i = 0; i < options.length; i++) {
			Option thisOption = options[i];
			JTextArea uInput = thisOption.getUserInputS();
			String cmdText = uInput.getText().trim();
			paramsList.add(thisOption.getFlag() + "$" + cmdText);
			if(cmdText.length()==0) {
				// Nothing to do if no text
				continue;
			}
			cmdList.add(thisOption.getFlag()); 
			if(cmdText.compareTo("Y")==0 || cmdText.compareTo("y")==0) {
				cmdDisplay.add(thisOption.getFlag());
				continue;
			} else {
				cmdList.add(cmdText);           
				cmdDisplay.add(thisOption.getFlag() + " " + cmdText);
			}        
		}
		paramsList.add("$>$");

		String str [] = new String [cmdList.size ()];
		cmdList.toArray (str);

		String dispStr [] = new String [cmdDisplay.size ()];
		cmdDisplay.toArray (dispStr);

		String paramsStr [] = new String [paramsList.size ()];
		paramsList.toArray (paramsStr);

		String[] helpCmd = {java, "-jar", jarFile.toString()};

		CommandThread newCommand = new CommandThread(this);
		newCommand.setFileName(currentCommand);
		newCommand.setRunStartTime(new java.text.SimpleDateFormat("MM/dd/yyyy hh:mm:ss a").format(new java.util.Date()));
		newCommand.setCommandThreadID(new java.text.SimpleDateFormat("MMddyyhhmmssSSS").format(new java.util.Date()));
		newCommand.setCommand(str);
		newCommand.setCommandDisplay(dispStr);
		newCommand.setSaveParams(paramsStr);
		newCommand.setHelpCommand(helpCmd);
		newCommand.start();

		jobs.put(newCommand.getCommandThreadID(), newCommand);
		jobsDialog.refreshJobsList();
	}

	public void setParams(String[] paramVals) {
		for(int i = 0; i < options.length; i++) {
			Option thisOption = options[i];
			JTextArea uInput = thisOption.getUserInputS();
			uInput.setText(paramVals[i]);
		}
	}

	public void clearParams() {
		for(int i = 0; i < options.length; i++) {
			Option thisOption = options[i];
			JTextArea uInput = thisOption.getUserInputS();
			uInput.setText("");
		}
	}

	private void showAppInstructions() {
		Container pane = getContentPane();
		TableLayout layout = (TableLayout) pane.getLayout();
		layout.deleteRow (2);
		layout.insertRow (2, 600);
		//splash screen with info and instructions
		InstuctionsPanel appPanel = new InstuctionsPanel();
		JScrollPane applicationScrollPane;
		applicationScrollPane = new JScrollPane(appPanel);
		this.addComponentListener(appPanel);
		pane.add (applicationScrollPane, "0, 2");
		setVisible(false);  
		setSize (799, 699);
		setSize (800, 700);
		setVisible(true);
	}       

	/**Called when an app is selected from the menu*/
	public void processJarFile(String fileName, boolean initParams) {

		jarFile = new File(appsDirectory, fileName);
		currentCommand = fileName;

		if (jarFile.canRead() == false) {
			Swing.throwError(null, "Aborting: cannot find java jar application -> "+fileName);
			System.exit(0);
		}

		//launch jar application and fetch output
		String[] command = {java, "-jar", jarFile.toString()};
		menu = IO.executeCommandLineReturnAll(command);

		//Misc.printArray(menu);

		//check menu
		checkMenu(command);

		//parse menu
		parseMenu();

		//check parsing
		checkMenuParsing();

		//print to debug
		if (debug){
			System.out.println("\n"+ appDescription+ "\n");
			Misc.printArray(options);
			System.out.println("\n");
			Misc.printArray(menu);
		}

		//make swing bits
		makeSwingComponents();

		Container pane = getContentPane();
		TableLayout layout = (TableLayout) pane.getLayout();
		layout.deleteRow (2);
		layout.insertRow (2, TableLayout.PREFERRED);
		ApplicationPanel appPanel = new ApplicationPanel(this);


		addComponentListener(appPanel);

		pane.add(appPanel, "0,2");
		setVisible(false);  
		pack();
		Dimension dim = appPanel.getPreferredSize();
		dim.height = dim.height+90; /////////////////////Fix this!
		setSize(dim);

		setVisible(true);
		toolPanel.enableActionButtons();
		resultsDialog.refreshResultsWindow(currentCommand);
		if (initParams) {
			resultsDialog.setClickedParameters(0);
		}
	}


	/**
	 * Adds menus to top of frame
	 */
	public void addMenu (Container pane) {
		JMenuBar menuBar = new JMenuBar();
		appsMenuAL = new JMenu("Apps (A-L)");
		appsMenuAL.addMenuListener(new MainMenuListener());
		menuBar.add (appsMenuAL);
		appsMenuMZ = new JMenu("Apps (M-Z)");
		appsMenuMZ.addMenuListener(new MainMenuListener());
		menuBar.add (appsMenuMZ);
		setAppsMenusItems();
		prefsMenu = new JMenu("Preferences");
		prefsMenu.addMenuListener(new MainMenuListener());
		menuBar.add (prefsMenu);
		JMenuItem item = new JMenuItem("View/Edit Preferences");
		item.addActionListener(new PrefsMenuItemListener());
		prefsMenu.add (item);          
		resultsMenu = new JMenu("Results");
		resultsMenu.addMenuListener(new MainMenuListener());
		menuBar.add (resultsMenu);
		item = new JMenuItem("Show Results");
		item.addActionListener(new ResultsMenuItemListener());
		resultsMenu.add (item);          
		jobsMenu = new JMenu("Jobs");
		jobsMenu.addMenuListener(new MainMenuListener());
		menuBar.add (jobsMenu);
		item = new JMenuItem("Show Job Queue");
		item.addActionListener(new JobsMenuItemListener());
		jobsMenu.add (item);          
		helpMenu = new JMenu("Help");
		helpMenu.addMenuListener(new MainMenuListener());
		menuBar.add (helpMenu);
		setHelpMenuItems();
		pane.add (menuBar, "0, 0");
	}

	public void setAppsMenusItems() {
System.out.println("Adding Apps menu items "+appsDirectory);
		appsMenuAL.removeAll();
		appsMenuMZ.removeAll();
		if(appsDirectory.exists()) {
			
			File[] listOfFiles = appsDirectory.listFiles();
			ArrayList<String> appList=new ArrayList<String>(); 
			for (int i = 0; i < listOfFiles.length; i++) {
				if (listOfFiles[i].isFile()) {
					appList.add(listOfFiles[i].getName());
				}
			}
			// Convert to String [] so list can be sorted alphabetically
			String str [] = new String [appList.size ()];
			appList.toArray (str);
			Arrays.sort(str);
			Misc.printArray(str);
			//for each app
			for(int i = 0; i < str.length; i++) {
				String itemName = str[i];
				JMenuItem item = new JMenuItem(itemName);
				item.addActionListener(new AppsMenuItemListener());
				if(itemName.toLowerCase().compareTo("m") < 0) {
					appsMenuAL.add (item); 
				} else {
					appsMenuMZ.add (item);                      
				}
			}         
		}    
	}    


	public void setHelpMenuItems() {
		helpMenu.removeAll();
		if(docsDirectory.exists()) {
			File[] listOfFiles = docsDirectory.listFiles();
			ArrayList<String> docList=new ArrayList<String>(); 
			for (int i = 0; i < listOfFiles.length; i++) {
				if (listOfFiles[i].isFile()) {
					String fileName = listOfFiles[i].getName();
					int htmlIndx = fileName.indexOf(".html");
					if(htmlIndx > 0) {
						String menuLabel = fileName.substring(0, htmlIndx);
						menuLabel = menuLabel.substring(0, 1).toUpperCase() + menuLabel.substring(1);
						helpMenuFiles.put(menuLabel, fileName);
						docList.add(menuLabel);              
					}
				}
			}

			// Convert to String [] so list can be sorted alphabetically
			String str [] = new String [docList.size ()];
			docList.toArray (str);
			Arrays.sort(str);
			for(int i = 0; i < str.length; i++) {
				JMenuItem item = new JMenuItem(str[i]);
				item.addActionListener(new HelpMenuItemListener());
				helpMenu.add (item);          
			} 
			helpMenu.addSeparator();
			JMenuItem item = new JMenuItem("About");
			item.addActionListener(new HelpMenuItemListener());
			helpMenu.add(item);
		}    
	}

	/**Will load preferences or set defaults in PrefsDialog*/
	public void showPrefsDialog() {
		prefsDialog.initTextFields();
		prefsDialog.setLocationRelativeTo(bw);
		prefsDialog.setBounds(20, 40, 600, 250);
		prefsDialog.setVisible(true);     
	}

	private class AppsMenuItemListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			String actionCommand = e.getActionCommand();
			processJarFile(actionCommand, true);
		}
	}

	private class PrefsMenuItemListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			showPrefsDialog();
		}
	}    

	private class ResultsMenuItemListener implements ActionListener{
		public void actionPerformed(ActionEvent e){
			resultsDialog.setLocationRelativeTo(bw);
			resultsDialog.refreshResultsWindow(currentCommand);           
			resultsDialog.setVisible(true);
		}
	}    

	private class JobsMenuItemListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			jobsDialog.setLocationRelativeTo(bw);
			jobsDialog.refreshJobsList();
			jobsDialog.setVisible(true);
		}
	}    

	private class HelpMenuItemListener implements ActionListener
	{

		public void actionPerformed(ActionEvent e)
		{
			String actionCommand = e.getActionCommand();

			if(actionCommand.compareTo("About")==0) {
				AboutWindow aw = new AboutWindow();
				aw.setLocationRelativeTo(bw);
				aw.setBounds(20, 40, 800, 600);
				aw.setVisible(true);   
				return;
			}

			String url = docsDirectory+"/"+helpMenuFiles.get(actionCommand);

			boolean browserLaunched = false;
			if (Desktop.isDesktopSupported()) {
				Desktop desktop = Desktop.getDesktop();
				URI uri = null;
				File helpFile = new File(url);
				if (desktop.isSupported(Desktop.Action.BROWSE)) {
					browserLaunched = true;
					try {
						uri = helpFile.toURI();
						desktop.browse(uri);
					} catch(IOException ioe) {
						System.out.println("The system cannot find the " + uri + " file specified");
						//ioe.printStackTrace();
					} 
				}
			}
			if(!browserLaunched) {
				HelpWindow hw = new HelpWindow(url);
				hw.setLocationRelativeTo(bw);
				hw.setBounds(20, 40, 800, 600);
				hw.setVisible(true);              
			}
		}
	}



	private class MainMenuListener implements MenuListener
	{
		public void menuSelected(MenuEvent e) {
			//JMenu menu = (JMenu) e.getSource();
			//String menuText = menu.getText();
		}

		public void menuCanceled(MenuEvent e) {
		}

		public void menuDeselected(MenuEvent e) {
		}

	}



	public void addCommandButtons (Container pane, TableLayout layout)
	{
		JPanel buttonPanel = new JPanel();
		pane.add (buttonPanel, "1, 4, 5, 4");

		for (int i = 1; i <= 5; i++)
		{
			JButton button = new JButton("Button " + i);
			buttonPanel.add (button);
		}
	}



	public void allowClosing ()
	{
		addWindowListener
		(new WindowAdapter()
		{
			public void windowClosing (WindowEvent e)
			{
				System.exit (0);
			}
		}
		);
	}

	public void makeSwingComponents(){
		appName = new JTextArea(jarFile.getName());
		appName.setEditable(false);
		appName.setFont(ApplicationPanel.FONT_BOLD);
		appName.setBackground(Color.WHITE);
		appName.setForeground(Color.BLUE);

		appDescriptionS = new JTextArea(appDescription);
		appDescriptionS.setEditable(false);
		appDescriptionS.setWrapStyleWord(true);
		appDescriptionS.setLineWrap(true);
		appDescriptionS.setFont(ApplicationPanel.FONT_BOLD);
		appDescriptionS.setBackground(Color.WHITE);

		for (int i=0; i< options.length; i++) options[i].makeSwingComponents(i);
	}

	public void checkMenuParsing(){
		if (appDescription == null || options == null || options.length == 0){
			Swing.throwError(null, "Problem parsing comand line menu? No appDescription? No options?\n"+Misc.stringArrayToString(menu, "\n"));
			System.exit(0);
		}
	}
	public void checkMenu(String[] command){
		//executed?
		if (menu == null){
			Swing.throwError(null, "Problem executing java command -> '"+Misc.stringArrayToString(command, " ")+ "'");
			System.exit(0);
		}
		//jar file OK?
		if (menu[0].contains("corrupt")){
			Swing.throwError(null, "Problem executing java jar command. "+menu[0]);
			System.exit(0);
		}
	}


	public void parseMenu(){
		StringBuilder descriptionSB = new StringBuilder();

		ArrayList<Option> parametersAL = new ArrayList<Option>();
		boolean pastDescription = false;
		for (int i=0; i < menu.length; i++){
			String line = menu[i].trim();
			//empty?
			if (line.length() == 0) continue;
			//at stop?
			if (line.startsWith("Example")) break;
			//skip?
			if (line.startsWith("Option") || line.startsWith("Param") || line.startsWith("*")) continue;
			//parameter?
			Matcher mat = flag.matcher(line);
			if (mat.find()) {
				pastDescription = true;
				while (i < menu.length){
					i++;
					//look at next line, at stop?
					if (i == menu.length) {
						parametersAL.add(new Option(line, this));
						break;
					}
					else {
						//does the next line start with a flag
						String following = menu[i];
						mat = flag.matcher(following);
						if (mat.find()) {
							parametersAL.add(new Option(line, this));
							i--;
							break;
						}

						//does the next line start with white space
						mat = startSpace.matcher(following);
						if (mat.find()){
							//and is not empty
							if (following.trim().length() !=0) line = line +"\n"+following;
							else {
								parametersAL.add(new Option(line, this));
								break;
							}
						}
						else {
							i--;
							parametersAL.add(new Option(line, this));
							break;
						}
					}
				}
			}
			//add to appDescription?
			if (pastDescription) continue;
			else {
				descriptionSB.append(line);
				descriptionSB.append(" ");
			}
		}
		//set fields
		appDescription = descriptionSB.toString();
		options = new Option[parametersAL.size()];
		parametersAL.toArray(options);

	}

	public File getJarFile() {
		return jarFile;
	}
	public void setJarFile(File jarFile) {
		this.jarFile = jarFile;
	}
	public String[] getMenu() {
		return menu;
	}
	public void setMenu(String[] menu) {
		this.menu = menu;
	}
	public String getDescription() {
		return appDescription;
	}
	public void setDescription(String description) {
		this.appDescription = description;
	}
	public Option[] getOptions() {
		return options;
	}
	public void setOptions(Option[] parameters) {
		this.options = parameters;
	}
	public String getJava() {
		return java;
	}
	public void setJava(String java) {
		this.java = java;
	}

	public JTextArea getAppName() {
		return appName;
	}

	public void setAppName(JTextArea appName) {
		this.appName = appName;
	}

	public String getAppDescription() {
		return appDescription;
	}

	public void setAppDescription(String appDescription) {
		this.appDescription = appDescription;
	}

	public JTextArea getAppDescriptionS() {
		return appDescriptionS;
	}

	public void setAppDescriptionS(JTextArea appDescriptionS) {
		this.appDescriptionS = appDescriptionS;
	} 

	class JobsComparator implements Comparator<Object>  {  
		public int compare ( Object o1, Object o2 )   {  
			String s1 =  ( String ) o1; 
			String s2 =  ( String ) o2; 
			return s2.compareTo ( s1) ; 
		}  
	}
	public File getSaveDirectory() {
		return saveDirectory;
	}

	public PrefsDialog getPrefsDialog() {
		return prefsDialog;
	}

	public JobsDialog getJobsDialog() {
		return jobsDialog;
	}

	public ResultsDialog getResultsDialog() {
		return resultsDialog;
	}

	public File getAppsDirectory() {
		return appsDirectory;
	}

	public File getDocsDirectory() {
		return docsDirectory;
	}  

}
