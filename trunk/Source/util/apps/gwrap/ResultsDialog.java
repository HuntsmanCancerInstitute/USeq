package util.apps.gwrap;

import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeMap;

import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;


class ResultsDialog extends JDialog implements ActionListener,PropertyChangeListener {

	private static final long serialVersionUID = 1L;

	private JOptionPane optionPane;
	private JTextArea jTxtScroll;;

	private String btnClear = "Clear Saved Results";
	private String btnClose = "Close";
	private GWrap_GUI_ClickMe mainApp;
	private String currentCommand = "?";
	protected TreeMap<Integer, String[]> params;


	/** Creates the reusable dialog. */
	public ResultsDialog(GWrap_GUI_ClickMe aFrame) {
		super();
		mainApp = aFrame; 

		currentCommand = "?";
		setTitle("Results for " + currentCommand);


		jTxtScroll = new JTextArea(6, 20);
		

		JScrollPane sp = new JScrollPane(jTxtScroll);

		MouseListener mouseListener = new MouseAdapter() {
			public void mouseClicked(MouseEvent e) {
				if (e.getClickCount() == 2) {
					int position = jTxtScroll.viewToModel( e.getPoint() );
					int index = jTxtScroll.getDocument().getDefaultRootElement().getElementIndex(position)+1;
					setClickedParameters(index);
				}
			}
		};        
		jTxtScroll.addMouseListener(mouseListener);

		JLabel jp = new JLabel("Double click anywhere on a run to copy its parameters into main window for next run.");        
		jp.setBounds(10, 10, 300, 20);
		Rectangle bounds = new Rectangle(10, 20, 300, 300);
		sp.setBounds(bounds);

		this.getContentPane().setLayout(null);
		bounds = new Rectangle(20, 20, 600, 600);
		this.setBounds(bounds);

		//Create an array of the text and components to be displayed.
		Object[] array = {jp, sp};

		//Create an array specifying the number of dialog buttons
		//and their text.
		Object[] options = {btnClear, btnClose};

		//Create the JOptionPane.
		optionPane = new JOptionPane(array,
				JOptionPane.PLAIN_MESSAGE,
				JOptionPane.YES_NO_OPTION,
				null,
				options,
				options[0]);

		//Make this dialog display it.
		setContentPane(optionPane);


		//Handle window closing correctly.
		setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent we) {
				/*
				 * Instead of directly closing the window,
				 * we're going to change the JOptionPane's
				 * value property.
				 */
				optionPane.setValue(new Integer(
						JOptionPane.CLOSED_OPTION));
			}
		});

		//Register an event handler that reacts to option pane state changes.
		optionPane.addPropertyChangeListener(this);
	}

	public void setClickedParameters(int index) {
		Iterator<Integer> it = params.keySet().iterator();
		Integer lastKey = null;
		Integer thisKey = null;
		String[] paramVals = null;
		boolean paramsSet = false;
		while (it.hasNext())
		{
			thisKey = (Integer) it.next();
			if (index < thisKey.intValue()) {
				if(lastKey != null) {
					paramVals = params.get(lastKey);                         
				} else {
					paramVals = params.get(thisKey);                          
				}
				mainApp.setParams(paramVals);
				paramsSet = true;
				break;
			}
			lastKey = new Integer(thisKey);
		}   
		if(thisKey != null && !paramsSet) {
			paramVals = params.get(thisKey);
			mainApp.setParams(paramVals);                   
		}      
	}

	public void refreshResultsWindow(String currCommand) {
		FileReader fr = null;
		BufferedReader br = null;

		jTxtScroll.setText("");
		params = new TreeMap<Integer, String[]>();
		ArrayList<String> paramsList=new ArrayList<String>(); 

		if(currCommand != null) {
			currentCommand = currCommand;
		}
		setTitle("Results for " + currentCommand);

		if(currentCommand.compareTo("?") != 0) {
			File file = new File(mainApp.getPrefsDialog().getHistoryDirectory(), currentCommand+".txt");
			boolean paramsOn = false;
			if(file.exists()) {
				try {
					fr = new FileReader(file);
					br = new BufferedReader(fr);
					String line = null;
					int paramLine = 0;
					int currentLine = 0;
					while ((line=br.readLine()) != null) {
						if(line.compareTo("$<$")== 0) {
							paramsOn = true;
							paramLine = currentLine;
							paramsList=new ArrayList<String>();
							continue;
						}
						if(line.compareTo("$>$")== 0) {
							paramsOn = false;
							String str [] = new String [paramsList.size ()];
							paramsList.toArray (str);
							params.put(new Integer(paramLine), str);
							continue;
						}
						if(!paramsOn) {
							jTxtScroll.append(line+"\n");
							currentLine++;
						} else {
							line=line.substring(line.indexOf("$")+1);
							paramsList.add(line);
						}

					}
					paramsList = null;
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				} catch (IOException e) {
					e.printStackTrace();
				}
				finally {
					if(br != null) {
						try {
							br.close();
							fr.close();
						} catch (IOException e) {
							e.printStackTrace();
						}
					}

				}
			} else {
				jTxtScroll.append("The command has not yet been run.");
			}

		} else {      
			jTxtScroll.append("Please first select a command.");
		}
		jTxtScroll.setCaretPosition(0);
	}    

	/** This method handles events for the text field. */
	public void actionPerformed(ActionEvent e) {
		optionPane.setValue(btnClear);
	}


	/** This method reacts to state changes in the option pane. */
	public void propertyChange(PropertyChangeEvent e) {
		String prop = e.getPropertyName();

		if (isVisible()
				&& (e.getSource() == optionPane)
				&& (JOptionPane.VALUE_PROPERTY.equals(prop) ||
						JOptionPane.INPUT_VALUE_PROPERTY.equals(prop))) {
			Object value = optionPane.getValue();

			if (value == JOptionPane.UNINITIALIZED_VALUE) {
				//ignore reset
				return;
			}

			//Reset the JOptionPane's value.
			//If you don't do this, then if the user
			//presses the same button next time, no
			//property change event will be fired.
			optionPane.setValue(
					JOptionPane.UNINITIALIZED_VALUE);

			if (btnClear.equals(value)) {      
				File file = new File(mainApp.getPrefsDialog().getHistoryDirectory(), currentCommand+".txt");
				if(file.exists()) {
					file.delete();
					refreshResultsWindow(currentCommand);
				}
			} else { //user closed dialog or clicked cancel
				clearAndHide();
			}
		}
	}

	/** This method clears the dialog and hides it. */
	public void clearAndHide() {
		//textField.setText(null);
		setVisible(false);
	}

}

