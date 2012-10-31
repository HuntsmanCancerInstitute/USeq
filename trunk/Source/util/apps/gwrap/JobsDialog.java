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
import java.util.ArrayList;
import java.util.Iterator;

import javax.swing.DefaultListModel;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;

/* Used by BioWrapper.java. */
class JobsDialog extends JDialog
implements ActionListener,
PropertyChangeListener {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private JOptionPane optionPane;

	private String btnDelete = "Delete Selected Job";
	private String btnClear = "Kill Selected Job";
	private String btnClose = "Close";
	private GWrap_GUI_ClickMe mainFrame;
	private DefaultListModel listModel;
	private JList jlistScroll; 


	/** Creates the reusable dialog. */
	public JobsDialog(GWrap_GUI_ClickMe aFrame) {
		//super(aFrame, false);
		super();
		mainFrame = aFrame;       

		setTitle("Job Queue");

		listModel = new DefaultListModel();  
		jlistScroll = new JList(listModel);        
		JScrollPane sp = new JScrollPane(jlistScroll);

		MouseListener mouseListener = new MouseAdapter() {
			public void mouseClicked(MouseEvent e) {
				if (e.getClickCount() == 2) {
					String selectedVal = (String) jlistScroll.getSelectedValue();
					if(selectedVal != null) {
						String jobID = selectedVal.substring(selectedVal.indexOf("ID=")+3);
						CommandThread ct = mainFrame.jobs.get(jobID);
						mainFrame.processJarFile(ct.getFileName(), false);
						String [] savedParams = ct.getSaveParams();
						ArrayList<String> paramsList=new ArrayList<String>();
						for(int i = 0; i < savedParams.length; i++) {
							String line = savedParams[i];
							if(line.compareTo("$<$")== 0 || line.compareTo("$>$")== 0) {
								continue;
							} else {
								line=line.substring(line.indexOf("$")+1);
								paramsList.add(line);                      
							}
						}
						String str [] = new String [paramsList.size ()];
						paramsList.toArray (str);
						mainFrame.setParams(str);   
					}           
				}
			}
		};        
		jlistScroll.addMouseListener(mouseListener);        

		JLabel jp = new JLabel("Double click on a job to copy its parameters into main window.");        
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
		Object[] options = {btnClear, btnDelete, btnClose};

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

	public void refreshJobsList() {
		listModel.clear();
		Iterator<String> it = mainFrame.jobs.keySet().iterator();

		while (it.hasNext()) { 
			String commandThreadID = (String) it.next();
			CommandThread ct = mainFrame.jobs.get(commandThreadID);
			listModel.addElement(ct.getFileName() + " " + ct.getRunStartTime() 
					+ " " +  ct.getStatus() + " ID=" + commandThreadID);
		}
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
				String selectedVal = (String) this.jlistScroll.getSelectedValue();
				if(selectedVal != null) {
					String jobID = selectedVal.substring(selectedVal.indexOf("ID=")+3);
					int confirm = JOptionPane.showConfirmDialog(this,"Are you sure you want to kill process ID #" 
							+ jobID+"?","Kill Confirm",JOptionPane.OK_CANCEL_OPTION);
					if (confirm == JOptionPane.YES_OPTION) {
						CommandThread ct = mainFrame.jobs.get(jobID);
						ct.killProcess();
						ct.setStatus("Process Killed");
						refreshJobsList();
					}
				} else {               
					JOptionPane.showMessageDialog(this, "Please select an item from the list.");  
				}

			} else { //user closed dialog or clicked cancel

				if (btnDelete.equals(value)) {              

					String selectedVal = (String) this.jlistScroll.getSelectedValue();
					if(selectedVal != null) {
						String jobID = selectedVal.substring(selectedVal.indexOf("ID=")+3);
						int confirm = JOptionPane.showConfirmDialog(this,"Are you sure you want to delete process ID #" 
								+ jobID+"?","Delete Confirm",JOptionPane.OK_CANCEL_OPTION);
						if (confirm == JOptionPane.YES_OPTION) {
							CommandThread ct = mainFrame.jobs.get(jobID);
							ct.killProcess();
							mainFrame.jobs.remove(jobID);
							refreshJobsList();
						}
					} else {               
						JOptionPane.showMessageDialog(this, "Please select an item from the list.");  
					}

				} else { //user closed dialog or clicked cancel              
					clearAndHide();
				}
			}
		}
	}

	/** This method clears the dialog and hides it. */
	public void clearAndHide() {
		//textField.setText(null);
		setVisible(false);
	}
}

