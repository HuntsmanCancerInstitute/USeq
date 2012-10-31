package util.apps.gwrap;

import java.awt.Color;
import java.awt.Container;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import javax.swing.JButton;
import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;

class HelpWindow extends JFrame {
    /**
   * 
   */
  private static final long serialVersionUID = 1L;

    private Container pane;
    private JEditorPane dispText;
    
    private String closeActionCommand = "cancel";
    

    /** Creates the reusable dialog. */
    public HelpWindow(String helpFilePath) {
         
        setTitle("Help");
        pane = getContentPane();
        pane.setLayout(null);
        pane.setBounds(0, 0, 800, 600);
        
        dispText = new JEditorPane();
        dispText.setBounds(0, 0, 792, 530);
        
        String helpFileContent = readFile(helpFilePath);
        dispText.setContentType("text/html");
        dispText.setText(helpFileContent);
        dispText.setCaretPosition(0);
        dispText.setEditable(false);
        dispText.addHyperlinkListener(new HyperlinkListener() {
          public void hyperlinkUpdate(HyperlinkEvent e) {
            if (e.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
              // If we arrive here that means we weren't able to launch the help files
              // from the browser in the first place, so just indicate that hyperlinks won't work
              JOptionPane.showMessageDialog(null, "Sorry. The application was unable to interact with the destop to launch the URL. Desktop functionality may not be enabled on this computer.");           
            }
          }
        });        

        JScrollPane sp = new JScrollPane(dispText);
        sp.setBounds(0, 0, 792, 530);
        
        add(sp);
        
        HelpButtonListener prefsButtonListener = new HelpButtonListener();
        
        
        
        JButton closeBtn = new JButton("Close");
        closeBtn.setActionCommand(closeActionCommand);
        closeBtn.setBounds(360, 540, 80, 20);
        closeBtn.addActionListener(prefsButtonListener);
        add(closeBtn);
    }
    
    private String readFile(String fileName) {
      
      File file = new File(fileName);
      
      char[] buffer = null;
      
      try {
          BufferedReader bufferedReader = new BufferedReader(new FileReader(file));

          buffer = new char[(int)file.length()];

          int i = 0;
          int c = bufferedReader.read();

          while (c != -1) {
              buffer[i++] = (char)c;
              c = bufferedReader.read();
          }
      } catch (FileNotFoundException e) {
        e.printStackTrace();
      } catch (IOException e) {
        e.printStackTrace();
      }
      return new String(buffer);
    }
          
    private class HelpButtonListener implements ActionListener {
       public void actionPerformed(ActionEvent e) {
         if(e.getActionCommand().compareTo(closeActionCommand)==0) {
           clearAndHide(); 
           return;
         }
       }
    }    


    /** This method clears the dialog and hides it. */
    public void clearAndHide() {
        //textField.setText(null);
        setVisible(false);
    }
    
}

