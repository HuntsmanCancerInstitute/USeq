package edu.utah.ames.bioinfo;

import javax.mail.internet.AddressException;
import javax.mail.internet.InternetAddress;

/*
 * This class holds an email address for use by methods of other classes
 */

/**
 *
 * @author darren.ames@hci.utah.edu
 */
public class EmailAddress {
    
    //fields
    private String address;
    
    //constructor
    public EmailAddress(String email) {
       address = email;
    }
    
    //email address validation
    public boolean isValidEmailAddress(String em) {
        boolean result = true;
        try {
            InternetAddress emailAddr = new InternetAddress(em);
            emailAddr.validate();
        }
        catch (AddressException ex) {
            result = false;
        }
        return result;
    }
    
    //sets email address field
    public void setEmailAddress(String email) {
        address = email;
    }
    
    //returns the email address from address field
    public String getEmailAddress() {
        return address;
    }
    
}
