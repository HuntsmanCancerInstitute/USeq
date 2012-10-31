package util.gen;

/**
 * Emails messages with and without attachments, handles authentication too.
 * 
 * Adapted from http://www.javacommerce.com/displaypage.jsp?text=javamail.sql&id=18274 Author : Sudhir Ancha
*/

import javax.mail.*;
import javax.mail.internet.*;
import javax.activation.*;
import java.util.*;

public class Email{
	
	private String SMTP_HOST_NAME = "myserver.smtphost.com";
	private String SMTP_AUTH_USER = "myusername";
	private String SMTP_AUTH_PWD  = "mypwd";
	private String contentType = "text/plain"; //ie "text/html" or "text/plain"
	
	
	public Email(String smtpHostName, String smtpUser, String smtpPassword, String MIMEContentType){
		SMTP_HOST_NAME = smtpHostName;
		SMTP_AUTH_USER = smtpUser;
		SMTP_AUTH_PWD = smtpPassword;
		if (Misc.isNotEmpty(MIMEContentType)) contentType = MIMEContentType;
	}
	
	public boolean postMail( String recipients[ ], String subject, String message , String from) {
		try{
			//Set the host smtp address
			Properties props = new Properties();
			props.put("mail.smtp.host", SMTP_HOST_NAME);
			props.put("mail.smtp.auth", "true");
			props.put("mail.smtp.starttls.enable","true");
			props.put("mail.transport.protocol", "smtp");
			
			Authenticator auth = new SMTPAuthenticator();
			Session session;
			if (System.getSecurityManager() == null) session = Session.getInstance(props,auth);
			else session = Session.getDefaultInstance(props, auth);
			
			//change to true and then print messages using a stream
			session.setDebug(false);

			// create a message
			Message msg = new MimeMessage(session);
			
			// set the from and to address
			InternetAddress addressFrom = new InternetAddress(from);
			msg.setFrom(addressFrom);
			
			InternetAddress[] addressTo = new InternetAddress[recipients.length];
			for (int i = 0; i < recipients.length; i++)
			{
				addressTo[i] = new InternetAddress(recipients[i]);
			}
			msg.setRecipients(Message.RecipientType.TO, addressTo);
			
			
			// Setting the Subject and Content Type
			msg.setSubject(subject);
			msg.setContent(message, contentType);
			Transport.send(msg);
			return true;
		}catch (MessagingException e){
			System.err.println("Prob with postMail()");
			e.printStackTrace();
		}
		return false;
	}
	
	public boolean postWithAttachment( String recipients[ ], String subject, String textMessage, 
			String fullPathFileName, String attachmentFileName, String from) {
		try{
			//Set the host smtp address
			Properties props = new Properties();
			props.put("mail.smtp.host", SMTP_HOST_NAME);
			props.put("mail.smtp.auth", "true");
			
			Authenticator auth = new SMTPAuthenticator();
			Session session;
			if (System.getSecurityManager() == null) session = Session.getInstance(props,auth);
			else session = Session.getDefaultInstance(props, auth);
			
			//change to true and then print messages using a stream
			session.setDebug(false);

			//create a message
			Message msg = new MimeMessage(session);
			
			//set the from and to address
			InternetAddress addressFrom = new InternetAddress(from);
			msg.setFrom(addressFrom);
			
			InternetAddress[] addressTo = new InternetAddress[recipients.length];
			for (int i = 0; i < recipients.length; i++)
			{
				addressTo[i] = new InternetAddress(recipients[i]);
			}
			msg.setRecipients(Message.RecipientType.TO, addressTo);
			
			//Setting the Subject
			msg.setSubject(subject);
			
			//make body part, text message to show in email
			MimeBodyPart mbp1 = new MimeBodyPart();
			mbp1.setText(textMessage);
	   
	        //create the second body part, the email attachment
	        MimeBodyPart mbp2 = new MimeBodyPart();
	        DataSource source = new FileDataSource(fullPathFileName);
	        mbp2.setDataHandler(new DataHandler(source));
	        mbp2.setFileName(attachmentFileName);
	        
	        //create the Multipart and its parts to it
	        Multipart mp = new MimeMultipart();
	        mp.addBodyPart(mbp1);
	        mp.addBodyPart(mbp2);
	   
	        //add the Multipart to the message
	        msg.setContent(mp);
			//send
			Transport.send(msg);
			return true;
		}catch (MessagingException e){
			System.err.println("Prob with postWithAttachment()");
			e.printStackTrace();
		}
		return false;
	}
	
	/**Simple static method for sending email, no authentication, modified from http://www.javacommerce.com/articles/sendingmail.htm Sudhir Ancha.*/
	public static boolean postMailNoAuthentication(String[] recipients,  String subject, String message , String from, 
			String MIMEType, String SMTPAddress){
		try{
			// add handlers for main MIME types
			MailcapCommandMap mc = (MailcapCommandMap)CommandMap.getDefaultCommandMap();
			mc.addMailcap("text/html;; x-java-content-handler=com.sun.mail.handlers.text_html");
			mc.addMailcap("text/xml;; x-java-content-handler=com.sun.mail.handlers.text_xml");
			mc.addMailcap("text/plain;; x-java-content-handler=com.sun.mail.handlers.text_plain");
			mc.addMailcap("multipart/*;; x-java-content-handler=com.sun.mail.handlers.multipart_mixed");
			mc.addMailcap("message/rfc822;; x-java-content-handler=com.sun.mail.handlers.message_rfc822");
			CommandMap.setDefaultCommandMap(mc);
			
			//Set the host smtp address
			Properties props = new Properties();
			props.put("mail.smtp.host", SMTPAddress);//"mail.smtp.host", "smtp.jcom.net");
			
			// create some properties and get  the default Session
			Session session = Session.getDefaultInstance(props,  null);
			session.setDebug(false);
			
			// create a message
			Message msg = new MimeMessage(session);
			// set the from and to address
			InternetAddress addressFrom = new InternetAddress(from);
			msg.setFrom(addressFrom);
			
			InternetAddress[] addressTo = new  InternetAddress[recipients.length];
			for (int i = 0; i < recipients.length; i++){
				addressTo[i] = new  InternetAddress(recipients[i]);
			}
			msg.setRecipients(Message.RecipientType.TO, addressTo);
			
			
			//  Optional : You can also set your custom headers in the Email if you Want
			//msg.addHeader("xxx", "yyy");
			
			//  Setting the Subject and Content Type ie "text/html" or "text/plain"
			msg.setSubject(subject);
			msg.setContent(message, MIMEType);
			Transport.send(msg);
			return true;
		}catch (Exception e){
			System.err.println("Prob with postMailNoAuthentication()");
			e.printStackTrace();
		}
		return false;
	}	
	
	
	/**
	 * SimpleAuthenticator is used to do simple authentication
	 * when the SMTP server requires it.
	 */
	private class SMTPAuthenticator extends javax.mail.Authenticator
	{
		
		public PasswordAuthentication getPasswordAuthentication()
		{
			String username = SMTP_AUTH_USER;
			String password = SMTP_AUTH_PWD;
			return new PasswordAuthentication(username, password);
		}
	}
	
}


