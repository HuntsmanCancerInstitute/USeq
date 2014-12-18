package util.gen;

import java.util.Properties;
import javax.mail.Message;
import javax.mail.Session;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeMessage;

/**
 * @author Crunchify.com with mods by nix
 * Insecure!  You'll need to modify your google account to allow the messages to be sent
 *
 */
public class Gmail {

	static Properties mailServerProperties;
	static Session getMailSession;
	static MimeMessage generateMailMessage;
	
	public static void main (String[] args){
		Gmail.generateAndSendEmail(new String[]{"david.nix@hci.utah.edu"}, "Test2", "Hey Dave", "david.austin.nix@gmail.com", "");
	}

	public static boolean generateAndSendEmail(String[] recipients, String subject, String message, String googleEmailAccount, String password)  {
		try {
			//Step1
			//System.out.println("\n 1st ===> setup Mail Server Properties..");
			mailServerProperties = System.getProperties();
			mailServerProperties.put("mail.smtp.port", "587");
			mailServerProperties.put("mail.smtp.auth", "true");
			mailServerProperties.put("mail.smtp.starttls.enable", "true");
			//System.out.println("Mail Server Properties have been setup successfully..");

			//Step2
			//System.out.println("\n\n 2nd ===> get Mail Session..");
			getMailSession = Session.getDefaultInstance(mailServerProperties, null);
			generateMailMessage = new MimeMessage(getMailSession);
			for (String add: recipients){
				generateMailMessage.addRecipient(Message.RecipientType.TO, new InternetAddress(add));
			}
			//generateMailMessage.addRecipient(Message.RecipientType.CC, new InternetAddress("test2@crunchify.com"));
			generateMailMessage.setSubject(subject);
			generateMailMessage.setContent(message, "text/html");
			//System.out.println("Mail Session has been created successfully..");

			//Step3
			//System.out.println("\n\n 3rd ===> Get Session and Send mail");  
			Transport transport = getMailSession.getTransport("smtp");
			// Enter your correct gmail UserID and Password (XXXarpitshah@gmail.com)
			transport.connect("smtp.gmail.com", googleEmailAccount, password);
			transport.sendMessage(generateMailMessage, generateMailMessage.getAllRecipients());
			transport.close();
			return true;
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}
	}

}
