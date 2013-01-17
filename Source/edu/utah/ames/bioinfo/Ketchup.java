package edu.utah.ames.bioinfo;

import java.io.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.mail.Message;
import javax.mail.MessagingException;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeMessage;

/**
 * This class identifies users in the Tomato jobs directory, creates a User object for
 * each, and then for each user, finds and stores their associated email
 * addresses to be used for notification of impending file deletion. Files are grouped into
 * either warnFiles or deleteFiles based on size and lastModifiedDate. Emails users
 * who have warnFiles and/or deleteFiles that they have files pending deletion.
 * 
 * @author darren.ames@hci.utah.edu
 */

//TODO: add this class and dependencies to jar file

public class Ketchup {

	//fields
	private String path = "/tomato/job/";
	//regex to find email addresses
	static String REGEX = "^#e\\s+(.+@.+)";
	static HashMap<File, Long> fileDates = null;
	//set warn cutoff time: 5 days ago
	private long warnDate = System.currentTimeMillis()
			- (5L * 24 * 60 * 60 * 1000);
	//set delete cutoff time: 7 days ago
	private long deleteDate = System.currentTimeMillis()
			- (7L * 24 * 60 * 60 * 1000);
	private long fileSize = 1048576;
	private String smtpHostName = "mail.inscc.utah.edu";

	//constructor
	public Ketchup (String[] args) {
		
		processArgs(args);
		
		File tomatoDir = new File(path);
		File savedFileDates = new File(tomatoDir, "savedFileDates.ser");
		
		//fetch file dates and save in serialized file
		fileDates = fetchFileDates(savedFileDates);

		//create File object 
		File[] dirs = tomatoDir.listFiles();
		ArrayList<User> users = fetchUsers(dirs);
		
		for (User user : users) {
			
			//find email addresses for each user
			fetchEmailAddress(user);

			//start walking
			walk(user);
			
			//send some emails!
			Ketchup.sendMail(user);
		}
		Stuff.saveObject(savedFileDates, fileDates);
	}
	
	/**
	 * Descends through files in each user's job directory.
	 * @param user
	 * @throws IOException
	 */
	public void walk(User user) {
		
		
		ArrayList<File> alFiles = fetchAllFilesRecursively(user.getDirectory());
		
		for (File f : alFiles) {
			checkFile(user, f);
		}
	}
	
	/**
	 * Fetches files that don't start with a '.' from a directory recursing through sub directories
	 * @param directory
	 * @return
	 */
	public static ArrayList<File> fetchAllFilesRecursively (File directory){
		ArrayList<File> files = new ArrayList<File>(); 
		File[] list = directory.listFiles();
		if (list != null){
			for (int i = 0; i < list.length; i++){
				if (list[i].isDirectory()) {
					ArrayList<File> al = fetchAllFilesRecursively (list[i]);
					int size = al.size();
					for (int x = 0; x < size; x++){
						File test = al.get(x);
						if (test.getName().startsWith(".") == false) files.add(test);
					}
				}
				else if (list[i].getName().startsWith(".") == false) files.add(list[i]);				
			}
		}
		return files;
	}
	
	/**
	 * Compares file last modified dates to see if action is necessary.
	 * @param user
	 * @param fileObject
	 * @throws IOException
	 */
	public void checkFile(User user, File fileObject) {

		try {
			//is the file too big?
			if (fileObject.length() > fileSize) {
				//Long fileDate = fileDates.get(f);
				Long fileDate = Long.valueOf(fileObject.lastModified());

				//check for fileDate, if null, set to current time
				if (fileDate == null) {
					fileDate = System.currentTimeMillis();
					fileDates.put(fileObject, fileDate);
				} else {
					//warn cutoff? >=5 days ago && < 7 days ago
					if ((fileDate.longValue() <= warnDate)
							&& (fileDate.longValue() > deleteDate))
						//add fileObject to warnFiles ArrayList
						user.getWarnFiles().add(fileObject);
					else
						//delete cutoff >=7 days ago
						if (fileDate.longValue() <= deleteDate)
							//add fileObject to deleteFiles ArrayList
							user.getDeleteFiles().add(fileObject);
					//****uncomment this to start deleting offending files****
					//fileObject.deleteOnExit();
				}
			}
		}
		catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	/**
	 * Lists directory names, which will be used as user names.
	 * @param files
	 * @return
	 */
	public static ArrayList<User> fetchUsers(File[] files) {

		//create ArrayList to hold User names
		ArrayList<User> userNames = new ArrayList<User>();

		//create User object to hold User names
		for (File file : files) {
			if (file.isDirectory()) {
				User user = new User(file);

				//add user names to user objects
				userNames.add(user);
			}
		}
		return userNames;
	}

	/**
	 * Grabs email addresses from cmd.txt files found in user directories.
	 * @param user
	 * @throws IOException
	 */
	public static void fetchEmailAddress(User user) {

		try {
			//call pattern to find valid email addresses within a file
			Pattern pattern = Pattern.compile(REGEX, Pattern.CASE_INSENSITIVE);

			//find the cmd.txt file
			FileFind ff = new FileFind(user.getDirectory(), user);
			ff.find();

			//read in cmd.txt file
			user.getCmdFiles().get(0);

			for (File c : user.getCmdFiles()) {
				BufferedReader br = new BufferedReader(new FileReader(c));

				//search for email address string in file
				for (String line = br.readLine(); line != null; line = br
						.readLine()) {
					Matcher m = pattern.matcher(line);
					//System.out.println("test");
					if (m.find()) {

						//set email address in User HashSet
						user.getEmail().add(m.group(1));
					}
				}
				br.close();
			}
		}
		catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	/**
	 * Grabs file modified dates from each user's directories and saves in a serialized File object.
	 * @param serFile
	 * @return
	 */
	public static HashMap<File, Long> fetchFileDates(File serFile) {

		if (serFile.canRead() == false || !serFile.exists()) {
			return new HashMap<File, Long>();
		}
		return (HashMap<File, Long>) Stuff.fetchObject(serFile);
	}

	/**
	 * method to compose email body
	 * @param msg
	 * @param user
	 * @return
	 * @throws IOException
	 */
	public static String getMessage(String msg, User user) throws IOException {
		
		//instantiate new StringBuffer object for holding email message body
		StringBuffer message = new StringBuffer();
		message.append("Hello Tomato user,"
				+ "\n\n"
				+ "Our system indicates that you have data files in your Tomato jobs directory that are scheduled for deletion. " 
				+ "\nWe ask that you either move or delete these files before they are removed. Details are below. " 
				+ "\nPlease contact a member of the Bioinformatics Core if this message was received in error, or if you"
				+ "\nneed greater space allocation and also happen to make excellent brownies."
				+ "\n\nThank you," + "\n\nU of U Bioinformatics Core");

		message.append("\n\n\nWarn files (scheduled for deletion):\n");
		
		//check if warnFiles ArrayList is empty
		if (user.getWarnFiles().isEmpty() == true) {
			message.append("\nnone\n");
		}
		
		//for each user with offending files in the warnFiles ArrayList
		for (int i = 0; i < user.getWarnFiles().size(); i++) {
			
			//get files in warnFile ArrayList
			File warnFile = user.getWarnFiles().get(i);
			
			//format lastModifiedDate so it's readable
			long format = warnFile.lastModified();
			DateFormat sdf = SimpleDateFormat.getInstance();
			message.append("\n" + warnFile.getCanonicalPath() + "/" + warnFile.getName() + "\t" + (warnFile.length()/1000000) + "MB\t" + sdf.format(format));
		}
		message.append("\nDeleted files:\n");
		
		//check if deleteFiles ArrayList is empty
		if (user.getDeleteFiles().isEmpty() == true) {
			message.append("\nnone\n");
		}
		
		//for each user with offending files in the deleteFiles ArrayList
		for (int i = 0; i < user.getDeleteFiles().size(); i++) {
			
			//get files in the deleteFile ArrayList
			File deleteFile = user.getDeleteFiles().get(i);
			//format lastModifiedDate so it's readable
			long format = deleteFile.lastModified();
			DateFormat sdf = SimpleDateFormat.getInstance();
			message.append("\n" + deleteFile.getCanonicalPath() + "/" + deleteFile.getName() + "\t" + (deleteFile.length()/1000000) + "MB\t" + sdf.format(format));
		}
		return message.toString();
	}
	
	/**
	 * Emailing method that configures mailing params
	 * @param recipients
	 * @param subject
	 * @param message
	 * @param from
	 * @throws MessagingException
	 */
	public static void postMail(String recipients, String subject, String message, String from) throws MessagingException {
		
		//set the host smtp address
		Properties props = new Properties();
		props.put("mail.smtp.host", "hci-mail.hci.utah.edu");
		
		//create some properties and get the default Session
		javax.mail.Session session = javax.mail.Session.getDefaultInstance(props, null);
		
		//create message
		Message msg = new MimeMessage(session);
		
		//set the from and to address
		InternetAddress addressFrom = new InternetAddress(from);
		msg.setFrom(addressFrom);
		msg.setRecipients(Message.RecipientType.TO, InternetAddress.parse(recipients, false));
		
		//optional: can also set custom headers here in the Email if wanted
		//msg.addHeader("MyHeaderName", "myHeaderValue");
		
		//setting the Subject and Content type
		msg.setSubject(subject);
		msg.setContent(message, "text/plain");
		Transport.send(msg);
	}
	
	/**
	 * Magic. Do not touch.
	 * @param a
	 * @param separator
	 * @return
	 */
	public static String arrayToString(String[] a, String separator) {
		if (a == null || separator == null) {
			return null;
		}
		StringBuilder result = new StringBuilder();
		if (a.length > 0) {
			result.append(a[0]);
			for (int i = 1; i < a.length; i++) {
				result.append(separator);
				result.append(a[i]);
			}
		}
		return result.toString();
	}
	
	/**
	 * Method to send mail to users if they have warnFiles or deleteFiles.
	 * @param user
	 */
	public static void sendMail(User user) {
		String msg = null;
		
		//determine if there are files in the delete files ArrayList
		if (user.getDeleteFiles().isEmpty() != true) {
			String[] emails = new String[user.getEmail().size()];
			
			int i = 0;
			for (String e : user.getEmail()) {
				emails[i] = e;
				i++;
			}
			try {
				//try sending email to user with list of delete files
				Ketchup.postMail(Ketchup.arrayToString(emails, ","), "TOMATO JOB DIRECTORY CLEANUP!!!", Ketchup.getMessage(msg, user), "doNotReply@nobody.com");
			} catch (MessagingException e1) {
				e1.printStackTrace();
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		}
		else {
			//determine if there are files in the warn files ArrayList
			if (user.getWarnFiles().isEmpty() != true) {
				String[] emails = new String[user.getEmail().size()];
						
				int i = 0;
				for (String e : user.getEmail()) {
					emails[i] = e;
					i++;
				}
				//try sending email to user with list of warn files
				try {
					Ketchup.postMail(Ketchup.arrayToString(emails, ","), "TOMATO JOB DIRECTORY CLEANUP!!!", Ketchup.getMessage(msg, user), "doNotReply@nobody.com");
				} catch (MessagingException e1) {
					e1.printStackTrace();
				} catch (IOException e1) {
					e1.printStackTrace();
				}
			}
		}
	}
	
	/**
	 * This method will process each argument and assign new variables.
	 * @param args
	 */
	public void processArgs(String[] args) {
		Pattern pat = Pattern.compile("-[a-z]");
		String programArgs = Misc.stringArrayToString(args, ",");
		boolean verbose = false;
		if (verbose) System.out.println("\nArguments: " + programArgs + "\n");
		for (int i = 0; i < args.length; i++) {
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()) {
				char test = args[i].charAt(1);
				try {
					switch (test) {
					case 'p': path = new String(args[++i]); break;
					case 'h': smtpHostName = new String(args[++i]); break;
					case 's': fileSize = Long.parseLong(args[++i]); break; //file size cutoff
					case 'd': deleteDate = Long.parseLong(args[++i]); break;	//delete date cutoff
					case 'w': warnDate = Long.parseLong(args[++i]); break;	//warn date cutoff
					default: Misc.printErrAndExit("\nProblem--unknown option used!" + mat.group());
					}
				}
				catch (Exception e) {
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -" + test + "\n");
				}
			}
		}
	}
	
	/**
	 * @param args
	 * @throws IOException
	 * @throws MessagingException
	 */
	public static void main(String[] args) throws IOException, MessagingException {
		
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		new Ketchup(args);
	}
	
	public static void printDocs() {
		System.out.println("\n" +
				"**********************************************************************************\n" +
				"**                                Ketchup: Oct 2012                             **\n" +
				"**********************************************************************************\n" +
				"Ketchup identifies users in the Tomato jobs directory, creates a User object for \n" +
				"each, and then for each user, finds and stores their associated email addresses to \n" +
				"be used for notification of impending file deletion. Files are grouped into either \n" +
				"warn files (to be deleted in 2 days) or delete files (passing delete criteria and \n" +
				"removed from file system). Current criteria based upon file size (>1MB by default) \n" +
				"and file last modified date. Users are sent an email with a list of warn files or \n" +
				"those that have been deleted.\n" +

				"\nOptions:\n" +
				"-p Full path for directory where Ketchup will start walking.\n" +
				"   Default: /Users/darren/Desktop/testDir/\n" +
				"-h SMTP host name for emailing purposes.\n"+
				"	Default: mail.inscc.utah.edu\n" +
				"-d Delete date in number of days (int) old relative to current system time.\n" +
				"	Default: 7 days ago.\n"+
				"-w Warn date in number of days (int) old relative to current system time.\n" +
				"	Default: 5 days ago.\n"+
				"-s Minimum file size cutoff for identifying big files.\n"+
				"	Default: 1048576 bytes (1 MB)\n" +
				
				"\nExample: java pathToUSeq/Apps/Ketchup -p /tomato/job/\n" +
				"     -h mail.inscc.utah.edu\n\n" +

		"**************************************************************************************\n");
	}
}