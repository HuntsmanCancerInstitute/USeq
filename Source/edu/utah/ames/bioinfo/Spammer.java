package edu.utah.ames.bioinfo;

import java.io.*;
//import util.gen.*;

/**
 * For sending out email messages
 */
public class Spammer {
	private File list;
	
	public Spammer (String[] args){
		//load contacts
		list = new File ("/Users/nix/JavaStuff/BioRoot/WEB-INF/MiscStuff/UCBContacts.txt");
		String[] contacts = IO.loadFileIntoStringArray(list);
		
		//instantiate Email obj
		//Email email = new Email(Util.smtpHostName, Util.smtpUser, Util.smtpPassword, "text/plain");
		
		//for each contact send email
		for (int i=0; i<contacts.length; i++){
			//parse the last text and email
			String[] tokens = contacts[i].split("\\s+");
			String lastName = tokens[0];
			String[] address = {tokens[tokens.length-1]};
			String message = getMessage(lastName);
			//boolean mailed = email.postMail(address, lastName+" Lab BioRoot Collection!", message, "nix@bioroot.org");
			//if (mailed) System.out.println("Emailed-> "+contacts[i]);
			//else  System.out.println("Failed-> "+contacts[i]);
		}
		System.out.println("Done!");
		
	}
	
	public static void main(String[] args) {
		new Spammer(args);

	}
	
	public static String getMessage(String lastName){
		String message = 
			"Hello Dr. "+lastName+"!\n"+
			"\n"+
			"We have a favor to ask that we believe will actually be a benefit to you.\n"+
			"\n"+
			"How do you keep track of your lab's bio-reagents: oligos, plasmids, strains, and antibodies? Ever worry they are getting lost as your lab grows?  Ever feel that folks are wasting time and lab funds reinventing the bio-reagent wheel?  There is a free and simple solution.\n"+
			"\n"+
			"We have started a non-profit organization called Bioroot Bioinformatics ( http://bioroot.org/ ). As postdocs at the Gallo (Bonci Lab) and UC Berkeley (Drubin and Eisen labs), we know first hand the costs associated with poor molecular biology informatics. Over the past two years, we have developed a free, community service, bio-reagent Laboratory Information Management System (LIMS) to collect, store, and disseminate information about commonly used molecular biology reagents.\n"+  
			"\n"+
			"At its core, BioRoot is a free web driven database that helps labs organize and selectively publish bio-reagents.  The idea is to bring \"cloning by phone\" into the post genomic era by providing a centralized database where people can search for a particular reagent before spending time and money creating it. Every lab has its own collection, every lab member their own account. Each reagent is assigned, by its owner, a visibility setting restricting its view to themselves, their lab mates, or the WWW. Thus, your lab collection is as private or public as you wish. Without exaggeration, BioRoot is the best molecular biology LIMS available. See http://bioroot.org/Documentation.jsp for details.\n"+
			"\n"+
			"Use of the BioRoot LIMS will not only save your lab time and money, it will accelerate science on a global scale.\n"+
			"\n"+
			"One way you can help is to get the "+lastName+" lab to use the bio-reagent LIMS and give us your feedback.\n"+
			"\n"+
			"It's free. At present eight labs are enrolled, including Michael Eisen's, David Drubin's, and Georjana Barnes'.  Once we reach a critical mass, we will apply for grants to subsidize its further development and maintenance. Eventually, our goal is to fold this into the NCBI.\n"+
			"\n"+
			"Would your lab be interested in an account? If you like, we can come by to answer questions and give a demonstration.\n"+
			"\n"+
			"\n"+
			"Best regards,\n"+
			"\n"+
			"David Nix and Billy Chen\n"+
			"\n"+
			"\n"+
			"BioRoot Bioinformatics, Inc.\n"+
			"Berkeley, CA\n"+
			"510-215-1568\n"+
			"\n"+
			"http://bioroot.org\n"+
			"chen@bioroot.org\n"+
			"nix@bioroot.org\n";
		return message;
	}

}
 