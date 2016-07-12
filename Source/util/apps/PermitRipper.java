package util.apps;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Gmail;
import util.gen.IO;
import util.gen.Misc;

public class PermitRipper {
	
	//fields
	private File tempDir = null;
	private String curlWget = "curl "; //or replace with "wget -O "
	//date parsing and formatting
	private Pattern dateAndAdd = Pattern.compile(".+setNewArrivalDate\\(\"(.+)\", (\\d+)\\).+");
	private Pattern availibilityCode = Pattern.compile(".+class='status (\\w).+");
	private SimpleDateFormat sdf = new SimpleDateFormat("EEE MMM dd");
	private Calendar calendar = Calendar.getInstance();
	//email info
	private HashSet<String> emailedDates = new HashSet<String>();
	private String[] recipients = null;
	private String senderEmail = "david.austin.nix@gmail.com";
	private String numberRafts = null;
	//misc
	private LinkedHashMap<File, String> fileUrl = new LinkedHashMap<File, String>();
	private Random random = new Random();

	
	//URLs to call, must add a date to each
	//Name, URL
	private String[] calls = {
			//"Deso", "http://www.recreation.gov/entranceDetails.do?permitTypeId=1267601734&entrancePermitTypeId=1267621042&entranceId=331126&useTypeId=0&offset=1&entryType=1&contractCode=NRSO&parkId=72440&arvdate=",
			//"Yampa", "http://www.recreation.gov/entranceDetails.do?permitTypeId=2102428458&entrancePermitTypeId=1781929176&entranceId=356818&useTypeId=0&offset=2&entryType=1&contractCode=NRSO&parkId=115139&arvdate=",
			"GatesOfLodore", "http://www.recreation.gov/entranceDetails.do?permitTypeId=2102427966&entrancePermitTypeId=1781929176&entranceId=356817&useTypeId=0&offset=2&entryType=1&contractCode=NRSO&parkId=115139&arvdate=",
			"MiddleFrkSalmon", "http://www.recreation.gov/entranceDetails.do?permitTypeId=523879550&entrancePermitTypeId=524203764&entranceId=292685&useTypeId=0&offset=1&entryType=1&contractCode=NRSO&parkId=75534&arvdate=",
			//"MainSalmon", "http://www.recreation.gov/entranceDetails.do?permitTypeId=523898830&entrancePermitTypeId=524869336&entranceId=292735&useTypeId=0&offset=0&entryType=1&contractCode=NRSO&parkId=75533&arvdate=",
			//"Selway", "http://www.recreation.gov/entranceDetails.do?permitTypeId=523888682&entrancePermitTypeId=524874968&entranceId=292736&useTypeId=0&offset=6&entryType=1&contractCode=NRSO&parkId=75535&arvdate=",
			//"SanJuan", "http://www.recreation.gov/entranceDetails.do?permitTypeId=1782381178&entrancePermitTypeId=1782408288&entranceId=382420&useTypeId=0&offset=0&entryType=1&contractCode=NRSO&parkId=75510&pGroupSize=1&arvdate=",
			//"Snake", "http://www.recreation.gov/entranceDetails.do?permitTypeId=523907650&entrancePermitTypeId=524873263&entranceId=292785&useTypeId=0&offset=0&entryType=1&contractCode=NRSO&parkId=75536&arvdate=",
	};
	//other rivers
	//Salt
	//Dino http://www.freeyampa.org/dates/dinoupdates.js  rip it!
	
	//dates to search, 13 days apart   Aug 13th through 16th
	private String[] dates = {"9/12/2016", "9/25/2016"};
	private String[] datesUnderscore = {"9_12_2016", "9_25_2016"};
	private double maxDaysOld = 3;
	
	public PermitRipper(String[] args){
		long startTime = System.currentTimeMillis();
		System.out.println(Misc.getDateTime());
		
		//set params
		numberRafts = args[0];
		tempDir = new File (args[1]);
		tempDir.mkdirs();
		recipients = Misc.COMMA.split(args[2]);
		
		//load hash of prior emailedDates
		loadSaveEmailedDates(true);
		
		//build script
		System.out.println("Building shell script...");
		String cmd = buildShellScript();

		//execute
		System.out.println("Executing...");
		IO.executeShellScript(cmd, tempDir);
		
		//parse
		System.out.println("Parsing...\n");
		String table = parseEm();
		//email
		if (table.length() !=0){
			table = "<b>Launch Permits Availible as of "+Misc.getDateTime()+"</b><br>\n"+table;
			boolean sent = Gmail.generateAndSendEmail(recipients, "River Permit Notification", table, senderEmail, numberRafts);
			System.out.println("\nEmailing..."+sent);
			
			//save emailedDates
			loadSaveEmailedDates(false);
		}
		else System.out.println("Nada new to email....");
		
		
		
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" Sec\n");
	}
	
	private void loadSaveEmailedDates(boolean load) {
		File hashFile = new File (tempDir, "emailedDatesHash");
		if (load == false) {
			System.out.println("Saving emailed dates hash...");
			IO.saveObject(hashFile, emailedDates);
		}
		else if (hashFile.exists()) {
			//too old?
			double daysOld = IO.numberDaysOld(hashFile);
			System.out.println("Loading emailed dates hash, "+daysOld+" ...");
			if (daysOld < maxDaysOld) {
				emailedDates = (HashSet<String>)IO.fetchObject(hashFile);
			}
			else hashFile.delete();
		}
	}
	
	private String buildShellScript() {
		String[] cmds = new String[(calls.length/2) * dates.length];
		//for each Name, URL
		int index = 0;
		
		for (int i=0; i< calls.length; i++){
			String name = calls[i++];
			String url = calls[i];
			//for each date
			for (int j=0; j< dates.length; j++){
				StringBuilder sb = new StringBuilder();
				//build curl cmd
				File results = new File (tempDir, name+"-"+datesUnderscore[j]);
				results.deleteOnExit();
				String fullUrl = url+dates[j];
				fileUrl.put(results, fullUrl);
				sb.append(curlWget);
				sb.append("'");
				sb.append(fullUrl);
				sb.append("' | grep \"class='status\" > "); 
				sb.append(results.toString());
				sb.append("\n");
				//sleep 0-2sec, lame attempt to avoid blacklisting
				//int num = random.nextInt(3);
				//sb.append("sleep ");
				//sb.append(num);
				//sb.append("s\n");
				cmds[index++] = sb.toString();
			}
		}	
		Misc.randomize(cmds, System.currentTimeMillis());
		return Misc.stringArrayToString(cmds, "");
	}

	private String parseEm() {
		StringBuilder sb = new StringBuilder();
		//for each file
		for (File f: fileUrl.keySet()){
			String url = fileUrl.get(f);
			String name = f.getName();
			if (f.exists() == false) Misc.printErrAndExit("Error: couldn't find curl/ wget output file for "+f.getName());
			//parse dates
			String[] dates = parseDates(f);
			//print to screen
			System.out.println(name+"\t"+dates[0]);
			//any launches?
			if (dates.length == 1) continue;
			//look to see if already emailed
			ArrayList<String> datesToSend = new ArrayList<String>();
			datesToSend.add(dates[0]);
			for (int i=1; i< dates.length; i++){
				String date = dates[i];
				String id = name+date;
				if (emailedDates.contains(id) == false){
					emailedDates.add(id);
					datesToSend.add(date);
				}
			}
			//any dates?
			if (datesToSend.size() == 1) continue;
			String linkName = name.replace("-", " ");
			//build link 
			sb.append("<a href='");
			sb.append(url);
			sb.append("'>");
			sb.append(linkName);
			sb.append("</a>&nbsp;:&nbsp;");
			sb.append(Misc.stringArrayListToString(datesToSend, ", "));
			sb.append("<br>\n");
		}
		return sb.toString();
	}

	private String[] parseDates(File f) {
		/*
		A: Available for online reservation (click to book entry date)
	L: Accepting Lottery Application (click to apply for the lottery)
	W: Available at the Facility
	R: Reserved
	C: Closed
	X: Not available
		 */
		//code counters for Codes [r, a, l, x, w]
		int r = 0;
		int a = 0; //or w
		int l = 0;
		int x = 0;
		int w = 0;
		int c = 0;
		
		//load file
		String[] lines = IO.loadFile(f);
		ArrayList<String> al = new ArrayList<String>();
		
		//parse the availability code and date from each
		for (int i=0; i< lines.length; i++){
			//available?
			Matcher avail = availibilityCode.matcher(lines[i]);
			if (avail.matches() == false) Misc.printErrAndExit("Error: failed to parse an availability from line\n"+lines[i]+"\nin file\n"+f.getName());

			//parse date?
			if (avail.group(1).equals("a")){
				Matcher pat = dateAndAdd.matcher(lines[i]);
				if (pat.matches() == false) Misc.printErrAndExit("Error: failed to parse a date from line\n"+lines[i]+"\nin file\n"+f.getName());
				Date date = new Date (pat.group(1));
				int daysToAdd = Integer.parseInt(pat.group(2));
				al.add(addDays(date, daysToAdd));
				a++;
			}
			else if (avail.group(1).equals("l")) l++;
			else if (avail.group(1).equals("x")) x++;
			else if (avail.group(1).equals("r")) r++;
			else if (avail.group(1).startsWith("w")) w++;
			else if (avail.group(1).startsWith("c")) c++;
			else Misc.printErrAndExit("Error: unknown availibility code in line\n"+lines[i]+"\nin file\n"+f.getName());
		}
		//build return
		String codeCount = "A"+a+":W"+w+":L"+l+":R"+r+":X"+x+":C"+c;
		al.add(0, codeCount);
		String[] availDates = new String[al.size()];
		al.toArray(availDates);
		return availDates;
	}

	public static void main(String[] args) {
		if (args.length == 0) System.out.println("Password TempDir CommaListRecipientsEmail");
		else new PermitRipper(args);
	}
	
	public String addDays(Date date, int days){
        calendar.setTime(date);
        calendar.add(Calendar.DATE, days); 
        return sdf.format(calendar.getTime());
    }

}
