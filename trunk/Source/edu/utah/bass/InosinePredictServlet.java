package edu.utah.bass;
import java.io.*;
import java.util.*;
import javax.servlet.*;
import javax.servlet.http.*;
import org.apache.commons.fileupload.*;
import util.bio.parsers.MultiFastaParser;
import util.gen.IO;


/** Web app interface for the InosinePredict.java command line app.
 @author davidnix
 */
public class InosinePredictServlet extends HttpServlet {

	//fields
	public File matrixFileDirectory;
	public File[] matrixFiles;
	public String[] matrixFileNames;
	public File resultsDirectory;
	public long maxFileAge = 1000 * 60 * 60 ;

	//methods


	public void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		doGet(request, response);
	}

	public void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		response.setContentType("text/html");
		PrintWriter out = response.getWriter();
		
		deleteOldFiles();

		//if empty request then just print form
		if (FileUpload.isMultipartContent(request) == false ) out.println(getForm(null, null));

		//process submission
		else {
			StringBuilder error = new StringBuilder();
			//make a results directory
			File tempResultsDirectory = new File (resultsDirectory, "IP_"+System.currentTimeMillis());
			tempResultsDirectory.mkdir();
			File sequenceFile = new File (tempResultsDirectory, "sequence.fasta");
			
			//will null sequenceFile if none was uploaded
			HashMap<String,String> params = processMultiPartFile(request, sequenceFile);

			//fetch matrix file
			File matrixFile = null;
			String matrixIndex = params.get("matrix");
			if (matrixIndex != null) matrixFile = matrixFiles[Integer.parseInt(matrixIndex)];
			else {
				error.append("<FONT COLOR='#B22222' SIZE=3>Failed to identify a matrix preference?!</FONT>\n");
				System.err.println("Failed to identify a matrix preference?!");
			}

			
			//if sequence file not found then look for sequence from text area
			if (sequenceFile.exists() == false){ 
				String sequence = params.get("sequence");
				if (isEmpty(sequence)) {
					error.append("<FONT COLOR='#B22222' SIZE=3>Please provide one or more fasta formatted sequences in the text box or as a file.</FONT>\n");
				}
				else writeString(sequence, sequenceFile);
			}
			
			//look for error and check sequence
			InosinePredict ip = null;
			if (error.length() == 0){
				MultiFastaParser mfp = new MultiFastaParser(sequenceFile);
				if (mfp.isFastaFound() == false){
					error.append("<FONT COLOR='#B22222' SIZE=3>Your (multi) fasta sequence is malformed, please correct and resubmit.</FONT>\n");
				}
				else {
					out.println("Launching...");
					File zipArchive = new File(resultsDirectory, tempResultsDirectory.getName()+".zip");
					ip = new InosinePredict(sequenceFile, mfp, matrixFile, tempResultsDirectory, zipArchive);
				}
			}
			
			if (error.length() != 0) out.println(getForm(error.toString(), null));
			else {
				//output resultsPage
				out.println(getResults(ip));
			}
			
			
			
			
			
		}
	}
	
 	public static boolean isEmpty(String s){
		if (s == null || s.trim().equals("")) return true;
		return false;
	}

	public static HashMap<String,String> processMultiPartFile(HttpServletRequest request, File file){
		try {
			FileItem fi = null;
			DiskFileUpload upload = new DiskFileUpload();
			//Process multipart params
			List items = upload.parseRequest(request);
			Iterator iter = items.iterator();
			HashMap<String,String> nameValue = new HashMap<String,String>();
			while (iter.hasNext()) {
				FileItem item = (FileItem) iter.next();
				if (item.isFormField()) {
					nameValue.put(item.getFieldName(), item.getString());
				}
				else{
					fi = item;
					if (fi.getName()!="") fi.write(file);
					else file = null;
				}
			}
			if (fi == null) file = null;
			return nameValue;
		} catch (Exception e) {
			System.err.println("Error writing file to disk!");
			e.printStackTrace();
			file = null;
			return null;
		}
	}

	/**Writes a String to disk. */
	public static void writeString(String data, File file) {
		try {
			PrintWriter out = new PrintWriter(new FileWriter(file));
			out.print(data);
			out.close();
		} catch (IOException e) {
			System.err.println("Problem writing sequence to disk.");
			e.printStackTrace();
			file = null;
		}
	}
	
	public String getResults(InosinePredict ip){
		String tempFolder = ip.getSaveDirectory().getName();
		String zipArchive = tempFolder+".zip";
		
		//the form
		String top = 
		"<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n"+
		"<html xmlns=\"http://www.w3.org/1999/xhtml\">\n"+
		"<head>\n"+
		"<meta http-equiv=\"Content-Type\" content=\"text/html; charset=UTF-8\" />\n"+
		"<title>Inosine Predict Application</title>\n"+
		"<link href=\"IPApp/iPredict_template.css\" rel=\"stylesheet\" type=\"text/css\" />\n"+
		"<script src=\"IPApp/SpryAssets/SpryValidationTextarea.js\" type=\"text/javascript\"></script>\n"+
		"<script src=\"IPApp/SpryAssets/SpryValidationTextField.js\" type=\"text/javascript\"></script>\n"+
		"<script src=\"IPApp/SpryAssets/SpryValidationRadio.js\" type=\"text/javascript\"></script>\n"+
		"<link href=\"IPApp/SpryAssets/SpryValidationTextarea.css\" rel=\"stylesheet\" type=\"text/css\" />\n"+
		"<link href=\"IPApp/SpryAssets/SpryValidationTextField.css\" rel=\"stylesheet\" type=\"text/css\" />\n"+
		"<link href=\"IPApp/SpryAssets/SpryValidationRadio.css\" rel=\"stylesheet\" type=\"text/css\" />\n"+
		"</head>\n"+
		"\n"+
		"<body>\n"+
		"\n"+
		"<div id=\"container\">\n"+
		"	\n"+
		"<!-- Begin Masthead -->\n"+
		"	\n"+
		"<div id=\"header\">\n"+
		"<p><span class=\"deemphasized\">InosinePredict</span><strong></strong></p>\n"+
		"</div>\n"+
		"\n"+
		"<!-- End Masthead, Begin mainMenu -->\n"+
		"\n"+
		"<div id=\"navigation\">\n"+
		"<ul>\n"+
		"<li><a href=\"http://medicine.utah.edu/\" target=\"_blank\"><img src=\"IPApp/Images/Ulogo.png\" alt=\"\" width=\"26\" height=\"30\" hspace=\"5\" vspace=\"0\" border=\"0\" align=\"middle\" />   University of Utah School of Medicine</a>\n"+
		"</li></ul>\n"+
		"</div>\n"+
		"\n"+
		
		"<!-- Main content -->\n"+
		"\n"+
		"<div id=\"content\">\n"+
		"\n"+
		"<h2>Results</h2>\n"+
		"<p>InosinePredict derives from experimentally determined editing sites within a perfectly base-paired dsRNA reacted to 20% overall adenosine to inosine editing, using 4 different versions of human ADARs.</p>\n"+
		"<p>Actual editing sites may vary from those predicted, especially under different reaction conditions, including:</p>\n"+
		"<ul type='circle'>\n"+
			"<li>Reactions with ADARs from different organisms</li>\n"+
			"<li>Reactions performed to &lt; 20% overall editing or &gt; 20% overall editing</li>\n"+
			"<li>Reactions using dsRNA substrates containing mismatches, bulges and loops</li>\n"+
		"</ul>\n"+

		"<p><strong><font color='#CCCC99'>Important and of Interest:</font></strong> Comparing results with InosinePredict to experimentally determined editing sites may reveal interesting features that affect editing, in addition to neighboring bases. For more details, please see Eggington et al.</p>\n"+

		"<p>Right click the following to download a <a href='Results/"+zipArchive+"'><FONT SIZE='+1'>zip archive of your results</FONT></a>. Hint: Import the eps files into " +
				"Word/ Excel/ PowerPoint/ Illustrator for downstream use. Grayed bases are impossible to calculate due to end proximity or an ambigous base in the submission.</p>\n"+
		"\n";
		
		//for each sequence
		File[] pngs = ip.getPngFiles();
		int[] widths = ip.getWidths();
		StringBuilder sb = new StringBuilder();
		for (int i=0; i< pngs.length; i++){
			String name = pngs[i].getName().replace(".png", "");
			sb.append("<h2>>");
			sb.append(name);
			sb.append("</h2>\n<img ");
			if (widths[i] > 500) sb.append("width=800 ");
			sb.append("src=Results/");
			sb.append(tempFolder);
			sb.append("/");
			sb.append(pngs[i].getName());
			sb.append(">\n<p>&nbsp;</p>\n");
		}
		
		String bottom =
		"<a href=\"\"><FONT COLOR='#B22222' SIZE='+1'>Click to return to Submission Form</FONT></a><p>&nbsp;</p>\n"+
		"<div id=\"footer\">\n"+
		"  <ul>\n"+
		"<li><a href=\"http://medicine.utah.edu/\" target=\"_blank\">University of Utah SOM</a></li>\n"+
		"<li><a href=\"http://medicine.utah.edu/biochemistry/\" target=\"_blank\">Dept of Biochemistry</a></li>\n"+
		"<li><a href=\"http://www.biochem.utah.edu/bass/\" target=\"_blank\">Brenda L. Bass Lab</a></li>\n"+
		"<li><a href=\"http://www.utah.edu/portal/site/uuhome/menuitem.4694b7a3dd66f40516df1210d1e916b9/?vgnextoid=960992d315bb3110VgnVCM1000001c9e619bRCRD\" target=\"_blank\">Disclaimer</a></li>\n"+
		"<li><script type=\"text/javascript\">\n"+
		"var e = unescape(\"%3Db%21isfg%3E%23nbjmup%3BscpplnboAcjpdifn/vubi/fev%23%3FXfcnbtufs%3D0b%3F\");\n"+
		"var i,p='';for(i=0;i<e.length;i++){p+=String.fromCharCode(((e.charCodeAt(i)-33)%240)+32);}\n"+
		"document.write(p);\n"+
		"</script><noscript>david<em>dot</em>nix <em>-at-</em> hci <em>dot</em> utah <em>dot</em> edu</noscript></li> \n"+
		"  </ul>\n"+
		"</div>\n"+
		"<!-- End Footer -->\n"+
		"</div>\n"+
		"</body>\n"+
		"</html>\n";
		
		return top+sb.toString()+bottom;
	}

	
	public String getForm(String errorMessage, String sequence){
		String message = "";
		if (errorMessage != null){
			message = "<h2>Error! </h2>\n"+
			"<p>"+errorMessage+"</p>\n";
		}
		if (sequence == null) sequence = ">SequenceName\nGAAUUAACCAAGGAAAAUAACAAGGACAGGGACCAGG";
		//the form
		return 
		"<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n"+
		"<html xmlns=\"http://www.w3.org/1999/xhtml\">\n"+
		"<head>\n"+
		"<meta http-equiv=\"Content-Type\" content=\"text/html; charset=UTF-8\" />\n"+
		"<title>Inosine Predict Application</title>\n"+
		"<link href=\"IPApp/iPredict_template.css\" rel=\"stylesheet\" type=\"text/css\" />\n"+
		"<script src=\"IPApp/SpryAssets/SpryValidationTextarea.js\" type=\"text/javascript\"></script>\n"+
		"<script src=\"IPApp/SpryAssets/SpryValidationTextField.js\" type=\"text/javascript\"></script>\n"+
		"<script src=\"IPApp/SpryAssets/SpryValidationRadio.js\" type=\"text/javascript\"></script>\n"+
		"<link href=\"IPApp/SpryAssets/SpryValidationTextarea.css\" rel=\"stylesheet\" type=\"text/css\" />\n"+
		"<link href=\"IPApp/SpryAssets/SpryValidationTextField.css\" rel=\"stylesheet\" type=\"text/css\" />\n"+
		"<link href=\"IPApp/SpryAssets/SpryValidationRadio.css\" rel=\"stylesheet\" type=\"text/css\" />\n"+
		"</head>\n"+
		"\n"+
		"<body>\n"+
		"\n"+
		"<div id=\"container\">\n"+
		"	\n"+
		"<!-- Begin Masthead -->\n"+
		"	\n"+
		"<div id=\"header\">\n"+
		"<p><span class=\"deemphasized\">InosinePredict</span><strong></strong></p>\n"+
		"</div>\n"+
		"\n"+
		"<!-- End Masthead, Begin mainMenu -->\n"+
		"\n"+
		"<div id=\"navigation\">\n"+
		"<ul>\n"+
		"<li><a href=\"http://medicine.utah.edu/\" target=\"_blank\"><img src=\"IPApp/Images/Ulogo.png\" alt=\"\" width=\"26\" height=\"30\" hspace=\"5\" vspace=\"0\" border=\"0\" align=\"middle\" />   University of Utah School of Medicine</a>\n"+
		"</li></ul>\n"+
		"</div>\n"+
		"\n"+
		"<!-- Main content -->\n"+
		"\n"+
		"<div id=\"content\">\n"+
		"\n"+
		"<h2>Description</h2>\n"+
		"<p><em>InosinePredict</em> uses an algorithm based on an experimentally determined <em>in vitro</em> editing of a long, perfectly base-paired dsRNA, as described in:</p>\n"+
		"\n"+
		"<p><strong><em>Predicting sites of ADAR editing in double-stranded RNA</em></strong><br />\n"+
		"Julie M. Eggington, Thomas Greene, Brenda L. Bass, <em>2011 Nat Commun 2:319</em>.</p>\n"+
		"\n"+
		"<p>Predictions take into account the identity of the four bases 5&#8217; and four\n"+
		"    bases 3&#8217; of an adenosine, and are based on experimentally determined editing\n"+
		"    site preferences for human ADAR1, human ADAR2, and proteins truncated to\n"+
		"    only contain the catalytic domain (hADAR1-D and hADAR2-D). </p>\n"+
		"  \n"+
		"<h2>Instructions</h2>\n"+
		" \n"+
		" <p>To apply the algorithm to your RNA sequence, upload\n"+
		"   a file in <a href=\"http://en.wikipedia.org/wiki/FASTA_format\" target=\"_blank\">FASTA or Multi-FASTA\n"+
		"     format</a> <strong>OR</strong> cut and paste your FASTA sequence in the box\n"+
		"     provided in step 1.&nbsp; Observe the following:</p>\n"+
		" <ul>\n"+
		"   <li>A, C, U, G bases only</li>\n"+
		"   <li>Enter only one strand of your dsRNA (5&#8217; to 3&#8217;)</li>\n"+
		"</ul>\n"+
		"</p>\n"+
		"\n"+
		message +
		"<h2>Step 1</h2>\n"+
		"<p>Cut and paste your FASTA formatted sequence(s) here:</p>\n"+
		"<form action=\"\" method=\"post\" enctype=\"multipart/form-data\">\n"+
		"  <span id=\"sprytextarea1\">\n"+
		"     <textarea name=\"sequence\" cols=\"80\" rows=\"6\" id=\"sequence\" tabindex=\"1\">"+ sequence +"</textarea>\n"+
		"</span>\n"+
		"  <p>&nbsp;<strong>OR</strong></p>\n"+
		"  <span id=\"sprytextfield1\">\n"+
		"  <p>Upload a file in FASTA format: &nbsp;&nbsp;&nbsp;<input name=\"fileName\" type=\"file\" id=\"fileName\" tabindex=\"2\" size=\"40\" /><p/>\n"+
		"</span>\n"+
		"\n"+
		"  <h2>Step 2</h2>  \n"+
		"  <p>Select an ADAR neighbor preference table to predict editing.\n"+
		"  <p><span id=\"spryradio1\">\n"+
		"    <label>\n"+
		"      &nbsp; \n"+
		"      &nbsp; \n"+
		"      <input name=\"matrix\" type=\"radio\" id=\"0\" value=\"0\" checked=\"checked\" />\n"+
		matrixFileNames[0] + "</label>\n"+
		"    <br />\n"+
		"    <label>\n"+
		"      &nbsp; \n"+
		"      &nbsp; \n"+
		"      <input type=\"radio\" name=\"matrix\" value=\"1\" id=\"1\" />\n"+
		matrixFileNames[1] + "</label>\n"+
		"    <br />\n"+
		"    <label>\n"+
		"      &nbsp; \n"+
		"      &nbsp; \n"+
		"      <input type=\"radio\" name=\"matrix\" value=\"2\" id=\"2\" />\n"+
		matrixFileNames[2] + "</label>\n"+
		"    <br />\n"+
		"    <label>\n"+
		"      &nbsp; \n"+
		"      &nbsp; \n"+
		"      <input type=\"radio\" name=\"matrix\" value=\"3\" id=\"3\" />\n"+
		matrixFileNames[3] + "</label>\n"+
		"    <br />\n"+
		"    <span class=\"radioRequiredMsg\">Please make a selection.</span></span></p>\n"+
		"<p>&nbsp;</p>\n"+
		"\n"+
		"<input type=\"submit\" name=\"submit\" value=\"Launch InosinePredict!\">\n"+
		"</form>\n"+
		"\n"+
		"<p>&nbsp;</p>\n"+
		"  <h2>Considerations</h2>\n"+
		"  <p><em>InosinePredict</em> will predict editing sites in the RNA sequence uploaded,\n"+
		"      as well as the complementary strand, assuming Watson-Crick base-pairing.\n"+
		"      At present, GU base-pairs are not considered, but can be accounted for by\n"+
		"      entering each strand separately.</p>\n"+
		"<p> <em>InosinePredict</em> requires at least four neighboring bases on each side\n"+
		"  of the relevant adenosine. Adenosines less than 4 nucleotides from 5&#8217; or\n"+
		"  3&#8217; termini are not considered. </p>\n"+
		" \n"+
		"</div>\n"+
		"\n"+
		"<div id=\"footer\">\n"+
		"  <ul>\n"+
		"<li><a href=\"http://medicine.utah.edu/\" target=\"_blank\">University of Utah SOM</a></li>\n"+
		"<li><a href=\"http://medicine.utah.edu/biochemistry/\" target=\"_blank\">Dept of Biochemistry</a></li>\n"+
		"<li><a href=\"http://www.biochem.utah.edu/bass/\" target=\"_blank\">Brenda L. Bass Lab</a></li>\n"+
		"<li><a href=\"http://www.utah.edu/portal/site/uuhome/menuitem.4694b7a3dd66f40516df1210d1e916b9/?vgnextoid=960992d315bb3110VgnVCM1000001c9e619bRCRD\" target=\"_blank\">Disclaimer</a></li>\n"+
		"<li><script type=\"text/javascript\">\n"+
		"var e = unescape(\"%3Db%21isfg%3E%23nbjmup%3BscpplnboAcjpdifn/vubi/fev%23%3FXfcnbtufs%3D0b%3F\");\n"+
		"var i,p='';for(i=0;i<e.length;i++){p+=String.fromCharCode(((e.charCodeAt(i)-33)%240)+32);}\n"+
		"document.write(p);\n"+
		"</script><noscript>david<em>dot</em>nix <em>-at-</em> hci <em>dot</em> utah <em>dot</em> edu</noscript></li> \n"+
		"  </ul>\n"+
		"</div>\n"+	
		"<!-- End Footer -->\n"+
		"<!-- End Container -->\n"+
		"</div>\n"+
		
		"</body>\n"+
		"</html>\n";
	}
	
	public void deleteOldFiles(){
		File[] files = IO.fetchFilesRecursively(resultsDirectory);
		long currentTime = System.currentTimeMillis();
		for (int i=0; i< files.length; i++){
			long fileTime = files[i].lastModified();
			if ((currentTime- fileTime) > maxFileAge) {
				if (files[i].isFile() && files[i].exists()) files[i].delete();
				else if (files[i].isDirectory() && files[i].exists()) {
					IO.deleteDirectory(files[i]);
				}
			}
		}
	}

	/*initialize the servlet*/
	public void init(){
		System.out.println("############### Initializing the InosinePredictServlet.java ################");

		//look for directory
		matrixFileDirectory = new File (getServletContext().getRealPath("/")+"/MatrixFiles");
		if (matrixFileDirectory.exists() == false) System.err.println("\nError: matrixFileDirectory doesn't exist? -> "+matrixFileDirectory);

		//look for results, delete it if it exists, for some reason, dirs don't get deleted
		resultsDirectory = new File (getServletContext().getRealPath("/")+"/Results");
		if (resultsDirectory.exists()) {
			IO.deleteDirectory(resultsDirectory);
			resultsDirectory = new File (getServletContext().getRealPath("/")+"/Results");
		}
		resultsDirectory.mkdir();

		//look for matrix files
		matrixFileNames = new String[] {
				"hADAR1", 
				"hADAR2",
				"hADAR1-D",
				"hADAR2-D"
		};
		matrixFiles = new File[4];
		for (int i=0; i< matrixFiles.length; i++) {
			matrixFiles[i] = new File (matrixFileDirectory, matrixFileNames[i]);
			if (matrixFiles[i].exists() == false) System.err.println("\nError: matrixFile doesn't exist? -> "+matrixFiles[i]);
		}
	}

	/*/Applications/apache-tomcat-5.5.28/bin/bootstrap.jar
nix:~ davidnix$ cp /Users/davidnix/Code/BioTools/util/gen/IO.class /Applications/apache-tomcat-5.5.28/webapps/Bass/WEB-INF/classes/util/gen/IO.class ; cp /Users/davidnix/Code/BioTools/util/gen/Misc.class /Applications/apache-tomcat-5.5.28/webapps/Bass/WEB-INF/classes/util/gen/Misc.class ; cp /Users/davidnix/Code/BioTools/util/gen/Num.class /Applications/apache-tomcat-5.5.28/webapps/Bass/WEB-INF/classes/util/gen/Num.class ; cp /Users/davidnix/Code/BioTools/util/gen/Swing.class /Applications/apache-tomcat-5.5.28/webapps/Bass/WEB-INF/classes/util/gen/Swing.class ; cp /Users/davidnix/Code/BioTools/util/bio/parsers/MultiFastaParser.class /Applications/apache-tomcat-5.5.28/webapps/Bass/WEB-INF/classes/util/bio/parsers/MultiFastaParser.class ; cp /Users/davidnix/Code/BioTools/edu/utah/bass/*class /Applications/apache-tomcat-5.5.28/webapps/Bass/WEB-INF/classes/edu/utah/bass ; /Applications/apache-tomcat-5.5.28/bin/shutdown.sh ; rm /Applications/apache-tomcat-5.5.28/logs/catalina.out ; touch /Applications/apache-tomcat-5.5.28/logs/catalina.out ; sleep 3s; /Applications/apache-tomcat-5.5.28/bin/startup.sh*/

}
