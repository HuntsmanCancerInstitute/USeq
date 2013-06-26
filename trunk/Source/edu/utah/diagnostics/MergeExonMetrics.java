package edu.utah.diagnostics;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;


public class MergeExonMetrics {
	private File metricsDir = null;
	private ArrayList<File> metricsList = null;
	private String prefix = null;
	

	public MergeExonMetrics(String[] args) {
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		
		processArgs(args);
		
		mergeFiles();
		
		for (File m: metricsList) {
			m.delete();
		}
		
		System.out.println("Done!");
		
		
	}
	
	private void mergeFiles() {
		try {
					
			ArrayList<String> sectionNames = new ArrayList<String>();
			ArrayList<ArrayList<ArrayList<String>>> sectionData = new ArrayList<ArrayList<ArrayList<String>>>();
			
			Pattern sectionStart = Pattern.compile("<section id='(.+?)'>");
			Pattern sectionEnd = Pattern.compile("</section>");
			Pattern headingPatt = Pattern.compile("<h2>.+?</h2>");
			
			for (int i=0; i<metricsList.size();i++) {
				BufferedReader br = new BufferedReader(new FileReader(metricsList.get(i)));
				String line = null;
				boolean slurp = false;
			
				ArrayList<ArrayList<String>> sample = new ArrayList<ArrayList<String>>();
				ArrayList<String> section = null;
				
				while ((line = br.readLine()) != null) {
					Matcher matcherStart = sectionStart.matcher(line);
					Matcher matcherEnd = sectionEnd.matcher(line);
					Matcher headingMatch = headingPatt.matcher(line);
					if (matcherStart.matches()) {
						slurp = true;
						if (i==0) {
							
							
							sectionNames.add(matcherStart.group(1));
						}
						section = new ArrayList<String>();
					} else if (matcherEnd.matches()) {
						slurp = false;
						sample.add(section);
						
					} else if  (headingMatch.matches()) {
						continue;
					} else if (slurp) {
						section.add(line);
					} 
				}
				
				sectionData.add(sample);
				br.close();
			}
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(this.metricsDir,prefix + ".metrics.html")));
			
			bw.write("<!DOCTYPE html>\n");
			bw.write("<html>\n");
			bw.write("<head>\n");
			bw.write("<title>Sample Metrics: " + prefix + "</title>\n");
			bw.write("<meta name='author' content='USeq ParseExonMetrics'>\n");
			bw.write("</head>\n");
			bw.write("<body>\n");
			bw.write("<h1>Sample Metrics: " + this.prefix + "</h1>\n");
			
			for (int i=0; i<sectionNames.size(); i++) {
				bw.write("<section id='" + sectionNames.get(i) + "'>\n");
				bw.write("<h2>" + sectionNames.get(i) + "</h2>\n");
				for (ArrayList<ArrayList<String>> sample: sectionData) {
					for (String line: sample.get(i)) {
						bw.write(line);
					}
					bw.write("<br />\n");
				}
				bw.write("</section>\n");
			}
			
			
			bw.write("</body>\n");
			bw.write("</html>\n");
			bw.close();
		} catch (IOException ioex) {
			System.out.println("Problem reading/writing metrics files");
			ioex.printStackTrace();
			System.exit(1);
		}
		
	}
	
	public static void main(String[] args) {
		new MergeExonMetrics(args);
	}
	
	private void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          MergeExonMetrics : June 2013                              **\n" +
				"**************************************************************************************\n" +
				"This app simply merges the output from several metrics html files.\n\n\n" +

				"Required:\n"+
				"-f Directory containing metrics html files and a image directory\n"+
				"-o Name of the combined metrics file\n\n" +

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/MergeExonMetrics -f metrics -o 9908_metrics \n" +
		        "**************************************************************************************\n");

	}
	
	private void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
	
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': this.metricsDir = new File(args[++i]); break;
					case 'o': this.prefix = args[++i]; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		if (this.prefix == null) {
			System.out.println("Must specify an output name");
			System.exit(1);
		}
		
		if (this.metricsDir == null) {
			System.out.println("Must specify and metrics directory");
			System.exit(1);
		}
		
		if (!metricsDir.exists()) {
			System.out.println("The specified metrics directory does not exist");
			System.exit(1);
		}
		
		File imageDir = new File(metricsDir,"images");
		
		if (!imageDir.exists()) {
			System.out.println("An image directory does not exist in the metrics directory folder");
			System.exit(1);
		}
		
		File[] mFiles = metricsDir.listFiles();
		ArrayList<String> names = new ArrayList<String>();
		Pattern p = Pattern.compile("(.+?).html");
		metricsList = new ArrayList<File>();
		
		
		//Grab all metrics files in the directory and store their names
		for (File file: mFiles) {
			Matcher m = p.matcher(file.getName());
			if (m.matches()) {
				this.metricsList.add(file);
				names.add(m.group(1));
			}
		}
		
		//Make sure html files are found
		if (names.size() == 0) {
			System.out.println("No files matching *html were found");
			System.exit(1);
		}
		
		
		//Create patterns for each prefix
		ArrayList<Pattern> namePatterns = new ArrayList<Pattern>();
		ArrayList<Integer> nameCount = new ArrayList<Integer>();
		
		for (int i=0;i<names.size();i++) {
			namePatterns.add(Pattern.compile(names.get(i) + ".+?.png"));
			nameCount.add(0);
		}
		
		
		//Check the images directory to see if there are matching image files
		//Check that the number of images files is equal for all samples. 
		//Don't want to specify the exact image name in case metrics change in the future.
		File[] iFiles = imageDir.listFiles();
		
		for (File file: iFiles) {
			for (int i=0;i<namePatterns.size();i++) {
				Pattern np = namePatterns.get(i);
				Matcher m = np.matcher(file.getName());
				if (m.matches()) {
					nameCount.set(i, nameCount.get(i) + 1);
				}
			}
		}
		
		Integer size = null;
		
		for (int i=0;i<nameCount.size();i++) {
			if (nameCount.get(i) == 0) {
				System.out.println("No image files matching prefix: " + names.get(i));
				System.exit(1);
			}
			
			if (size == null) {
				size = nameCount.get(i);
			} else if (size != nameCount.get(i)) {
				System.out.println("Number of image files for each prefix doesn't match");
				System.exit(1);
			}
			
		}
		
		
		
		
			
			
			
		
	}
		

}
