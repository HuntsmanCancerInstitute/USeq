package edu.utah.seq.qc;

import java.io.File;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;

public class AggregateQCStats {

	//fields
	private File saveDirectory;
	private File[] jsonFiles;

	private String fastqMatch = "(.+)_fastqCount.json.gz";
	private String saeMatch = "(.+)_samAlignmentExtractor.json.gz";
	private String mpaMatch = "(.+)_mergePairedAlignments.json.gz";
	private String s2uMatch = "(.+)_sam2USeq.json.gz";

	private TreeMap<String, SampleQC> samples;
	private Pattern fastqPattern;
	private Pattern saePattern;
	private Pattern mpaPattern;
	private Pattern s2uPattern;


	//constructors
	public AggregateQCStats(String[] args){

		long startTime = System.currentTimeMillis();
		processArgs(args);
		
		loadSamples();
		checkSamples();

		printStatsTxtSheet();
		printReadCoverageTxtSheet();
		printTargetsTxtSheet();
		
		printHtmlTable();
		printReadCoverageHtml();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" Sec\n");
	}
	
	private void printHtmlThresholdsDescriptions(){
		SampleQC sqc = samples.values().iterator().next();
		StringBuilder sb = new StringBuilder();
		sb.append("<html>\n");
		sb.append("<script type=\"text/javascript\" src=\"https://www.gstatic.com/charts/loader.js\"></script>\n");
		sb.append("<body>\n");

		sb.append("<div id=\"table_legend\"></div><body>\n");
		
		
		
		sb.append("<script>\n");
		sb.append("google.charts.load('current', {'packages':['table']});\n");
		sb.append("google.charts.setOnLoadCallback(drawLegend);\n");
		sb.append("\n");
		sb.append("function drawLegend() {\n");
		sb.append("	var data = new google.visualization.DataTable();\n");
		
		//add columns
		sqc.appendHtmlColumns(sb);
		
		sb.append("	data.addRows([\n");
		
		//for each sample
		for (SampleQC s: samples.values()){
			s.appendHtmlDataRow(sb);
		}
		sb.append("	]);\n");
		sb.append("var table = new google.visualization.Table(document.getElementById('table_main'));\n");
		sb.append("table.draw(data, {title:'Summary Stats', showRowNumber: false, width:'100%', cssClassNames:{headerCell: 'googleHeaderCell'}});\n");
		sb.append("}\n");
		sb.append("</script> \n");
		sb.append("</body>\n");
		sb.append("</html> \n");
		
		File html = new File(saveDirectory, "qcReport.html");
		IO.writeString(sb.toString(), html);
	}
	
	private void printHtmlTable(){
		SampleQC sqc = samples.values().iterator().next();
		StringBuilder sb = new StringBuilder();
		sb.append("<html>\n");
		sb.append("<script type=\"text/javascript\" src=\"https://www.gstatic.com/charts/loader.js\"></script>\n");
		sb.append("<body>\n");

		sb.append("<div id=\"table_main\"></div><body>\n");
		
		sb.append("<script>\n");
		sb.append("google.charts.load('current', {'packages':['table']});\n");
		sb.append("google.charts.setOnLoadCallback(drawTable);\n");
		sb.append("\n");
		sb.append("function drawTable() {\n");
		sb.append("	var data = new google.visualization.DataTable();\n");
		
		//add columns
		sqc.appendHtmlColumns(sb);
		
		sb.append("	data.addRows([\n");
		
		//for each sample
		for (SampleQC s: samples.values()){
			s.appendHtmlDataRow(sb);
		}
		sb.append("	]);\n");
		sb.append("var table = new google.visualization.Table(document.getElementById('table_main'));\n");
		sb.append("table.draw(data, {title:'Summary Stats', showRowNumber: false, width:'100%', cssClassNames:{headerCell: 'googleHeaderCell'}});\n");
		sb.append("}\n");
		sb.append("</script> \n");
		sb.append("</body>\n");
		sb.append("</html> \n");
		
		File html = new File(saveDirectory, "qcReport.html");
		IO.writeString(sb.toString(), html);
	}
	
	private void printReadCoverageHtml(){
		StringBuilder sb = new StringBuilder();
		sb.append("<html>\n");
		sb.append("<script type=\"text/javascript\" src=\"https://www.gstatic.com/charts/loader.js\"></script>\n");
		sb.append("<body>\n");

		sb.append("<div id=\"chart_read_coverage\"></div><body>\n");

		sb.append("<script>\n");
		sb.append("google.charts.load('current', {'packages':['corechart', 'line']});\n");
		sb.append("google.charts.setOnLoadCallback(draw);\n");
		sb.append("\n");
		sb.append("function draw() {\n");
		sb.append("	var rcData = new google.visualization.DataTable();\n");
		
		//how many X axis numbers 0x, 1x, 2x, 3x, 4x,...
		int max25 = 0;
		//for each sample
		for (SampleQC s: samples.values()){
			int hit25 = s.whenHit25ReadCoverage();
			if (hit25> max25) max25= hit25;
		}
		if (max25 < 10) max25 = 10;
		
		//add 1st index X column
		sb.append("rcData.addColumn('number','X')\n");
		//add one for each sample
		for (SampleQC s: samples.values()){
			sb.append("rcData.addColumn('number','"+ s.getSampleName() +"')\n");
		}
		
		//add the data rows
		sb.append("	rcData.addRows([\n");
		int last = max25-1;
		for (int i=0; i< max25; i++){
			sb.append("[");
			sb.append(i);
			for (SampleQC s: samples.values()){
				sb.append(",");
				sb.append(s.getFractionTargetBpsWithIndexedCoverage()[i]);
			}
			if (i != last )sb.append("],");
			else sb.append("]");
		}
		sb.append("	]);\n");
		
		//add options
		sb.append("var options={ title:'Coverage Over Target BPs', hAxis:{title:'Unique Observation Fold Coverage'}, vAxis:{title:'Fraction Target BPs'}, height:600 };\n");
		
		sb.append("var chart = new google.visualization.LineChart(document.getElementById('chart_read_coverage'));\n");
		sb.append("chart.draw(rcData, options);\n");
		sb.append("}\n");
		sb.append("$(window).resize(function(){");
		sb.append("draw()\n");
		sb.append("}\n");
		sb.append("</script> \n");
		sb.append("</body>\n");
		sb.append("</html> \n");
		
		File html = new File(saveDirectory, "qcReportReadCoverage.html");
		IO.writeString(sb.toString(), html);
	}


	private void printStatsTxtSheet() {
		StringBuilder sb = new StringBuilder();
		
		//add header
		SampleQC sqc = samples.values().iterator().next();
		sb.append(sqc.fetchTabbedHeader()); sb.append("\n");
		
		//for each sample
		for (SampleQC s: samples.values()){
			sb.append(s.fetchTabbedLine());
			sb.append("\n");
		}
		sb.append("\n");
		
		//add thresholds
		sb.append("Thresholds:\n");
		sb.append(sqc.fetchThresholds("\t", "\n"));
		sb.append("\n\n");
		
		//add descriptions 
		sb.append("Descriptions:\n");
		sb.append(sqc.fetchDescriptions("", "\t", "\n"));

		File txt = new File(saveDirectory, "qcReport_Stats.xls");
		IO.writeString(sb.toString(), txt);
	}

	private void printReadCoverageTxtSheet() {
		StringBuilder sb = new StringBuilder();

		//read coverage
		sb.append("Coverage\t");
		sb.append(fetchTabbedSampleNames());
		sb.append("\n");

		//how many X axis numbers 0x, 1x, 2x, 3x, 4x,...
		int max25 = 0;
		//for each sample
		for (SampleQC s: samples.values()){
			int hit25 = s.whenHit25ReadCoverage();
			if (hit25> max25) max25= hit25;
		}
		if (max25 < 11) max25 = 11;

		//start at 1 to skip 0x coverage and graph better in excel
		for (int i=1; i< max25; i++){
			sb.append(i);
			for (SampleQC s: samples.values()){
				sb.append("\t");
				sb.append(s.getFractionTargetBpsWithIndexedCoverage()[i]);
			}
			sb.append("\n");
		}

		File txt = new File(saveDirectory, "qcReport_ReadCoverage.xls");
		IO.writeString(sb.toString(), txt);
	}
	
	private void printTargetsTxtSheet() {
		SampleQC sqc = samples.values().iterator().next();
		StringBuilder sb = new StringBuilder();

		//read coverage
		sb.append("Target\t");
		sb.append(fetchTabbedSampleNames());
		sb.append("\n");

		//check the samples and make sure the regions are identical, except for the coverage count
		String[] regions = splitOutCoor(sqc.getExonicMedianPerBpCoverage())[0];
		String[][] counts = new String[samples.size()][regions.length];
		int counter = 0;
		for (SampleQC s: samples.values()){
			String[][] testRegions = splitOutCoor(s.getExonicMedianPerBpCoverage());
			if (regions.length != testRegions[0].length) Misc.printErrAndExit("\nERROR: the number of regions differ between samples?! "+testRegions.length);
			for (int i=0; i< regions.length; i++) if (regions[i].equals(testRegions[0][i])== false) Misc.printErrAndExit("\nERROR: the regions differ between samples?! "+testRegions[0][i]);
			//add em
			counts[counter++] = testRegions[1];
		}
		
		//write em out
		for (int i=0; i< regions.length; i++){
			sb.append(regions[i]);
			//for each sample
			for (int j=0; j< counts.length; j++){
				sb.append("\t");
				sb.append(counts[j][i]);
			}
			sb.append("\n");
		}

		File txt = new File(saveDirectory, "qcReport_PerRegionCoverage.xls");
		IO.writeString(sb.toString(), txt);
	}


	private String[][] splitOutCoor(String[] regions){
		String[] coor = new String[regions.length];
		String[] values = new String[regions.length];
		for (int i=0; i< regions.length; i++){
			String[] tokens = Misc.WHITESPACE.split(regions[i]);
			coor[i] = tokens[0];
			values[i] = tokens[1];
		}
		return new String[][]{coor, values};
	}


	private String fetchTabbedSampleNames() {
		ArrayList<String> al = new ArrayList<String>();
		for (SampleQC s: samples.values()){
			al.add(s.getSampleName());
		}
		return Misc.stringArrayListToString(al, "\t");
	}

	private void loadSamples() {
		//make hash of name and sampleQC
		samples = new TreeMap<String, SampleQC>();

		//make patterns
		fastqPattern = Pattern.compile(fastqMatch);
		saePattern = Pattern.compile(saeMatch);
		mpaPattern = Pattern.compile(mpaMatch);
		s2uPattern = Pattern.compile(s2uMatch);

		//for each json file
		for (File j: jsonFiles){
			//fetch name and type
			String[] nameType = parseSampleName(j.getName());
			//fetch SampleQC
			SampleQC sqc = samples.get(nameType[0]);
			if (sqc == null){
				sqc = new SampleQC(nameType[0]);
				samples.put(nameType[0], sqc);
			}
			sqc.loadJson(j, nameType[1]);
		}
	}

	private void checkSamples() {
		SampleQC test = null;
		//for every sample make sure all four json files were parsed
		for (SampleQC s : samples.values()){
			if (test == null) test = s;
			else {
				//check thresholds
				if (test.checkThresholds(s)){
					Misc.printErrAndExit("\nERROR: this sample's AS, MQ, or AS proc settings differ? Are you changing thresholds between runs? "+s.getSampleName());
				}
				//check json files
				if (test.checkJsonFiles(s)) {
					Misc.printErrAndExit("\nERROR: this sample is missing one or more of the required four xxx.json.gz files from the FastaCounter, SAE, MPA, S2U apps. "+s.getSampleName());
				}
			}
		}
		
	}


	private String[] parseSampleName(String name) {
		
		Matcher mat = fastqPattern.matcher(name);
		if (mat.matches()) return new String[]{mat.group(1), "fastq"};
		
		mat = saePattern.matcher(name);
		if (mat.matches()) return new String[]{mat.group(1), "sae"};
		
		mat = mpaPattern.matcher(name);
		if (mat.matches()) return new String[]{mat.group(1), "mpa"};
		
		mat = s2uPattern.matcher(name);
		if (mat.matches()) return new String[]{mat.group(1), "s2u"};
		
		Misc.printErrAndExit("\nERROR: failed to parse the sample name from "+name +"\nLooking for "+
				saePattern.pattern()+" or "+mpaPattern.pattern()+" or "+s2uPattern.pattern());
		return null;
	}

	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new AggregateQCStats(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		File forExtraction = null;
		System.out.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'j': forExtraction = new File(args[++i]); break;
					case 'r': saveDirectory = new File(args[++i]); break;
					case 'f': fastqMatch = args[++i]; break;
					case 's': saeMatch = args[++i]; break;
					case 'm': mpaMatch = args[++i]; break;
					case 'u': s2uMatch = args[++i]; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//save directory
		if (saveDirectory == null) Misc.printErrAndExit("\nPlease provide a directory to save the results.\n");
		saveDirectory.mkdirs();

		//pull json files
		if (forExtraction == null || forExtraction.canRead() == false) Misc.printExit("\nError: please provide a directory to recurse through to find xxx.json.gz files.\n");
		jsonFiles = IO.fetchFilesRecursively(forExtraction, ".json.gz");
		if (jsonFiles == null || jsonFiles.length ==0 || jsonFiles[0].canRead() == false) Misc.printExit("\nError: cannot find any xxx.json.gz file(s) to parse!\n");

	}	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Aggregate QC Stats: June 2016                         **\n" +
				"**************************************************************************************\n" +
				"Parses and aggregates alignment quality statistics from json files produced by the\n"+
				"SamAlignmentExtractor, MergePairedAlignments, Sam2USeq, and fastq counter.\n"+

				"\nOptions:\n"+
				"-j Directory containing xxx.json.gz files for parsing. Recurses through all other\n"+
				"      directories contained within.\n" +
				"-r Results directory for writing the summary xls spreadsheets for graphing.\n"+

				"\nDefault Options:\n"+
				"-f FastqCount regex for parsing sample name, note the name must be identical across\n"+
				     "the json files, defaults to (.+)_fastqCount.json.gz\n"+
				"-s SAE regex, defaults to (.+)_samAlignmentExtractor.json.gz\n"+
				"-m MPA regex, defaults to (.+)_mergePairedAlignments.json.gz\n"+
				"-u S2U regex, defaults to (.+)_sam2USeq.json.gz\n"+
				"\n"+

				"Example: java -Xmx1G -jar pathToUSeq/Apps/AggregateQCStats -j . -s QCStats/ \n\n" +

				"**************************************************************************************\n");

	}	
}
