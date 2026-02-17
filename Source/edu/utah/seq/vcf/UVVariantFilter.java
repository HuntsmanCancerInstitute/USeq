package edu.utah.seq.vcf;

import java.io.*;
import java.util.regex.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import util.gen.*;
import java.util.*;

public class UVVariantFilter {

	//fields
	private File[] vcfFiles = null;
	private File indexedFasta = null;
	private File saveDirectory = null;
	private HashSet<String> ctContexts = null;
	private HashSet<String> gaContexts = null;

	private IndexedFastaSequenceFile fasta = null;
	//constructors
	public UVVariantFilter(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			makeUvContexts();

			splitVcfs();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem running UvVariantFilter app!");
		}
	}


	public void splitVcfs() throws Exception {
		IO.pl("FileName\t#UV\t#NonUV");
		for (File vcfFile: vcfFiles) {
			//make IO
			File uvVariants = new File(saveDirectory, Misc.removeExtension(vcfFile.getName())+".uv.vcf.gz");
			File cleanVariants = new File(saveDirectory, Misc.removeExtension(vcfFile.getName())+".noUv.vcf.gz");
			Gzipper uv = new Gzipper(uvVariants);
			Gzipper noUv = new Gzipper(cleanVariants);
			BufferedReader in = IO.fetchBufferedReader(vcfFile);
			String[] f;
			String line;
			int numUv = 0;
			int numNonUv = 0;

			while ((line = in.readLine())!=null) {
				line = line.trim();
				if (line.length()==0) continue;

				//header?
				if (line.startsWith("#")) {
					uv.println(line);
					noUv.println(line);
				}
				else {
					f = Misc.TAB.split(line);
					//CHROM	POS	ID	REF	ALT
					//  0    1   2   3   4
					if (f[3].equals("C") && f[4].equals("T")) {
						String tri = fetchTri(f);
						//IO.pl("\n"+line+"\n"+tri);
						if (ctContexts.contains(tri)) {
							uv.println(line);
							numUv++;
						}
						else {
							noUv.println(line);
							numNonUv++;
						}
					}
					else if (f[3].equals("G") && f[4].equals("A")) {
						String tri = fetchTri(f);
						//IO.pl("\n"+line+"\n"+tri);
						if (gaContexts.contains(tri)) {
							uv.println(line);
							numUv++;
						}
						else {
							noUv.println(line);
							numNonUv++;
						}
					}
					else {
						noUv.println(line);
						numNonUv++;
					}
				}

			}
			IO.pl(Misc.removeExtension(vcfFile.getName())+"\t"+numUv+"\t"+numNonUv);
			in.close();
			uv.close();
			noUv.close();
		}
	}


	private String fetchTri(String[] f) {
		long pos = Long.parseLong(f[1]);
		ReferenceSequence rs = fasta.getSubsequenceAt(f[0], pos-1l, pos+1l);
		String tri = new String(rs.getBases());
		return tri;
	}


	private void makeUvContexts() {
		ctContexts = new HashSet<String>();
		gaContexts = new HashSet<String>();

		// TCX: XCC: CCX: XCT: 
		String[] ct = new String[]{"TCA","TCG","TCC","TCT","GCC","ACC","TCC","CCC","CCG","CCA","CCT","CCC","GCT","ACT","TCT","CCT"};
		for (String s: ct) ctContexts.add(s);

		// AGX: XGG: GGX: XGA:   
		String[] ga = new String[]{"AGG","AGA","AGT","AGC","GGG","AGG","TGG","CGG","GGG","GGA","GGT","GGC","GGA","AGA","TGA","CGA"};
		for (String s: ga) gaContexts.add(s);


	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new UVVariantFilter(args);
	}		

	/**This method will process each argument and assign new varibles
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException{
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		File forExtraction = null;
		String source = useqVersion+" Args: USeq/UVVariantFilter "+ Misc.stringArrayToString(args, " ");
		System.out.println("\n"+ source +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 'f': indexedFasta = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;

					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//pull vcf files
		File[][] vcfs = new File[3][];
		vcfs[0] = IO.extractFiles(forExtraction, ".vcf");
		vcfs[1] = IO.extractFiles(forExtraction,".vcf.gz");
		vcfs[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(vcfs);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");

		//check save dir
		if (saveDirectory == null) Misc.printErrAndExit("\nError: cannot find your save directory "+saveDirectory);
		saveDirectory.mkdirs();
		if (saveDirectory.isDirectory() == false) Misc.printErrAndExit("\nError: your save directory does not appear to be a directory");

		//Create fasta fetcher
		if (indexedFasta != null) {
			fasta = new IndexedFastaSequenceFile(indexedFasta);
			if (fasta.isIndexed() == false) Misc.printErrAndExit("\nError: cannot find your xxx.fai fasta index "+ indexedFasta);
		}


	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               UVVariantFilter: Nov 2025                          **\n" +
				"**************************************************************************************\n" +
				"Splits a vcf file into variants with a UV mutation signature and those without. \n"+
				"C->T variants in the context of TCX,XCC,CCX,XCT and G->A in AGX,XGG,GGX,XGA contexts.\n"+

				"\nOptions:\n"+
				"-v Path to your vcf to filter, xxx.vcf(.gz/.zip OK)\n"+
				"-s Path to a directory for saving the parsed results.\n"+
				"-f Path to the reference fasta with xxx.fai index.\n"+


				"\nExample: java -Xmx2G -jar pathToUSeq/Apps/UVVariantFilter -v melanoma.vcf\n" +
				"     -f /Ref/human_g1k_v37.fasta -s ParsedVars\n\n" +

				"**************************************************************************************\n");

	}


}
