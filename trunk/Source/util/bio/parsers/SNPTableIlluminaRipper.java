package util.bio.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

public class SNPTableIlluminaRipper {

	//fields
	private File cgFile;
	private double minimGCScore = 0.9;
	private Pattern tab = Pattern.compile("\t");

	//constructor
	public SNPTableIlluminaRipper (String[] args){
		//rip dbSNPs for uniques
		//parseKeyFilesForUniques(new File(args[0]));
	
		//rip ucsc bed file of snps for those that are unique
		fetchCoordinates(new File(args[0]), new File(args[1]));
		
		//make file object
		//cgFile = new File (args[0]);

		//parse
		//parseFile();
	}

	/*
	 * SNP Name        Sample ID       Allele1 - Forward       Allele2 - Forward       Allele1 - Top   Allele2 - Top   GC Score
rs1000000	28114	C	C	G	G	0.8228
rs1000002	28114	G	G	G	G	0.8822
rs10000023	28114	T	T	A	A	0.8360
rs1000003	28114	A	G	A	G	0.8951
rs10000030	28114	G	G	G	G	0.6725
rs10000037	28114	A	G	A	G	0.9332
rs10000041	28114	T	T	A	A	0.7791
rs10000042	28114	C	C	G	G	0.8377
rs10000049	28114	A	C	A	C	0.9424
rs1000007	28114	A	G	A	G	0.9490
rs10000073	28114	T	T	A	A	0.9508

	 */
	
	private void fetchCoordinates (File key, File dbSNPBed){
		System.out.println("Loading key...");
		HashSet<String> goodRsIDs = IO.loadFileIntoHashSet(key);
		String line = null;
		try {
			System.out.println("Scanning snps...");
			BufferedReader in = IO.fetchBufferedReader(dbSNPBed);
			String name = dbSNPBed.getName();
			name = Misc.removeExtension(name);
			name = name+"_NoDups.txt.gz";
			Gzipper out = new Gzipper(new File (dbSNPBed.getParentFile(), name));
			String[] tokens;
			
			while ((line = in.readLine()) !=null){
				line = line.trim();
				if (line.length() == 0 || line.startsWith("#")) continue;
				tokens = tab.split(line);
				//is it in the good ids set?
				if (goodRsIDs.contains(tokens[3])) out.println(line);
			}
			out.close();
		} catch (IOException e) {
			System.err.println("BadLine -> "+line);
			e.printStackTrace();
		}
		
	}

	private void parseFile() {
		try {
			BufferedReader in = IO.fetchBufferedReader(cgFile);
			HashMap<String, PrintWriter> sampleIDs = new HashMap<String, PrintWriter>();
			String line;
			String[] tokens;
			
			
			//advance to "SNP Name"
			while ((line = in.readLine()) !=null){
				if (line.startsWith("SNP Name")){
					System.out.println(line);
					break;
				}
			}
			
			PrintWriter out = null;
			while ((line = in.readLine()) !=null){
				line = line.trim();
				tokens = tab.split(line);
				String sample = tokens[1];
				
				//sampleIDs.add(sample);
				
			}

	
			//print varTypes
			System.out.println("Found the following samples:");
			System.out.println(sampleIDs);
			
			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}



	}
	
	
	
	private void parseKeyFilesForUniques(File dir){
		File[] files = IO.extractFiles(dir, "gz");
		
		for (File f: files){
			System.out.println("\t"+f.getName());
			String name = f.getName();
			name = Misc.removeExtension(name);
			name = name+"_NoDups.txt.gz";
			try {
				BufferedReader in = IO.fetchBufferedReader(f);
				Gzipper out = new Gzipper(new File (f.getParentFile(), name));
				String oldRS = "";
				String oldChrom = "";
				int oldStart = -1;
				String line;
				String[] tokens;
				boolean oldIsDup = false;
				
				//load first data line
				while ((line = in.readLine()) != null){
					if (line.startsWith("#") || line.length() ==0) continue;
					tokens = tab.split(line);
					//new line chrY	19096363	19096363	rs3894		
					oldRS = tokens[3];
					oldChrom = tokens[0];
					oldStart = Integer.parseInt(tokens[1]) -1; // convert to interbase
					break;
				}
				
				while ((line = in.readLine()) != null){
					tokens = tab.split(line);
					//new rs id?
					if (tokens[3].equals(oldRS) == false){
						//is it not a dup?
						if (oldIsDup == false){
							//print it
							//out.print(oldChrom);
							//out.print("\t");
							//out.print(new Integer (oldStart).toString());
							//out.print("\t");
							out.println(oldRS);
							
						}
						oldRS = tokens[3];
						oldChrom = tokens[0];
						oldStart = Integer.parseInt(tokens[1]) -1; // convert to interbase
						oldIsDup = false;
					}
					//nope it's a dup
					else oldIsDup = true;
				}
				//print last?
				if (oldIsDup == false){
					//out.print(oldChrom);
					//out.print("\t");
					//out.print(new Integer (oldStart).toString());
					//out.print("\t");
					out.println(oldRS);
				}
				
				out.close();
				in.close();
			} catch (Exception e) {
				
				e.printStackTrace();
			} 
		}
		
	}


	// main method called with this class
	public static void main (String[] arg){
		new SNPTableIlluminaRipper (arg);
	}

}
