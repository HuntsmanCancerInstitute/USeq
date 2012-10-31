package util.bio.parsers;

import java.io.*;
import java.util.regex.*;
import util.gen.*;


public class MaqSnps2Bed {
	//fields
	private File[] macSnpFiles;
	
	
	public MaqSnps2Bed (String[] args){
		//process args
		processArgs(args);
		
		//for each file
		System.out.println("\nProcessing...");
		for (int i=0; i< macSnpFiles.length; i++){
			System.out.println("\t"+macSnpFiles[i].getName());
			convert(macSnpFiles[i]);
		}
		
		System.out.println("\nDone!\n");
	}
	
	
	public void convert(File file){
		try{
			BufferedReader in = IO.fetchBufferedReader(file);
			File bedFile = new File (file.getParentFile(), Misc.removeExtension(file.getName())+".bed");
			PrintWriter out = new PrintWriter ( new FileWriter(bedFile));
			out.println("#Parsed File = "+file.getCanonicalPath());
			out.println("#Name = Reference base -> Consensus base _ # Reads _ 2nd Best Call _ 3rd Best Call");
			out.println("#Score = Consensus quality -10Log10(p-value)");
			File allelerFile = new File (file.getParentFile(), Misc.removeExtension(file.getName())+"_Alleler.txt");
			PrintWriter outAlleler = new PrintWriter ( new FileWriter(allelerFile));
			outAlleler.println("#Parsed File = "+file.getCanonicalPath());
			outAlleler.println("#Chr Start Stop ReplaceWith #Score #Name");
			outAlleler.println("#Score = Consensus quality -10Log10(p-value)");
			outAlleler.println("#Name = Reference base -> Consensus base _ # Reads _ 2nd Best Call _ 3rd Best Call");
			
			String line;
			String[] tokens;
			Pattern tab = Pattern.compile("\t");
			while ((line=in.readLine()) != null){
				tokens = tab.split(line);
				if (tokens.length != 12) {
					System.out.println("Too few columns. Skipping:\n\t"+line);
					continue;
				}
				String chrom = tokens[0];
				//bowtie -> maqs -> 1 based coordinates
				int start = Integer.parseInt(tokens[1]) -1;
				int stop = start + 1;
				//score = consensus quality
				String score = tokens[4];
				//Name ref->cons_#Reads_2ndBestCall_3rdBestCall 
				String name = tokens[2]+"->"+tokens[3]+"_"+tokens[5]+"_"+tokens[9]+"_"+tokens[11];
				out.println(chrom +"\t"+ start +"\t"+ stop +"\t"+ name +"\t"+ score+ "\t.");
				String name2 = score +" "+name;
				outAlleler.println(chrom +"\t"+ start +"\t"+ stop +"\t"+ tokens[3] +"\t"+ name2);
			}
			in.close();
			out.close();
			outAlleler.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	

	public static void main(String[] args){
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new MaqSnps2Bed(args);
	}

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': macSnpFiles = IO.extractFiles(new File(args[++i])); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e) { Misc.printExit("\nSorry, something doesn't look quite right with this parameter: -"+test+"\n");}
			}
		}
	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              MaqSnps2Bed: June 2009                              **\n" +
				"**************************************************************************************\n" +
				"Converts a Maq snp text file (1 based coordinates) into a bed file (interbase \n" +
				"      coordinates).  Also writes out an Alleler formated text file.\n\n"+

				"-f Full path file text to the file or directory containing Maq snp txt files.\n"+


				"\nExample: java -Xmx1000M -jar path2/USeq/Apps/MaqSnps2Bed -f /data/maqSnpFile.txt\n\n" +

				"\n" +
		"**************************************************************************************\n");		
	}

}
