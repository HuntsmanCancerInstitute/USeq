package trans.roc;
import java.io.*;
import trans.misc.*;
import util.gen.*;

/**
 * Returns oligo values that overlap particular regions.
 */
public class ParseSgrsForParticularScores {
	private double minimum;
	private double maximum;
	private int minRun = 20;
	private int maxGap = 7;
	private File[] files;
	private String chromosome;
	
	public ParseSgrsForParticularScores(String[] args){
		
		//look for arguments
		if (args.length==0) {
			Misc.printExit("\nTo parse xxx.(s)gr files for particular scores enter a min score, max score, and a full path file or directory text.\n");
		}
		
		//parse scores and files
		minimum = Double.parseDouble(args[0]);
		maximum = Double.parseDouble(args[1]);
		System.out.println("\tMin "+minimum+"  Max "+maximum);
		files = IO.extractFiles(new File(args[2]));
		
		//run thru files
		for (int x=0; x<files.length; x++){
			//check if it's a .gr or .sgr file
			String fileName = files[x].getName();
			if (fileName.endsWith(".gr") == false && fileName.endsWith(".sgr") == false) {
				System.out.println("\tSkipping "+fileName);
			}
			else {
				System.out.println("\tParsing "+fileName);
				//assign chromosome
				chromosome= Util.parseChromosomeName(fileName);
				if (chromosome != null){
					System.out.println("\tChromosome "+chromosome);
				}
				try{
					String line;
					BufferedReader in = new BufferedReader(new FileReader(files[x]));
					File parsedFile = new File (files[x].getParentFile(), "parsed_"+fileName);
					PrintWriter out = new PrintWriter(new FileWriter(parsedFile));
					File bedFile;
					PrintWriter bed = null;
					if (chromosome != null){
						bedFile= new File (files[x].getParentFile(), fileName+".bed");
						bed = new PrintWriter(new FileWriter(bedFile));
					}
					File runFile = new File (files[x].getParentFile(), "run_"+fileName);
					PrintWriter outRun = new PrintWriter(new FileWriter(runFile));
					StringBuffer sb = null;
					int runCounter = 0;
					int lastPosition = 0;
					int startBase = 0;
					int endBase = 0;
					//print lines that meet min max
					while ((line=in.readLine())!=null){
						if (line.trim().length()==0) continue;
						String[] tokens = line.split("\\s+");
						double score = Double.parseDouble(tokens[tokens.length-1]);
						//System.out.print(score+" ");
						//good score
						if (score >= minimum && score <= maximum){
							out.println(line);
							//new run?
							if (runCounter == 0){
								sb = new StringBuffer(line +"\n");
								runCounter = 1;
								lastPosition = Integer.parseInt(tokens[tokens.length-2]);
								startBase = lastPosition;
							}
							//old run
							else{
								//check gap
								int position = Integer.parseInt(tokens[tokens.length-2]);
								//max gap OK, add
								if ((position - lastPosition) <= maxGap){
									sb.append(line+ "\n");
									runCounter++;
									lastPosition = position;
								}
								//max gap exceeded, print or dump
								else{
									//minrun met, print
									if (runCounter >= minRun){
										outRun.print(sb);
										//print bed
										if (chromosome != null) {
											endBase = lastPosition;
											bed.println(chromosome+"\t"+startBase+"\t"+endBase);
										}
									}
									startBase = position;
									sb = new StringBuffer(line +"\n");
									runCounter = 1;
									lastPosition = position;
								}
							}
						}
					}
					in.close();
					out.close();
					outRun.close();
					if (chromosome != null) bed.close();
				} catch (Exception e){
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void main(String[] args){
		new ParseSgrsForParticularScores(args);
	}
	
}
