package meme;
import java.io.*;
import java.util.*;
import java.util.regex.*;

import util.gen.*;

/**Parses meme output into MemeMotif objects and stores some generalize info re the meme run.
     * Must supply the full command line to fire meme, use full paths for meme and the seq file.*/
public class MemeParser {
    
    //fields
    private String cmdLn;       //one line command used to fire the meme program, with all the -switches
    private String parsedFileName;  //file text stripped of any path info
    private int numMotifs = 0;  //number of motifs pulled by meme
    private MemeMotif[] motifs; //array containing MemeMotif objects
    
    //constructor
    /**For single processor meme*/
    public MemeParser(String cndLine) {
        cmdLn = cndLine;
        parsedFileName = UtilMeme.extractFileName(cmdLn);
        String[] output = fireMeme(cmdLn);  // single processor meme runs
        parseMemeText(output);
    }
    /**For parallel processor meme on Sapo*/
    public MemeParser(String fullPathMeme, String fullPathFile, String memeParams){
    	File dataFile = new File(fullPathFile);
    	parsedFileName = dataFile.getName();  //to get text of multiFASTA file processed by MEME
    	cmdLn = fullPathMeme+ " "+fullPathFile+" "+memeParams;
    	String[] output = fireMemeParallel(fullPathMeme, fullPathFile, memeParams,
    		dataFile.getParent());  //note temp dir is the same as the data directory
		parseMemeText(output);
    }
    //primary methods
	/**takes DNA seq and reverse comps it ambiguous symbols OK.
	Will warn if it finds an unrecognized base.  Works with ' GATCRYWSKMBDHVNX .- '
	upper or lower case*/
     public static String reverseCompDNA(String seq){
        int seqLen = seq.length();
        StringBuffer rcSeq = new StringBuffer(seqLen);
        char test;
         for (int i=seqLen-1; i>=0; i--){
                test = seq.charAt(i);
                switch (test){
                    case 'a': rcSeq.append('t'); break;
                    case 'A': rcSeq.append('T'); break;
                    case 'c': rcSeq.append('g'); break;
                    case 'C': rcSeq.append('G'); break;
                    case 'g': rcSeq.append('c'); break;
                    case 'G': rcSeq.append('C'); break;
                    case 't': rcSeq.append('a'); break;
                    case 'T': rcSeq.append('A'); break;
                    case 'n': 
                    case 'N': 
                    case 'x': 
                    case 'X': 
                    case '-': 
                    case '.': 
                    case ' ':
                    case 'S': 
                    case 's':
                    case 'w':    
                    case 'W': rcSeq.append(test); break;
                    case 'r': rcSeq.append('y'); break;
                    case 'R': rcSeq.append('Y'); break;
                    case 'y': rcSeq.append('r'); break;
                    case 'Y': rcSeq.append('R'); break;
                    case 'm': rcSeq.append('k'); break;
                    case 'M': rcSeq.append('K'); break;
                    case 'k': rcSeq.append('m'); break;
                    case 'K': rcSeq.append('M'); break;
                    case 'b': rcSeq.append('v'); break;
                    case 'B': rcSeq.append('V'); break;
                    case 'd': rcSeq.append('h'); break;
                    case 'D': rcSeq.append('H'); break;
                    case 'h': rcSeq.append('d'); break;
                    case 'H': rcSeq.append('D'); break;
                    case 'v': rcSeq.append('b'); break;
                    case 'V': rcSeq.append('B'); break;
                    default: rcSeq.append (test); System.out.println("\nWarning: odd base in revComp-> '"+test+
                        "' Reverse Complement possibly incorrect!\n");
                }
            }
        return rcSeq.toString();
        }
    
    public String toString(){
        StringBuffer datadump = new StringBuffer(
        "Command Line: "+cmdLn+"\n"+
        "Number of motifs identified: "+numMotifs+"\n\n");
        if (numMotifs == 0) return datadump.toString();
        for (int i=0; i<motifs.length; i++) datadump.append(motifs[i].toString());
        return datadump.toString();
    }
    
    public void parseMemeText(String[] data){
        ArrayList alMotifs = new ArrayList();
        int len = data.length;
        //patterns
        Pattern sum = Pattern.compile("^MOTIF.+width.+sites.+");
        Pattern patHit = Pattern.compile("[GATCNX\\.]+\\s(\\S+)");
        Pattern patPSPM = Pattern.compile("\\s*(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s*");
        //run through data
        for (int i=60; i<len; i++){
            Matcher mat = sum.matcher(data[i]);
            if (mat.matches()){
                boolean go = true;
                //new MemeMotif!
                numMotifs++;
                MemeMotif mm = new MemeMotif();
                //set number
                mm.setNumberName(numMotifs);
                //set Summary line
                mm.setMotifSumLn(data[i]);
                //set motif description
                StringBuffer sb = new StringBuffer();
                i+=5;
                while (go){
                    if (data[i].startsWith("-")) go=false;
                    else sb.append(data[i]+"\n");
                    i++;
                }
                mm.setMotifDesc(sb.toString());
                sb=null; //needed for some reason
                //set motif sites
                go=true;
                i+=4;
                sb= new StringBuffer();
                ArrayList hits = new ArrayList();
                while (go){
                    if (data[i].startsWith("---------------------")) go = false; //need long -- to avoid short -- line
                    else {
                        sb.append(data[i]+"\n");
                        //set hits
                        Matcher mat2 = patHit.matcher(data[i]);
                        if (mat2.find()) hits.add(mat2.group(1));
                    }
                    i++;
                }
                mm.setMotifSites(sb.toString());
                //set hits
                String[] x =new String[hits.size()];
                hits.toArray(x);
                mm.setHits(x);
                i+=24;
                //set memePSPM
                ArrayList pspm = new ArrayList();
                go=true;
                while (go){
                    if (data[i].startsWith("letter-probability matrix")){
                        i++;
                        while (go){
                            Matcher mat3 = patPSPM.matcher(data[i]);
                            if (mat3.matches()) {
                                double[] ACGT = {Double.parseDouble(mat3.group(1)), Double.parseDouble(mat3.group(2)),
                                Double.parseDouble(mat3.group(3)),Double.parseDouble(mat3.group(4))};
                                pspm.add(ACGT);
                                i++;
                            }
                            else go = false;
                        }
                    }
                    i++;
                }
                double[][] memePSPM = new double[pspm.size()][4];
                for (int j=0; j<memePSPM.length; j++){
                    memePSPM[j] = (double[])pspm.get(j);
                }
                mm.setMemePSPM(memePSPM); 
                i+=5;
                //set additional params
                mm.setMemeParserRef(this);
                alMotifs.add(mm);
            }
        }
        //set array of motifs in MemeParser
        if (numMotifs >0){
            motifs = new MemeMotif[alMotifs.size()];
            alMotifs.toArray(motifs);
        }
    }
    
    public String[] fireMeme(String aCmd){
        //method to launch meme, need to have -nostatus option to get straight output
        ArrayList dataArrayList = new ArrayList(1000);
        try {
            Runtime rt = Runtime.getRuntime();
            Process p = rt.exec(aCmd);
            BufferedReader data = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String line;
            while ((line = data.readLine()) != null){
                dataArrayList.add(line);
                //System.out.println(line);
            }
            data.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        dataArrayList.trimToSize();
        String[] x = new String[dataArrayList.size()];
        dataArrayList.toArray(x);
        return x;
    }
 
	/**Fires meme using a shell script and qsub, checks every 15 sec to see if done, returns unparsed MEME output.
	 * To be run on Sapo cluster.
	 * @param fullPathMeme full path to the MEME program /home/sapo/software/seqanal/motifs/meme/meme.3.0.4/bin/meme
	 * @param fullPathFile full path to the multi FASTA for MEME to run on /home/sapo/nix/parMeme/kr4truk
	 * @param memeParams flags and options for the meme program -nmotifs 10 -evt 0.1 -minw 4 -maxw 12 -mod zoops -dna -revcomp -text -nostatus
	 * @param fullPathToTempDir full path to a directory where temporary files will be written and then deleted.*/
	public String[] fireMemeParallel(String fullPathMeme, String fullPathFile, String memeParams, String fullPathToTempDir){
		//create shell script wrapper file
		String command = fullPathMeme+" "+ fullPathFile+" -p $NSLOTS "+memeParams;
		String name = "MemeR"+new Date().getTime();
		String shell = "#!/bin/sh \n#$ -N "+name+" \n#$ -pe lam 10-19 \n#$ -j y \n#$ -o "+fullPathToTempDir+" \n"+command+ " \necho \"Finnished!\"\n";
		File shellFile = new File(fullPathToTempDir,name +".sh");
		IO.writeString(shell, shellFile);

		//fire shell script
		try {
			Runtime rt = Runtime.getRuntime();
			Process p = rt.exec("qsub "+shellFile);
			System.out.println("Fired qsub");
		}
		catch (IOException e) {
			e.printStackTrace();
		}

		//Wait 5min then check for file and for Finnished!
		boolean wait = true;
		File dir = new File(fullPathToTempDir);
		File outputFile = null;
		long counter =0;
		long milSec = 15000;  //checking every 15 sec
		ArrayList dataArrayList = new ArrayList(1000);
		try{
			while (wait){
				Thread.sleep(milSec);
				//find file
					if (outputFile==null){				
						//fetch contents of directory
						String[] fileNames = dir.list();
						//run through each file
						for (int i= fileNames.length-1; i>=0; i--){
							if (fileNames[i].startsWith(name+".o")) { //the qsub output file MemeR1081900495429.o50151
								outputFile = new File(fullPathToTempDir, fileNames[i]);
								System.out.println("Found results file: "+outputFile.getName());
								break;
							}
						}
					}
					if (outputFile!=null){ //check if null again
						//check if the last line is "Finnished!"
						String lastLine = "";
						String line;
						dataArrayList.clear();
						BufferedReader in = new BufferedReader(new FileReader(outputFile));
						while ((line = in.readLine()) !=null) {
							lastLine= line;
							dataArrayList.add(line);
						}
						in.close();
						//System.out.println("\nLast line is : "+lastLine);			
						if (lastLine.startsWith("Finnished!")) wait=false;  //stop found exit wait loop
					
					}
					
				//put break in so thread will stop after a very long time.
				counter++;	
				System.out.print(".");

				if (counter>57600000) { //16hrs 57600000
					System.out.println("\n    Error: MemeR shell script failed to return from qsub after "+counter/60000+
						" minutes.\nFind and kill the job: (text MemeR, process number is the last digits after the .o in -> "+outputFile.getName());
					System.exit(1);	
				}
				
			}
		} catch (InterruptedException e){
			e.printStackTrace();
		} catch (IOException ex){
			ex.printStackTrace();
		}
		System.out.println("Results are ready!");
		
		//delete files
		String[] fileNames = dir.list();
		//run through each file
		for (int i= fileNames.length-1; i>=0; i--){
			if (fileNames[i].startsWith(name)) { //the qsub output file MemeR1081900495429.o50151
				 new File(fullPathToTempDir, fileNames[i]).delete();
				System.out.println("Deleting-> "+fileNames[i]);
			}
		}
		System.out.println("Time for run: "+(counter*(milSec/1000))+" seconds\n");
		//return meme data
		dataArrayList.trimToSize();
		String[] x = new String[dataArrayList.size()];
		dataArrayList.toArray(x);
		return x;
	}
 

 
 
     //getter methods
    public String getCmdLn(){return cmdLn;}
    public int getNumMotifs(){return numMotifs;}
    public MemeMotif[] getMemeMotifs(){return motifs;}
    public String getParsedFileName(){return parsedFileName;}  
}



