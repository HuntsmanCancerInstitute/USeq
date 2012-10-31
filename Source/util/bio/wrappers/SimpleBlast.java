package util.bio.wrappers;

import java.util.*;
import java.io.*;
import java.util.regex.*;

/**
 * Wrapper for bl2seq blast.
 * */
public class SimpleBlast {
	String bl2seq; //full path filename to bl2seqprogram
    
    public SimpleBlast(String bl2seq){
    	this.bl2seq = bl2seq;
    }
    
    public String[] blastIt(String seq1File, String seq2File){
        //uses BLAST bl2seq program to align two sequences found in the seq files, won't work on windows
        
		//convert to String[], cannot use a text since exec uses as stringTokenizer to bust up command, if spaces exist in the folder or  file names then this creates havoc
		String[] commandArray ={bl2seq,"-p","blastn","-i", seq1File,"-j", seq2File,"-e","0.1","-W","7"};
                               
        ArrayList dataArrayList = new ArrayList(500);
        try {
        	//System.out.println("launching blast...");
            Runtime rt = Runtime.getRuntime();
            //rt.traceInstructions(true); //for debugging
            //rt.traceMethodCalls(true); //for debugging
            Process p = rt.exec(commandArray);
            BufferedReader data = new BufferedReader(new InputStreamReader(p.getInputStream()));
			//BufferedReader data = new BufferedReader(new InputStreamReader(p.getErrorStream())); //for debugging
            String line;
            while ((line = data.readLine()) != null){
                dataArrayList.add(line);
                //System.out.println("X: "+line);
            }
            data.close();   //this is close/null stuff is needed to invoke the garbage collector
            p.waitFor();
            p=null;
            data = null;
            rt = null;
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        dataArrayList.trimToSize();
        String[] x = new String[dataArrayList.size()];
        dataArrayList.toArray(x);
        return x;
    }
    
    public static ArrayList parseBlast(String[] data) {
        //this takes an array of BLAST bl2seq results and parses the info into an array
        //note, if the seq returned is too big, must combine it with previous line
        // 0=score, 1=start, 2=stop, 3=seq, 4=start, 5=stop, 6=seq, 7=orientation (0 for +/+, 1 for +/-);
        
        //make some precompiled matchers
        Pattern patScore = Pattern.compile("Score = .+\\((\\d+)\\)");
        Pattern patSSS = Pattern.compile("(\\d+)\\s+([a-zA-Z-]+)\\s+(\\d+)");
        Pattern patLambda = Pattern.compile("^Lambda.+");

        int size = data.length;
        ArrayList parsedData = new ArrayList(size/8);
        try {
            for (int i=6; i< size; i++){
                Matcher w = patScore.matcher(data[i]);
                Matcher l = patLambda.matcher(data[i]);
                
                if (w.find()){
                    //get score
                    int score = new Integer((w.group(1))).intValue();
                    //get orientation
                    i+=2;
                    String o = data[i].substring(17,18);
                    int ori;
                    if (o.equals("P")) ori =0;
                    else ori=1;
                    
                    //get everything else
                    i+=3;
                    String seq1Ln = data[i];
                    i+=2;
                    String seq2Ln = data[i];
                    
                    //extract start stop seq
                    ArrayList datSeq1 = parseLine(seq1Ln, patSSS);
                    ArrayList datSeq2 = parseLine(seq2Ln, patSSS);
                    
                    //look ahead to see if there are any more lines
                    boolean inSeq = true;
                    while (inSeq){
                        i+=3;
                        if (data[i].startsWith("Query")){
                            //get additional lines
                            ArrayList new1 = parseLine(data[i], patSSS);
                            i+=2;
                            ArrayList new2 = parseLine(data[i], patSSS);

                            //join with previous seqs and change stop numbers
                            datSeq1.set(2,((String)datSeq1.get(2)).concat((String)new1.get(2)));
                            datSeq2.set(2,((String)datSeq2.get(2)).concat((String)new2.get(2)));
                            datSeq1.set(1, new1.get(1));
                            datSeq2.set(1, new2.get(1));
                        }
                        else inSeq = false;
                    }
                    //add data to arraylist 0=score, 1=start, 2=stop, 3=seq, 4=start, 5=stop, 6=seq, 7=orientation (0 for +/+, 1 for +/-);
                    parsedData.add(new Integer(score));
                    parsedData.addAll(datSeq1);
                    parsedData.addAll(datSeq2);
                    parsedData.add(new Integer(ori));
                }
            }
        } catch (NumberFormatException n){
            n.printStackTrace();
        }
        parsedData.trimToSize();
        return parsedData;
    }
    
    public static ArrayList parseLine(String line, Pattern pat){
        //returns start, stop, sequence
        ArrayList dat = new ArrayList(3);
        Matcher mat = pat.matcher(line);
        mat.find();
        dat.add(new Integer(mat.group(1)));
        dat.add(new Integer (mat.group(3)));
        dat.add(mat.group(2));      
        return dat;        
    }
    
    /**Parses blast output from blastall for dna into BlastResult objects.*/
	public static BlastResult[] parseBlastGeneral(String[] data) {
		String line;
		ArrayList blastResultsAL = new ArrayList();
		for (int x=0; x< data.length; x++){
			line = data[x];
			//skip to first >
			if (line.startsWith(">")){
				boolean go = true;
				while (go){
					BlastResult res = new BlastResult();
					//parse match id and text assumes divided by a :
					String[] idName = line.substring(1).split("\\s*:\\s*");
					if (idName.length==2){
						res.setMatchId(idName[0]);
						res.setMatchName(idName[1]);
					}
					else res.setMatchName(idName[0]);
					
					//parse match length
					line = data[++x];
					res.setMatchLength(line.replaceAll("\\D", ""));
					x+=2;
					// parse score and expect
					line = data[x];
					String[] scoreExpect = line.split(",");
					res.setScore(scoreExpect[0].split("\\s*=\\s*")[1]);
					res.setExpect(scoreExpect[1].split("\\s*=\\s*")[1]);
					//parse identities and gaps
					line = data[++x];
					String[] identitiesGaps = line.split(",");
					res.setIdentities(identitiesGaps[0].split("\\s*=\\s*")[1]);
					if (identitiesGaps.length==2) res.setGaps(identitiesGaps[1].split("\\s*=\\s*")[1]);
					//parse strand
					line = data[++x];
					res.setStrand(line.split("\\s*=\\s*")[1]);
					x+=2;
					//parse alignment
					ArrayList al = new ArrayList();
					boolean parse = true;
					while (parse){
						line = data[++x];
						boolean newMatch = line.startsWith(">");
						boolean end = line.indexOf("Database:")!=-1;
						if (newMatch == false && end == false) al.add(line);
						else {
							parse = false;
							if (end) go = false;
						}
					}
					//convert aligment to a String[] skip last two blank lines
					int num = al.size()-2;
					String[] alignment = new String[num];
					for (int i=0; i<num; i++) alignment[i] = (String)al.get(i);
					res.setAlignment(alignment);
					blastResultsAL.add(res);
				}
			}
		}
		BlastResult[] parsed = new BlastResult[blastResultsAL.size()];
		blastResultsAL.toArray(parsed);
		return parsed;
	}
}
