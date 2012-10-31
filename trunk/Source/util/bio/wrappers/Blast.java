package util.bio.wrappers;

/*import java.util.*;
import java.io.*;
import java.util.regex.*;
//import gata.*;
*/
/**
	Template to fire and fetch results from the ncbi bl2seq program.  Not really useful since it's part of GATA but good as a
	start.
 * */
public class Blast {}  /*
    public AlignParams ap;
    public Blast(AlignParams AP) {
        ap = AP;
   }
    public String[] blastIt(){
        //uses BLAST bl2seq program to align two sequences found in the seq files
        String[] dustOpt = {"-F","F"};
        if (ap.getDUST()) dustOpt = new String[]{"-F","T"};
        //check operating system and alter paths in windows
        String bl2seq = ap.getBl2seq();
        String ref = ap.getRefSeqFile();
        String comp = ap.getCompSeqFile();
        String os = System.getProperty("os.name").toLowerCase();
        if (os.matches(".*windows.*")){        	
        	bl2seq = "\""+bl2seq+ "\"";
			ref = "\""+ref+ "\"";
			comp = "\""+comp+ "\"";
        }
        
		//convert to String[], cannot use a text since exec uses as stringTokenizer to bust up command, if spaces exist in the folder or  file names then this creates havoc
		String[] commandArray ={
		bl2seq,
		"-p","blastn",
		"-i", ref, //seq1
		"-j", comp,    //seq2
		"-e","50",         //expectation score cut off
		"-q", Integer.toString(ap.getMISMATCH()),  //mismatch penalty (ie -4)
		"-r", Integer.toString(ap.getMATCH()),      //match reward (ie 5)
		"-G", Integer.toString(ap.getGAP_CREATE()*-1),  //gap creation penalty (ie 10), not neg!
		"-E", Integer.toString(ap.getGAP_EXT()*-1),     //gap extension penalty (ie 4), not neg!
		dustOpt[0],dustOpt[1],           //DUST repetative seq filtering
		"-W","7"};        //word size, only returns alignments with 7 consec identical matches
        
                               
        ArrayList dataArrayList = new ArrayList(2000);
        try {
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
    
    public ArrayList parseBlast(String[] data) {
        //this takes an array of BLAST bl2seq results and parses the info into an array
        //note, if the seq returned is too big, must combine it with previous line
        //0=score, 1=start, 2=stop, 3=seq, 4=orientation (0 for +/+, 1 for +/-);
        
        //make some precompiled matchers
        Pattern patScore = Pattern.compile("Score = .+\\((\\d+)\\)");
        Pattern patSSS = Pattern.compile("(\\d+)\\s+([a-zA-Z-]+)\\s+(\\d+)");
        Pattern patLambda = Pattern.compile("^Lambda.+");

        int size = data.length;
        ArrayList parsedData = new ArrayList(size/8);
        int scoreCutOff = ap.getMIN_SCORE();
        try {
            for (int i=6; i< size; i++){
                Matcher w = patScore.matcher(data[i]);
                Matcher l = patLambda.matcher(data[i]);
                
                if (w.find()){
                    //get score
                    int score = new Integer((w.group(1))).intValue();
                    if (score<scoreCutOff) continue;
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
				else if (l.find()){
				//get stats from blast report
				String[] goodies = data[i+1].trim().split("\\s+");			
				ap.setLAMBDA(Double.parseDouble(goodies[0]));
				double k =Double.parseDouble(goodies[1]);
				ap.setK(k);
				double H = Double.parseDouble(goodies[2]);
				ap.setH(H);
				//calculate eff m and n, blastn is getting funny between versions, the latest 2.2.8 has a bug in the eff n and m calculations so i'd deriving them
				int[] mn = GATAUtil.calculateEffectiveMandN((double)ap.getLENGTH_REFSEQ(), (double)ap.getLENGTH_COMPSEQ(),k, H);
				ap.setEFF_M(mn[0]);
				ap.setEFF_N(mn[1]);				
				break;
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
}
*/