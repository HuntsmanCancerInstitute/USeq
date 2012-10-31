package gata.aligner;

import java.io.*;

/**
 * @author nix
 *
 * A class to represent alignment parameters entered for a particular alignment.  This is
 * saved and made availible to GATAPlotter.
 */
public class AlignParams implements Serializable {
    
    //defaults
    private int WIN_SIZE =24 ;      // window size (24)
    private int MATCH = 5;          //bonus for a match (5)
    private int MISMATCH =-4;       //penalty for a mismatch (-4)
    private int GAP_CREATE= -10;    //gap creation penalty (-10)
    private int GAP_EXT =-4;        //gap extention penalty (-4)
    private int MIN_SCORE= 80;      // will only align those seqs with scores higher that (80)
    private int START_INDEX_REFSEQ = 0; //where the first base is 0, yet in the output its 1
    private int START_INDEX_COMPSEQ =0 ;//where the first base is 0
    private int LENGTH_REFSEQ;
    private int LENGTH_COMPSEQ;
    private String refSeqFile;
	private String compSeqFile;
    private String NAME_REFSEQ;
    private String NAME_COMPSEQ;
    private String bl2seq; //path to Blast bl2seq program
    private String baseName; //text for saving alignment and param object files
    private String PATH_RESULTS; //will default to where DNAFILE comes from, path to where user wants results placed
    private boolean DUST = false;
    private boolean EXTRACT = true;
    private boolean compVisible = true;
    private double LAMBDA;
    private double K;
    private double H;
    private double EFF_M;
    private double EFF_N;
    
    public AlignParams(){} //needed for retreiving object from disk;
    public AlignParams(AlignerPreferences ap) {
        WIN_SIZE = Integer.parseInt(ap.getWindow());
        MATCH = ap.getMatch().intValue(); 
        MISMATCH = ap.getMisMatch().intValue(); 
        GAP_CREATE = ap.getCreate().intValue(); 
        GAP_EXT = ap.getExtend().intValue();
        MIN_SCORE = (int)Math.round(Double.parseDouble(ap.getScore()));  
        refSeqFile= ap.getRefSeq();
        compSeqFile = ap.getCompSeq();
        bl2seq = ap.getBl2SeqProg();
        baseName = ap.getBaseName();
		PATH_RESULTS = ap.getStorageLoc(); 
        START_INDEX_REFSEQ = (int)Double.parseDouble(ap.getRefStart()); 
        START_INDEX_COMPSEQ = (int)Double.parseDouble(ap.getCompStart());
        if (ap.getMask().equals("Yes")) DUST = true;
        if (ap.getExtract().equals("No")) EXTRACT = false;
    }
    
    public int getWIN_SIZE(){ return WIN_SIZE;}
    public int getMATCH(){ return MATCH;}
    public int getMISMATCH(){ return MISMATCH;}
    public int getGAP_CREATE(){ return GAP_CREATE;}
    public int getGAP_EXT(){ return GAP_EXT;}
    public int getMIN_SCORE(){ return MIN_SCORE;}
    public int getMAX_SCORE(){ return MATCH * WIN_SIZE;}
    public int getSTART_INDEX_REFSEQ(){ return START_INDEX_REFSEQ;}
    public int getSTART_INDEX_COMPSEQ(){ return START_INDEX_COMPSEQ;}
    public int getLENGTH_REFSEQ(){ return LENGTH_REFSEQ;}
    public int getLENGTH_COMPSEQ(){ return LENGTH_COMPSEQ;}
    public String[] getDNAFILES(){ return new String[] {refSeqFile, compSeqFile};}
    public String getNAME_REFSEQ(){ return NAME_REFSEQ;}
    public String getNAME_COMPSEQ(){ return NAME_COMPSEQ;}
    public String getPathToResults() { return PATH_RESULTS;}
    public boolean getDUST() {return DUST;}
    
    public void setLENGTH_REFSEQ(int l){ LENGTH_REFSEQ = l;}
    public void setLENGTH_COMPSEQ(int l){ LENGTH_COMPSEQ = l;}
    public void setNAME_REFSEQ(String s){ NAME_REFSEQ = s;}
    public void setNAME_COMPSEQ(String s){ NAME_COMPSEQ = s;}
    
    public String toHTMLString(){
    	String alignerFile = "";
		try{ alignerFile = (new File(PATH_RESULTS,baseName)).getCanonicalPath();}
		catch (IOException e){e.printStackTrace();}
    	return
    	"<body bgcolor='#FFFFFF' link='#0000FF' font face='arial'>"
    	+"Reference Seq File = "+refSeqFile
    	+"<br>Comparative Seq File = "+compSeqFile
    	+"<br>GATAligner File = "+alignerFile
		+"<br>Name Reference Sequence = "+NAME_REFSEQ
		+"<br>Name Comparative Sequence = "+NAME_COMPSEQ
    	+"<br>"
        +"<br>Window Size = "+WIN_SIZE
        +"<br>Match = "+MATCH
        +"<br>MisMatch = "+MISMATCH
        +"<br>Gap Creation = "+GAP_CREATE
        +"<br>Gap Extension = "+GAP_EXT
        +"<br>Raw Score Cut Off = "+MIN_SCORE
		+"<br>Mask Low Complexity Sequences = "+DUST
		+"<br>"
		+"<br>Extract Start Indexes from FASTA  = "+EXTRACT
        +"<br>Start Index Reference Sequence = "+START_INDEX_REFSEQ
        +"<br>Start Index Comparative Sequence = "+START_INDEX_COMPSEQ
        +"<br>Length Reference Sequence = "+LENGTH_REFSEQ
        +"<br>Length Comparative Sequence = "+LENGTH_COMPSEQ
        +"</body>"
        
        ;
    }
    public double getEFF_M() {
		return EFF_M;
	}
	public double getEFF_N() {
		return EFF_N;
	}
	public double getH() {
		return H;
	}
	public double getK() {
		return K;
	}
	public double getLAMBDA() {
		return LAMBDA;
	}
	public void setEFF_M(double d) {
		EFF_M = d;
	}
	public void setEFF_N(double d) {
		EFF_N = d;
	}
	public void setH(double d) {
		H = d;
	}
	public void setK(double d) {
		K = d;
	}
	public void setLAMBDA(double d) {
		LAMBDA = d;
	}
	public String getBaseName() {
		return baseName;
	}
	public String getBl2seq() {
		return bl2seq;
	}
	public String getCompSeqFile() {
		return compSeqFile;
	}
	public String getRefSeqFile() {
		return refSeqFile;
	}

	public void setSTART_INDEX_COMPSEQ(int start_index_compseq) {
		START_INDEX_COMPSEQ = start_index_compseq;
	}
	public void setSTART_INDEX_REFSEQ(int start_index_refseq) {
		START_INDEX_REFSEQ = start_index_refseq;
	}
	public boolean isEXTRACT() {
		return EXTRACT;
	}
	public void setEXTRACT(boolean extract) {
		EXTRACT = extract;
	}
	public boolean isCompVisible() {
		return compVisible;
	}
	public void setCompVisible(boolean compVisible) {
		this.compVisible = compVisible;
	}
}