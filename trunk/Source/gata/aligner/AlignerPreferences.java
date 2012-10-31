package gata.aligner;

import java.io.*;
/**
 * @author nix 
 *
 Container class to hold aligment preferences, this is saved after running GATAligner as a
 serialized object and used in subsequent GATAligner runs so the user doesn't have to reset
 the defaults every run.
 */
public class AlignerPreferences implements Serializable {

	//fields, these are the defaults to be overriden by user
	private String refSeq = "";
	private String compSeq = "";
	private String bl2SeqProg = "";
	private String storageLoc = "";
	private String baseName = "";
	private Integer match = new Integer(5);
	private Integer misMatch = new Integer(-4);
	private Integer create = new Integer(-10);
	private Integer extend = new Integer(-4);
	private String mask = "No";
	private String extract = "No";
	private String window = "24";
	private String score = "80";
	private String refStart = "1";
	private String compStart = "1";
	
	public String toString(){
		return
		"refSeq: "+refSeq
		+"\ncompSeq: "+compSeq
		+"\nbl2SeqProg: "+bl2SeqProg
		+"\nstorageLoc: "+storageLoc
		+"\nbaseName: "+baseName
		+"\nmatch: "+match
		+"\nmisMatch: "+misMatch
		+"\ncreate: "+create
		+"\nextend: "+extend
		+"\nmask: "+mask
		+"\nextract: "+extract
		+"\nwindow: "+window
		+"\nscore: "+score
		+"\nrefStart: "+refStart
		+"\ncompStart: "+compStart
		;
	}

	public void setPreferences(
		String refSeq,
		String compSeq,
		String bl2SeqProg,
		String storageLoc,
		String baseName,
		Integer match,
		Integer misMatch,
		Integer create,
		Integer extend,
		String mask,
		String extract,
		String window,
		String score,
		String refStart,
		String compStart) {
		this.refSeq=refSeq;
		this.compSeq = compSeq;
		this.bl2SeqProg = bl2SeqProg;
		this.storageLoc = storageLoc;
		this.baseName = baseName;
		this.match = match;
		this.misMatch = misMatch;
		this.create = create;
		this.extend = extend;
		this.mask = mask;
		this.extract = extract;
		this.window = window;
		this.score = score;
		this.refStart = refStart;
		this.compStart = compStart;
	}

	public void saveAlignerPreferences() {
		try {
			ObjectOutputStream out =
				new ObjectOutputStream(
					new FileOutputStream("GATAlignerPreferences"));
			out.writeObject(this);
			out.close();
		} catch (Exception e) {
			System.out.println("\n\nCannot Save AlignerPreferences\n\n");
			e.printStackTrace();
		}
	}

	public String getBl2SeqProg() {
		return bl2SeqProg;
	}
	public String getCompSeq() {
		return compSeq;
	}
	public Integer getCreate() {
		return create;
	}
	public Integer getExtend() {
		return extend;
	}
	public String getMask() {
		return mask;
	}
	public Integer getMatch() {
		return match;
	}
	public Integer getMisMatch() {
		return misMatch;
	}
	public String getRefSeq() {
		return refSeq;
	}
	public String getStorageLoc() {
		return storageLoc;
	}
	public String getCompStart() {
		return compStart;
	}
	public String getRefStart() {
		return refStart;
	}
	public String getScore() {
		return score;
	}
	public String getWindow() {
		return window;
	}
	public String getBaseName() {
		return baseName;
	}
	public void setBl2SeqProg(String string) {
		bl2SeqProg = string;
	}

	public String getExtract() {
		return extract;
	}
}
