package util.bio.annotation;
import java.util.*;
import java.io.*;
/**
 * A TransGroup describes a product derived from a GeneGroup contains one transcript, possibly a translation, and a possibly a set of exons.
 */
public class TransGroup implements Serializable{
	//Fields
	Transcript transcript;
	Translation translation;			//might be null
	ExonIntron[] exons;					//might be null
	int[][] extractedStartEndExons;	//might be null
	
	public TransGroup (Transcript transcript, Translation translation, ExonIntron[] exons){
		this.transcript=transcript;
		this.translation=translation;
		this.exons=exons;
		//sort exons
		if (exons!=null){
			Arrays.sort(exons);
			extractedStartEndExons = GeneGroup.extractStartEndsInts(exons);
		}
		//set TransGroup reference which also fires the Translation.makeSegments() method
		if(translation!=null) translation.setTransGrp(this);
	}
	public String toString(){
		StringBuffer x = new StringBuffer();
		if (transcript!=null){
			x.append("TransGroup: ");
			x.append(transcript.getName());
			x.append("\n  ");
			x.append(transcript.toString());
			x.append("\n");
		}
		if (translation!=null) x.append("  "+translation.toString()+"\n");
		if (exons !=null){
			for (int i=0; i<exons.length; i++){
				x.append("  "+exons[i].toString()+"\n");
			}
		}
		return x.toString();
	}
	
	public ExonIntron[] getExons(){return exons;}
	public Transcript getTranscript(){return transcript;}
	public Translation getTranslation(){return translation;}
	public int[][] getExtractedExons(){
		return extractedStartEndExons;
	}
	
	
}
