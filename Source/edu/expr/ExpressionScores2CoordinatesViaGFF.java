package edu.expr;
import java.io.*;
import util.gen.*;
import java.util.*;
import util.bio.parsers.gff.*;

public class ExpressionScores2CoordinatesViaGFF {
	
	public static void main (String[] args){
		//load gff3 file
		Gff3Parser gffParser = new Gff3Parser();
		gffParser.setRegExTypes(".+");
		gffParser.parseIt(new File ("/Users/nix/HCI/PIs/Mango/GFF/genesWS175.gff3"));
		
		//make hash of gff on id
		Gff3Feature[] gffLines = gffParser.getFeatures();
		HashMap hash = new HashMap (gffLines.length);
		for (int i=0; i< gffLines.length; i++) hash.put(gffLines[i].getId(), gffLines[i]);
		
		//run through text, value file lines
		try {
			File nameValueFile = new File ("/Users/nix/HCI/PIs/Mango/GenLists/Dump/dataRatio.txt");
			BufferedReader in = new BufferedReader (new FileReader (nameValueFile));
			
			//load header of text, score file
			String line = in.readLine();
			String[] dataName = line.split("\\t");
			
			//set booleans for file score type
			boolean log2Transform = false;
			boolean log10Transform = true;
			if (dataName[1].endsWith("ratio")){
				log2Transform = true;
				log10Transform = false;
			}
			
			//make files and writers, skipping index 0
			File[] files = new File[dataName.length];
			PrintWriter[] outs =  new PrintWriter[dataName.length]; 
			for (int i=1; i< dataName.length; i++) {
				files[i] = new File (nameValueFile.getParentFile(), dataName[i]+".egr");
				outs[i] = new PrintWriter ( new FileWriter (files[i]));
				//print header
				outs[i].print("# genome_version = C_elegans_WS175\n# score0 = ");
				if (log2Transform) outs[i].println("log2("+dataName[i]+")\n");
				else if (log10Transform) outs[i].println("-10log10("+dataName[i]+")\n");
				else outs[i].println(dataName[i]+"\n");
			}
			
			//for each line in text, value file
			while ((line = in.readLine())!= null){
				String[] nameScores = line.split("\\t");
				//get gff feature
				Gff3Feature gff = (Gff3Feature) hash.get(nameScores[0]);
				if (gff == null) {
//System.out.println("ID not found in gff file?! "+nameScores[0]+" line "+line);
					continue;
				}
				String header = gff.getSeqId()+"\t"+gff.getStart()+"\t"+gff.getEnd()+"\t"+gff.getStrand()+"\t";
				String head2 = gff.getId()+"\t"+gff.getSeqId()+"\t"+gff.getStart()+"\t"+gff.getEnd();
				StringBuffer head2SB = new StringBuffer();
				boolean allZeros = true;
				System.out.println("# genome_version = C_elegans_WS175\n# score0 = coor");
				//for each score print those not == 0
				for (int i=1; i< dataName.length; i++){
					head2SB.append("\t");
					head2SB.append(nameScores[i]);
					
					if (nameScores[i].equals("0") == false){
						allZeros = false;
						//transform score?
						double score = Double.parseDouble(nameScores[i]);
						if (log2Transform){
							score = Num.log2(score);
						}
						else if (log10Transform){
							score = Num.minus10log10(score);
						}
						//print line
						outs[i].println(header+score);
					}
				}
				if (allZeros == false) System.out.println(head2+head2SB);
			}
			
			//close handles and zip files
			for (int i=1; i< dataName.length; i++){
				outs[i].close();
				IO.zipAndDelete(files[i]);
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		
	}
	
}
