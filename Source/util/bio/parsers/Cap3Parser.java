package util.bio.parsers;
import java.io.*;
import java.util.*;
import util.gen.*;

//NOT FUNCTIONAL
public class Cap3Parser {

	public Cap3Parser (String[] args){
		File data = new File (args[0]);
		System.out.println("\nParsing "+data.getName()+"\n");
		parseAndPrint(data);
	}

	public void parseAndPrint (File file){
		try {
			BufferedReader in = new BufferedReader( new FileReader (file));
			//skip header
			skipHeader(in);
			//parse out contig
			String line;
			ArrayList lines = new ArrayList();
			while ((line=in.readLine()) !=null){
				if (line.startsWith("*******************")){
					System.out.println(line);
					parseAndPrintContig(lines);
					lines = new ArrayList();
				}
				else lines.add(line);
			}
			//print last one
			parseAndPrintContig(lines);
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void parseAndPrintContig(ArrayList lines) {
		int num = lines.size();
		LinkedHashSet ests = new LinkedHashSet();
		for (int i=0; i<num; i++){
			String line = (String)lines.get(i);
			String[] tokens = line.split("\\s+");
		}
	}
	
	public void skipHeader (BufferedReader in) throws Exception {
			String line;
			boolean go = false;
			while ((line=in.readLine()) !=null){
				if (line.startsWith("DETAILED DISPLAY OF CONTIGS")){
					go = true;
					break;
				}
			}
			if (go == false) Misc.printExit("\nNo contigs found?\n");
			in.close();
	}

	
	
	public static void main(String[] args) {
		new Cap3Parser (args);
	}

}
