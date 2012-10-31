package util.bio.parsers;
import java.io.*;
import java.util.regex.*;

import util.gen.*;

public class RipMGIOrthologyResults {

	public static void main(String[] args) {

		//load lines
		File input = new File(args[0]);
		String[] lines = IO.loadFileIntoStringArray(input);

		//make patterns
		Pattern noResults = Pattern.compile("0 matching items");
		Pattern start = Pattern.compile("<TD CLASS=\"data\\d+leftMiddle\">human</TD>");
		Pattern end = Pattern.compile("<TD COLSPAN=\"4\" CLASS=\"data\\d+rightMiddle\"><B>References</B></TD>");
		Pattern source = Pattern.compile("<TD CLASS=\"data\\d+centerMiddle\">");
		Pattern tr = Pattern.compile("</*TR>");
		Pattern subTableStart = Pattern.compile("<TD CLASS=\"data\\dleftMiddle\"><TABLE BORDER=\"0\" CELLPADDING=\"0\" CELLSPACING=\"0\">");
		Pattern genBank = Pattern.compile("NM_\\d+");
		Pattern subTableEnd = Pattern.compile("</TABLE></TD>");
		Pattern mouseBuild = Pattern.compile("NCBI Mouse Build");


		//for each line
		boolean in = false;
		for (int i=0; i<lines.length; i++){
			//look for noResults
			Matcher mat = noResults.matcher(lines[i]);
			if (mat.find()) Misc.printExit("<TR><TD>No ortho<TD>"+Misc.removeExtension(input.getName())+"</TR>");
			//look for start
			mat = start.matcher(lines[i]);
			if (mat.find()){
				System.out.println("<TR>");
				//while within record
				in = true;
				while (in){

					//skip source and mouse build
					Matcher build = mouseBuild.matcher(lines[i]);
					mat = source.matcher(lines[i]);
					if (mat.find() || build.find()) {}
					else {
						//stop?
						mat = end.matcher(lines[i]);
						if (mat.find()) {
							in = false;
							System.out.println("</TR>");
						}
						else {
							//look for adjacent trs
							mat = tr.matcher(lines[i]);
							if (mat.matches()){
								mat = tr.matcher(lines[i+1]);
								if (mat.matches()) i++;
							}
							else{
								//look for sub table
								mat = subTableStart.matcher(lines[i]);
								if (mat.matches()){
									//within sub table
									System.out.print("<TD>");
									boolean go = true;
									while (go){
										//stop of table?
										mat = subTableEnd.matcher(lines[i]);
										if (mat.matches()) {
											go = false;
										}
										else {//GenBankRef?
											mat = genBank.matcher(lines[i]);
											if (mat.find()) {
												System.out.print(mat.group()+" ");
											}
											i++;
										}
									}
									System.out.println("</TD>");
								}
								else System.out.println(lines[i]);
							}

						}
					}
					i++;
				}
			}
		}
	}

}
