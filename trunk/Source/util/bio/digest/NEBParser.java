package util.bio.digest;
import java.util.regex.*;
import java.io.*;
import java.util.ArrayList;
/**
 * Parses a New England Biolabs restriction enzyme file into an array of {@link Enzyme}s.
 */
public class NEBParser {
	//Fields
	private String header;
	private Enzyme[] enzymes;
	
	public String getHeader (){			//must be called after the makeEnzymeArray method or returns null
		return header;
		}
        
	public Enzyme[] makeEnzymeArray(String gtype2c, int minLengthRecogSeq) {
		 ArrayList enzymesAL = new ArrayList (400);	//create new ArrayList
        
		Pattern patOne = Pattern.compile("<1>(\\S+)");	//compile patterns to parse enzy file
		Pattern patTwo = Pattern.compile("<2>(\\S+)");
		Pattern patThree = Pattern.compile("<3>(\\S+)");
		Pattern patFour = Pattern.compile("<4>(\\S+)");
		Pattern patFive = Pattern.compile("<5>(\\S+)");
        
		String line;	//initialize variables
		int Counter = 1;
		String name = "none";
		String iso = "none";
		String seq = "none";
		String meth = "none";
		String source = "none";

		try {	//Read in enzyme file line by line
			BufferedReader in = new BufferedReader (new FileReader ( gtype2c ));
			StringBuffer sb = new StringBuffer(); 
			Enzyme enzyme;
			while ((line = in.readLine()) !=null) {
				if (Counter < 41){ 	//parse header
					Counter++;
					sb.append(line + "\n");
					}					
				else {
					Matcher x = patOne.matcher ( line );
					if (x.lookingAt()) {
						name = x.group(1);
						continue;
						}
					 x = patTwo.matcher (line);
					 if (x.lookingAt()){
						iso = x.group(1);
						continue;
						}
					 x= patThree.matcher (line); 
					 if (x.lookingAt()){      
						seq = x.group(1);
						continue;   
						}
					 x= patFour.matcher (line); 
					 if (x.lookingAt()){      
						meth = x.group(1);
						continue;   
						}
					 x= patFive.matcher (line); 
					 if (x.lookingAt()){      
						source = x.group(1);
						enzyme = new Enzyme (name, seq, meth, source, iso);
						if (enzyme.getStrippedRecogSeq().replaceAll("[^GATCgatc]","").length()
								>=minLengthRecogSeq) enzymesAL.add(enzyme);		//add new Enzyme
						iso = "none";    //reset possible blank lines            
						meth = "none";
						source = "none";
						}
					}
				}
			header = sb.toString();
			in.close();
			}
		catch (IOException e) {
			e.printStackTrace();
			}
	enzymes = new Enzyme[enzymesAL.size()];		//convert ArrayList to a regular array
	enzymesAL.toArray(enzymes);
	return enzymes;
	}
	
	public Enzyme[] getEnzymes() {
		return enzymes;
	}

}    
    