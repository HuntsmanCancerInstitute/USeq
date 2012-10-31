package trans.misc;
import util.gen.*;
import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**Utility to convert simple one graph xxx.bar files to xxx.gr files.*/
public class Bar2Gr {
	
	private File[] files;
	
	public  Bar2Gr (String[] args) {
		processArgs(args);

		for (int x=0; x< files.length; x++){
			//attempt to parse xxx.bar file
			System.out.println("\tLoading... "+files[x].getName());
			GrGraph gr = Util.readSimpleGrBarFile(files[x]);
			if (gr == null) System.out.println("\t\tProblem parsing "+files[x]);
			else{
				File grFile = new File(files[x].getParent(), Misc.replaceEnd(files[x].getName(),".bar","") + 
						"_"+gr.getGenomeVersion()+"_"+gr.getChromosome()+".gr");
				System.out.println("\tWriting... "+grFile);
				gr.writeGrFile(grFile);
			}
			
		}
		
	}
	
	
	public static void main(String[] args) {
		if (args.length==0){
			printDocs();
			System.exit(0);
		}	
		new Bar2Gr(args);
	}		
	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String filesString = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': filesString = args[i+1]; i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nError: unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//parse files and genome version
		if (filesString == null ) Misc.printExit("\nError: cannot find your xxx.bar file(s)?");
		files = IO.extractFiles(new File(filesString), ".bar");
	}	
	
	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                 Bar2Gr: Nov 2006                                 **\n" +
				"**************************************************************************************\n" +
				"Converts xxx.bar to text xxx.gr files.\n\n" +

				"-f The full path directory/file text for your xxx.bar file(s).\n" +
				
				"\nExample: java -Xmx1500M -jar pathTo/T2/Apps/Bar2Gr -f /affy/BarFiles/ \n\n" + 
				
		"**************************************************************************************\n");		
	}
	
}
