package util.apps;
import java.io.*;
import util.gen.*;

/**Takes the output of RocWindowScanner from the cluster and renames the file */
public class RenameRWSResults {
	
	
	public static void main(String[] args) {
		File[] files = IO.extractFiles(new File(args[0]));
		for (int i=0; i< files.length; i++){
			String[] lines = IO.loadFile(files[i]);
			int minOligos = -1;
			String name = null;
			for (int j=0; j< lines.length; j++){
				if (lines[j].startsWith("Min Number Probes")){
					String[] tokens = lines[j].split(":");
					minOligos = Integer.parseInt(tokens[1].trim());
				}
				else if (lines[j].startsWith("Processing")){
					String[] tokens = lines[j].split("/");
					name = tokens[tokens.length-1].replaceAll("\\.\\.\\.","");
					break;
				}
			}
			if (minOligos == -1) {
				System.out.println("Error: could not parse a Min Number Probes! Skipping -> "+files[i]);
			}
			else {
				File newFile = new File (files[i].getParentFile(), name+"Min"+minOligos+".xls");
				if (newFile.exists()) Misc.printExit("Error: new file already exists?! "+files[i].getName()+" -> "+newFile.getName());
				files[i].renameTo(newFile);
				System.out.println(files[i].getName()+" -> "+newFile.getName());
			}
		}
		
	}
	
}
