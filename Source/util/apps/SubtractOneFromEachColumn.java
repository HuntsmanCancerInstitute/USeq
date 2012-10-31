package util.apps;
import java.io.*;

import util.gen.*;

/**Subtracts 1 from each column assumes integers. Tab delimited file or directory containing files to subtract.
 * Blank lines and # comments printed and skipped */
public class SubtractOneFromEachColumn {


	public static void main(String[] args) {
		try {
			File[] filesToParse = IO.extractFiles(new File (args[0]));
			
			//for each file
			for (int z =0; z<filesToParse.length; z++){
				File subFile = new File (filesToParse[z].getCanonicalPath()+".sub");
				BufferedReader in = new BufferedReader (new FileReader( filesToParse[z] ));
				PrintWriter out = new PrintWriter (new FileWriter(subFile));
				String line;
				while ((line = in.readLine()) != null){
					line = line.trim();
					//skip comments and blank lines
					if (line.startsWith("#") || line.length() ==0) {
						out.println(line);
						continue;
					}
					String[] columns = line.split("\\s");
					StringBuilder sb = new StringBuilder();
					int first = Integer.parseInt(columns[0]) - 1;
					sb.append(first);
					for (int i=1; i< columns.length; i++){
						int sub = Integer.parseInt(columns[i]) - 1;
						sb.append("\t");
						sb.append(sub);
					}
					out.println(sb);
				}
				in.close();
				out.close();
			}
		}
		catch (Exception e){
			e.printStackTrace();
		}

	}

}
