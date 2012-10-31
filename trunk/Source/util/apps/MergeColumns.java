package util.apps;
import java.io.*;

import util.gen.*;

public class MergeColumns {


	/**
	 * Merges the columns of two files using a tab, assumes the same length!
	 */
	public static void main(String[] args) {
		if (args.length == 1) mergeColumns(new File (args[0]));
		else if (args.length !=2) Misc.printExit("\nEnter two files of equal length to merge with a tab.\n");
		else {
			try{

				File file1 = new File (args[0]);
				File file2 = new File (args[1]);
				BufferedReader in1 = new BufferedReader(new FileReader(file1));
				BufferedReader in2 = new BufferedReader(new FileReader(file2));
				PrintWriter out = new PrintWriter(new FileWriter(new File(file1.getParent(), file1.getName()+"_"+file2.getName()+".merged")));
				String line1;
				String line2;

				while ((line1 = in1.readLine()) !=null  && (line2 = in2.readLine()) !=null){
					out.println(line1+"\t" +line2);
				}
				if (in1.readLine() != null || in2.readLine() !=null) System.out.println("\nWARNING: Files are of unequal length! Trunkated merged file.\n");
				in1.close();
				in2.close();
				out.close();
			}catch (Exception e){
				e.printStackTrace();
			}
		}
	}

	/**Assumes equal lenght files.*/
	public static void mergeColumns (File dir){
		File[] files = IO.extractOnlyFiles(dir);
		String[][] data = new String[files.length][];
		//load data
		for (int i=0; i< files.length; i++){
			data[i] = IO.loadFile(files[i]);
		}
		//print
		File outFile = new File (dir, "mergedColumns.txt");
		PrintWriter out;
		try {
			out = new PrintWriter( new FileWriter (outFile));

			//for each file print file name
			for (File f: files){
				out.print(Misc.removeExtension(f.getName()));
				out.print("\t");
			}
			out.println();

			//print data
			//for each row
			for (int r=0; r< data[0].length; r++){
				//for each file
				for (int i=0; i< files.length; i++){
					out.print(data[i][r]);
					out.print("\t");
				}
				out.println();
			}


			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}





	}
}
