package edu.utah.ames.bioinfo;
import java.io.*;

public class MyFileNameFilter implements FilenameFilter {

	static String PATH = "/Users/darren/Desktop/testDir/";
	String cmdFile = "cmd.txt";
	@Override
	public boolean accept(File arg0, String arg1) {
		// TODO Auto-generated method stub
		FileInfo fi = new FileInfo();
		
		boolean result = false;
		if (arg1.equals(cmdFile));
			result = true;
			fi.setCmdFile(arg0);
			//System.out.println(arg1);
			System.out.println(fi.getCmdFile());
		return result;
		
	}

	public static void main(String[] args) {
		File f = new File(PATH);
		String[] files = null;
		if (f.isDirectory()) {
			files = f.list(new MyFileNameFilter());
		}

	}
}
