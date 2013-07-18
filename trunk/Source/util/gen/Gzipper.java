package util.gen;

import java.io.*;
import java.util.ArrayList;
import java.util.zip.GZIPOutputStream;

/**Class for writing out compressed files. Be sure to close() this to clear the final buffer!*/
public class Gzipper {
	//fields
	private GZIPOutputStream out;
	private static final byte[] rtn = "\n".getBytes();
	private File gzipFile;
	
	//constructor
	public Gzipper (File gzipFile) throws FileNotFoundException, IOException{
		if (gzipFile.getName().endsWith(".gz") == false) gzipFile = new File (gzipFile+".gz");
		out = new GZIPOutputStream(new FileOutputStream(gzipFile));
		this.gzipFile = gzipFile;
	}
	
	public Gzipper(OutputStream gzipFile) throws IOException {
		out = new GZIPOutputStream(gzipFile);
	}
	
	/**Be sure to call this to clear the final buffer when done!*/
	public void close() throws IOException{
		out.close();
	}
	
	/**Adds a return onto the line*/
	public void println(String line) throws IOException{
		out.write(line.getBytes());
		out.write(rtn);
	}
	
	/**Adds a return onto the line*/
	public void println(Object line) throws IOException{
		out.write(line.toString().getBytes());
		out.write(rtn);
	}
	
	/**Adds a return onto each objects toString()*/
	public void println(Object[] obj) throws IOException{
		for (Object o: obj){
			out.write(o.toString().getBytes());
			out.write(rtn);
		}
	}
	
	public void print(String line) throws IOException{
		out.write(line.getBytes());
	}
	
	public void print(double line) throws IOException{
		out.write(new Double(line).toString().getBytes());
	}
	
	public void print(float line) throws IOException{
		out.write(new Float(line).toString().getBytes());
	}
	
	public void print(int line) throws IOException{
		out.write(new Integer(line).toString().getBytes());
	}
	
	public void println() throws IOException{
		out.write(rtn);
	}
	
	/**Adds a return onto each line*/
	public void println(ArrayList<String> lines) throws IOException{
		for (String line: lines){
			out.write(line.getBytes());
			out.write(rtn);
		}
	}
	
	/**Prints txt, txt.gz, or txt.zip*/
	public void print(File file) throws IOException{
		BufferedInputStream in = new BufferedInputStream(IO.fetchInputStream(file));
		//add file contents
		byte[] buf = new byte[1024];
		int i;
		while ((i = in.read(buf)) >= 0) out.write(buf, 0, i);
		//close in stream
		in.close();
	}

	public File getGzipFile() {
		return gzipFile;
	}
}
