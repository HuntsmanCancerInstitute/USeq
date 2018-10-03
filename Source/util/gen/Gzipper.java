package util.gen;

import java.io.*;
import java.util.ArrayList;
import java.util.zip.GZIPOutputStream;

/**Class for writing out compressed files. Be sure to close() this to clear the final buffer!*/
public class Gzipper {
	//fields
	private GZIPOutputStream out;
	private static final byte[] rtn = "\n".getBytes();
	private static final byte[] rtnComma = ",\n".getBytes();
	private File gzipFile;
	
	//constructor
	public Gzipper (File gzipFile) throws FileNotFoundException, IOException{
		if (gzipFile.getName().endsWith(".gz") == false) this.gzipFile = new File (gzipFile+".gz");
		else this.gzipFile = gzipFile;
		out = new GZIPOutputStream(new FileOutputStream(this.gzipFile));
	}
	
	public Gzipper(OutputStream gzipFile) throws IOException {
		out = new GZIPOutputStream(gzipFile);
	}
	
	/**Be sure to call this to clear the final buffer when done!*/
	public void close() throws IOException{
		out.close();
	}
	
	/**Flush the stream
	 * @throws IOException */
	public void flush() throws IOException{
		out.flush();
	}
	
	/**Attempts a close, no exceptions thrown or warnings.*/
	public void closeNoException(){
		try {
			out.close();
		} catch (IOException e) {}
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
	
	/**Adds a comma and return onto the line*/
	public void printcln(Object line) throws IOException{
		out.write(line.toString().getBytes());
		out.write(rtnComma);
	}
	
	/**Adds a return onto each objects toString()*/
	public void println(Object[] obj) throws IOException{
		for (Object o: obj){
			out.write(o.toString().getBytes());
			out.write(rtn);
		}
	}
	
	/**Adds a delimiter to each element in the array then a final return*/
	public void println(Object[] obj, String delimiter) throws IOException{
		byte[] db = delimiter.getBytes();
		out.write(obj[0].toString().getBytes());
		for (int i=1; i< obj.length; i++){
			out.write(db);
			out.write(obj[i].toString().getBytes());
		}
		out.write(rtn);
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
	
	//a bunch of methods for writing jsonified key value pairs
	public void printJson(String key, double[] values, boolean addComma) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append("\"");
		sb.append(key);
		sb.append("\"");
		sb.append(": [");
		sb.append(values[0]);
		for (int i=1; i< values.length; i++){
			sb.append(", ");
			sb.append(values[i]);
		}
		sb.append("]");
		out.write(sb.toString().getBytes());
		if (addComma) out.write(rtnComma);
		else out.write(rtn);
	}
	/**Calls toString on each object in the array, surrounds with "". No tabs or returns!*/
	public void printJson(String key, Object[] values, boolean addComma) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append("\"");
		sb.append(key);
		sb.append("\"");
		sb.append(": [\"");
		sb.append(values[0].toString());
		for (int i=1; i< values.length; i++){
			sb.append("\",\"");
			sb.append(values[i].toString());
		}
		sb.append("\"]");
		out.write(sb.toString().getBytes());
		if (addComma) out.write(rtnComma);
		else out.write(rtn);
	}

	public void printJson(String key, double value, boolean addComma) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append("\"");
		sb.append(key);
		sb.append("\"");
		sb.append(": ");
		sb.append(value);
		out.write(sb.toString().getBytes());
		if (addComma) out.write(rtnComma);
		else out.write(rtn);
	}
	public void printJson(String key, boolean value, boolean addComma) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append("\"");
		sb.append(key);
		sb.append("\"");
		sb.append(": ");
		sb.append(value);
		out.write(sb.toString().getBytes());
		if (addComma) out.write(rtnComma);
		else out.write(rtn);
	}
	public void printJson(String key, long value, boolean addComma) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append("\"");
		sb.append(key);
		sb.append("\"");
		sb.append(": ");
		sb.append(value);
		out.write(sb.toString().getBytes());
		if (addComma) out.write(rtnComma);
		else out.write(rtn);
	}
	public void printJson(String key, String value, boolean addComma) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append("\"");
		sb.append(key);
		sb.append("\"");
		sb.append(": ");
		sb.append("\"");
		sb.append(value);
		sb.append("\"");
		out.write(sb.toString().getBytes());
		if (addComma) out.write(rtnComma);
		else out.write(rtn);
	}
	public void printJson(String key, int value, boolean addComma) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append("\"");
		sb.append(key);
		sb.append("\"");
		sb.append(": ");
		sb.append(value);
		out.write(sb.toString().getBytes());
		if (addComma) out.write(rtnComma);
		else out.write(rtn);
	}
	public void printJson(String key, ArrayList<String> values, boolean addComma) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append("\"");
		sb.append(key);
		sb.append("\"");
		sb.append(": [\"");
		if (values.size() !=0){
			sb.append(values.get(0).toString());
			int len = values.size();
			for (int i=1; i< len; i++){
				sb.append("\", \"");
				sb.append(values.get(i).toString());
			}
		}
		sb.append("\"]");
		out.write(sb.toString().getBytes());
		if (addComma) out.write(rtnComma);
		else out.write(rtn);
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
	
	public static void main(String[] args) throws FileNotFoundException, IOException{
		String[] x = {"Go", "Fish", "Larry"};
		Gzipper g = new Gzipper(new File("/Users/u0028003/Downloads/delme.txt.gz"));
		g.printJson("ziper", x, false);
		g.close();
	}

	

}
