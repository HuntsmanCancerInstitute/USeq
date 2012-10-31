package trans.tpmap;
import java.io.*;
import java.util.*;

/**
 * Collapses a text version bpmap into a unique set of bpmap lines base on the oligo sequence.
 * */
public class CollapseTPMap {

	public static void main(String[] args) {
		File file = new File(args[0]);
		try{
			BufferedReader in = new BufferedReader(new FileReader (file));
			String[] tokens;
			String line;
			HashMap map = new HashMap(10000);
			while ( (line=in.readLine()) !=null){
				tokens = line.split("\\t");
				//check if it exists
				if (map.containsKey(tokens[0])){
					//get arraylist and add new
					ArrayList al = (ArrayList)map.get(tokens[0]);
					al.add(line);
				}
				else{
					ArrayList al = new ArrayList();
					al.add(line);
					map.put(tokens[0], al);
				}
			}
			in.close();
			
			PrintWriter out = new PrintWriter(new FileWriter(file.getCanonicalPath()+"c"));
			ArrayList lines = new ArrayList();
			lines.addAll(map.values());
			for (int i=0; i< lines.size(); i++){
				ArrayList al = (ArrayList)lines.get(i);
				for (int j=0; j<al.size(); j++){
					out.println(al.get(j));
				}
				out.println();
			}
			out.close();
			
		}catch (IOException e){
			e.printStackTrace();
		}

	}

}
