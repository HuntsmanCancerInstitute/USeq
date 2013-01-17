package edu.utah.ames.bioinfo;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.zip.*;


public class Stuff {

	 //create hashmap with String email as key and FileInfo as value
    /* HashMap<String, String> hm = new HashMap<String, String>();

     //put elements to the map
     hm.put(new String("email"), f1.getCmdFile());

     //get a set of the entries
     Set set = hm.entrySet();

     //get iterator
     Iterator i = set.iterator();

     //display the elements
     while (i.hasNext()) {
         Map.Entry me = (Map.Entry) i.next();
         //System.out.println(me.getKey());
         //System.out.println(me.getValue());
           
          
     }*/
	
	/**Extracts the full path file names of all the files and directories in a given directory. If a file is given it is
	* returned as the File[0].
	* Skips files starting with a '.'*/
	public static File[] extractFiles(File directory) {
		File[] files = null;
		String[] fileNames;
		if (directory.isDirectory()) {
			fileNames = directory.list();
			int num = fileNames.length;
			ArrayList<File> al = new ArrayList<File>();
			try {
				String path = directory.getCanonicalPath();
				for (int i = 0; i < num; i++) {
					if (fileNames[i].startsWith(".") == false)
						al.add(new File(path, fileNames[i]));
				}
				// convert arraylist to file[]
				if (al.size() != 0) {
					files = new File[al.size()];
					al.toArray(files);
				}
			} catch (IOException e) {
				System.out.println("Problem extractFiles() " + directory);
				e.printStackTrace();
				return null;
			}
		}
		if (files == null) {
			files = new File[1];
			files[0] = directory;
		}
		Arrays.sort(files);
		return files;
	}
	
	
	public static boolean saveObject(File file, Object ob) {
		try {
			ObjectOutputStream out =
				new ObjectOutputStream(new FileOutputStream(file));
			out.writeObject(ob);
			out.close();
			return true;
		} catch (Exception e) {
			System.out.println("There appears to be a problem with saving this file: "+ file);
			e.printStackTrace();
		}
		return false;
	}
	
	/**Fetches an Object stored as a serialized file.
	 * Can be zip/gz compressed too.*/
	public static Object fetchObject(File file) {
		Object a = null;
		try {
			ObjectInputStream in;
			if (file.getName().endsWith(".zip")){
				ZipFile zf = new ZipFile(file);
				ZipEntry ze = zf.entries().nextElement();
				in = new ObjectInputStream( zf.getInputStream(ze));
			}
			else if (file.getName().endsWith(".gz")) {
				in = new ObjectInputStream(new GZIPInputStream(new FileInputStream(file)));
			}
			else in = new ObjectInputStream(new FileInputStream(file));
			a = in.readObject();
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Problem fetchObject() "+file);
		}
		return a;
	}
}
