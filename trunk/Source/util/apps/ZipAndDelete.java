package util.apps;
import java.io.*;
import util.gen.*;
import java.util.*;

public class ZipAndDelete {

	public static void main(String[] args) {
		if (args.length ==0){
			Misc.printExit("\nProvide one or more directories to recurse through and zip files within.  " +
					"After successful compression, files are deleted. Compressed files (.bz2, .gz, .rar, " +
					".bz, .zip, .jpg, .png) will be skipped.\n");
		}
		System.out.println("\nArchiving...");
		for (int x=0; x< args.length; x++){
			File dir = new File(args[x]);
			if (dir.exists() == false || dir.isDirectory() == false || dir.canRead() == false ) Misc.printExit("Error: aborting, is this a directory? exits? read/ write permissions? -> "+args[x]);
			//fetch files recursively
			File[] files = IO.fetchFilesRecursively(dir);
			//filter to remove compressed files
			ArrayList<File> filteredFilesAL = new ArrayList<File>();
			for (int i =0; i< files.length; i++){
				String name = files[i].getName().toLowerCase();
				if (name.endsWith(".zip") || name.endsWith(".bz2") || name.endsWith(".gz") || name.endsWith(".bz") || name.endsWith(".rar") || name.endsWith(".jpg") || name.endsWith(".png")) continue;
				filteredFilesAL.add(files[i]);
			}
			files = new File[filteredFilesAL.size()];
			filteredFilesAL.toArray(files);
			//zip and delete
			System.out.println("\t"+files.length+ "\tFiles in "+dir);
			for (int i =0; i< files.length; i++){
				if (IO.zipAndDelete(files[i]) == false) Misc.printExit("\nZip failed, aborting, see file -> "+files[i]);
			}
		}
		System.out.println("\nDone!\n");
	}

}
