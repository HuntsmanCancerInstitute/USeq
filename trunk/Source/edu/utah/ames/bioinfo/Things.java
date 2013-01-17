package edu.utah.ames.bioinfo;

/**
 * Call TomatoMasher as wrapper (sends email?), FileWalker descends and deletes
 *
 * @author darren.ames@hci.utah.edu
 */

/*

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

public class FileWalker {

    //fields
    Calendar cToday = Calendar.getInstance();
    Calendar fiveDaysAgo = Calendar.getInstance();
    static SimpleDateFormat DT_FORMAT = new SimpleDateFormat("MMM dd, yyyy");
    //cToday.set(Calendar.YEAR, (cToday.get(Calendar.YEAR)) - 1);

    //constructor
    FileWalker() {
        cToday.add(Calendar.DATE, -30);
        fiveDaysAgo.add(Calendar.DATE, -5);
    }

    //TODO: change so only prints filename, size, mod date, then path all
    //in same line. Doesn't print subdirs. First to print is date, then size,
    //then path
    public void walk(String path, int nbrOfTabs) {

        //find cmd.txt

        //fetch email

        //put in hashmap
        //use FileMap.get (pass in email), get FileInfos back
        //if null, instantiate and put in hash
        //walk(String email, String path, Map fileMap)
        //fileMap.put()


        //System.out.println("Walking " + path);
        File root = new File(path);
        File[] list = root.listFiles();

        nbrOfTabs++;
        if (list != null) {

            //TODO: sort method here? 
            //TODO: do we need permissions info here?
            for (File f : list) {
                boolean canRead = f.canRead();
                boolean canWrite = f.canWrite();
                boolean canExecute = f.canExecute();

                //make this independent lines
                char readChar = ' ', writeChar = ' ', execChar = ' ';

                //create calendar object for later time comparison
                Calendar fileModifiedDate = Calendar.getInstance();

                fileModifiedDate.setTimeInMillis(f.lastModified());

                if (f.canRead()) {
                    readChar = 'r';
                }
                if (f.canExecute()) {
                    execChar = 'x';
                }
                if (f.canWrite()) {
                    writeChar = 'w';
                }

                StringBuilder permissions = new StringBuilder(3);
                permissions.append(readChar).append(writeChar).append(execChar);
                permissions.append(" ");

                //open new hidden? file
                //TODO: overwrite monitorFile if exists
                File monitorFile = new File(root, "monitorFile");
                PrintWriter out;

                try {
                    out = new PrintWriter(new FileWriter(monitorFile, true));
                    out.println(monitorFile);


                    //get directory data and write to hidden file
                    if (f.isDirectory()) {
                        for (int i = 0; i < nbrOfTabs; i++) {
                            out.print("  ");
                        }
                        out.print("d");
                        out.print(permissions.toString());
                        out.print(" ");
                        out.print(f.length());
                        out.print(" ");
                        out.print(DT_FORMAT.format(fileModifiedDate.getTime()));
                        out.print(" ");
                        out.println(f.getPath());

                        //use getCanonicalPath (resolved path)
                        walk(f.getAbsolutePath(), nbrOfTabs);

                    } else {
                        if (fileModifiedDate.before(cToday)) {
                            continue;
                        }
                        if (fileModifiedDate.after(fiveDaysAgo) && f.length() > 1000) { //set this param 
                            // out.println("*****************deleting" + f.getPath());
                            //f.delete()
                        }

                        //send warning day 5, delete three days later

                        //delete file and replace with _deletedFile.txt


                        for (int i = 0; i < nbrOfTabs; i++) {
                            out.print("  ");
                        }
                        out.print(permissions.toString());
                        out.print(" ");
                        out.print(f.length());
                        out.print(" ");
                        out.print(DT_FORMAT.format(fileModifiedDate.getTime()));
                        out.print(" ");
                        out.println(f.getPath());
                        out.print(" ");

                    }

                    //close output file
                    out.close();

                } catch (Exception ex) {
                    Logger.getLogger(FileWalker.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
    }

    //**TODO: call EmailRegexTextFile and EmailAddress**
    public static void main(String[] args) throws IOException {

        //TODO: specify this from command line?
        String path = "/Users/darren/Desktop/testDir/";

        //create FileInfo object
        FileInfo f1 = new FileInfo();

        //create new FileWalker object
        FileWalker fw = new FileWalker();

        //specify starting directory for search as command line arg for main
        fw.walk(path, 0);

        //search for cmd.txt
        String file = null;
        File folder = new File(path);
        File[] listOfFiles = folder.listFiles();
        for (File checkedFile : listOfFiles) {

            if (checkedFile.isFile()) {
                file = checkedFile.getName();
                //System.out.println(file);

                if (file.equals(path + "cmd.txt")) {
                    f1.setFileName(file);
                    //System.out.println(f1.getFileName());
                }
            }
        }


        //extract email address from cmd.txt
        String e = null;
        EmailAddress email = new EmailAddress(e);

        //instantiate f1 using regular FileInfo constructor
        //f1 = new FileInfo(file, false, false, 0, 0);
        //System.out.println(file);

        //create hashmap with String email as key and FileInfo as value
        HashMap<String, String> hm = new HashMap<String, String>();

        //put elements to the map
        //TODO: set fileName and email
        hm.put(new String("email"), f1.getFileName());

        //get a set of the entries
        Set set = hm.entrySet();

        //get iterator
        Iterator i = set.iterator();

        //display the elements
        while (i.hasNext()) {
            Map.Entry me = (Map.Entry) i.next();
            //System.out.println(me.getKey());
            //System.out.println(me.getValue());
        }
    }
}
**/