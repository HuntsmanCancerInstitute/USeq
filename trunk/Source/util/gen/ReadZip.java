package util.gen;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

public class ReadZip {
  public static void main(String args[]) {
    try {
      ZipFile zf = new ZipFile("/Users/dnix/Desktop/chrMT.gr.zip");
      
        ZipEntry ze = (ZipEntry) zf.entries().nextElement();
          
            BufferedReader br = new BufferedReader(
                new InputStreamReader(zf.getInputStream(ze)));
            String line;
            while ((line = br.readLine()) != null) {
              System.out.println(line);
            }
            br.close();
          
        
      
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
}