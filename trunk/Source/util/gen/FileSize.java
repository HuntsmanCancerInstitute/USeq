package util.gen;

import java.io.File;

public class FileSize implements Comparable{
		private String name;
		private long size;
		
		public FileSize (File f){
			name = f.getName();
			size = f.length();
		}
		public int compareTo (Object obj){
			FileSize other = (FileSize)obj;
			if (other.size< size) return 1;
			if (other.size> size) return -1;
			return 0;
		}
		public String getName() {
			return name;
		}
		public void setName(String name) {
			this.name = name;
		}
		public long getSize() {
			return size;
		}
		public void setSize(long size) {
			this.size = size;
		}
	
}
