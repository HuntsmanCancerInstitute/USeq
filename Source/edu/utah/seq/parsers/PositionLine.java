package edu.utah.seq.parsers;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

import edu.utah.seq.useq.apps.*;
import util.gen.IO;

public class PositionLine implements Comparable{

		private int position;
		private String line;
		
		public PositionLine(int position, String line){
			this.position = position;
			this.line = line;
		}
		
		/**Sorts by position.*/
		public int compareTo(Object other){
			PositionLine se = (PositionLine)other;
			if (position <se.position) return -1;
			if (position>se.position) return 1;
			return 0;
		}
		
		public static PositionLine[] load(File txtFile, int positionColumn, boolean sort){
			ArrayList<PositionLine> sss = new ArrayList<PositionLine>(10000);
			try {
				BufferedReader in = IO.fetchBufferedReader(txtFile);
				String line = null;
				String[] tokens = null;
				while ((line = in.readLine())!=null){
					tokens = Text2USeq.PATTERN_TAB.split(line);
					int position = Integer.parseInt(tokens[positionColumn]);
					sss.add(new PositionLine(position, line));
				}
			}
			catch (Exception e){
				e.printStackTrace();
				return null;
			} 
			PositionLine[] s = new PositionLine[sss.size()];
			sss.toArray(s);
			if (sort) Arrays.sort(s);
			return s;
		}

		
		/**Loads a binary file containing int,String 
		 * @return an UNSORTED array! or null if something bad happened.*/
		public static PositionLine[] loadBinary(File binaryFile){
			ArrayList<PositionLine> sss = new ArrayList<PositionLine>(10000);
			DataInputStream dis = null;
			try {
				dis = new DataInputStream(new BufferedInputStream(new FileInputStream(binaryFile)));
				while (true){
					int position = dis.readInt();
					//read line
					byte[] barray = new byte[dis.readInt()];
					dis.readFully(barray);
					String line = new String(barray);
					sss.add(new PositionLine(position, line));
				}
			} catch (EOFException eof){
				PositionLine[] s = new PositionLine[sss.size()];
				sss.toArray(s);
				return s;
			}
			catch (Exception e){
				e.printStackTrace();
				return null;
			} finally {
				if (dis != null) {
					try {
						dis.close();
					} catch (IOException ignored) {
						ignored.printStackTrace();
					}
				}
			}
		}

		public int getPosition() {
			return position;
		}

		public void setPosition(int position) {
			this.position = position;
		}

		public String getLine() {
			return line;
		}

		public void setLine(String line) {
			this.line = line;
		}
	}


