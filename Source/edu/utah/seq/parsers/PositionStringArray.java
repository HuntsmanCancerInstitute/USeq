package edu.utah.seq.parsers;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

import edu.utah.seq.useq.data.Region;
import edu.utah.seq.useq.data.RegionScoreText;


import util.gen.IO;

public class PositionStringArray implements Comparable{

		private int position;
		private String[] line;
		
		public PositionStringArray(int position, String[] line){
			this.position = position;
			this.line = line;
		}
		
		/**Sorts by position.*/
		public int compareTo(Object other){
			PositionStringArray se = (PositionStringArray)other;
			if (position <se.position) return -1;
			if (position>se.position) return 1;
			return 0;
		}

		public int getPosition() {
			return position;
		}

		public String[] getLine() {
			return line;
		}
	}


