package util.gen;
import java.io.File;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.temporal.ChronoUnit;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;
import java.util.regex.Matcher;

 
public class Delme {
	
	public static void main(String[] args) throws Exception {
		
		Double[] x = {1771022700000.0,
				1771023300000.0,
				1771023600000.0,
				1771024500000.0,
				1771025400000.0,
				1771025460000.0,
				1771025700000.0,
				1771026300000.0,
				1771026900000.0,
				1771027200000.0,
				1771027500000.0,
				1771028100000.0,
				1771029000000.0,
				1771029300000.0,
				1771029900000.0,
				1771030500000.0,
				1771030800000.0,
				1771031700000.0};
		
		Double[] y = {40.0,
				38.0,
				null,
				36.0,
				null,
				null,
				36.0,
				null,
				34.0,
				null,
				34.0,
				null,
				null,
				28.0,
				null,
				null,
				40.0,
				null};
		
		//check lengths are the same
		if (x.length!=y.length) throw new Exception("ERROR: array lengths differ.");
		
		//find first good index
		int goodIndex = -1;
		for (int i=0; i< y.length; i++) {
			if (y[i] != null) {
				goodIndex = i;
				break;
			}
		}
		IO.pl("First good index: "+goodIndex);
		
		//start with first index and start to walk
		for (int i=goodIndex+1; i< x.length; i++) {
			IO.pl("Looking at i "+i+"  x: "+x[i]+"  y: "+y[i]);
			//not null, just advance
			if (y[i] != null) {
				IO.pl("\tY not null advancing");
				goodIndex = i;
			}
			else {
				IO.pl("\tY is null");
				//find next good index
				int nextGoodIndex;
				for (int j=i+1; j< x.length; j++) {
					IO.pl("\t\tLooking for next non null j "+j+"  x: "+x[j]+"  y: "+y[j]);
					//not null
					if (y[j] != null) {
						IO.pl("\t\t\tNot null, interpolate");
						nextGoodIndex = j;
						//pull good values
						double x1 = x[goodIndex];
						double x2 = x[nextGoodIndex];
						double y1 = y[goodIndex];
						double y2 = y[nextGoodIndex];
						IO.pl("\t\t\t\t"+x1+" "+x2+" "+y1+" "+y2);
						for (int z=goodIndex+1; z<nextGoodIndex; z++) {
							
							double xTest = x[z];
							double yTest = Num.interpolateY(x1, y1, x2, y2, xTest);
							y[z] = yTest;
							IO.pl("\t\t\t\t\tInt xtest "+xTest+" yInt "+yTest+" setting in index "+z);
						}
						//reset
						i = nextGoodIndex;
						goodIndex = nextGoodIndex;
						break;
					}
					else IO.pl("\t\t\tStill null, advance!");
				}
			}
		}
		for (Double yf: y) {
			if (yf!=null)IO.pl(yf);
			else IO.pl(".");
		}
		
		
		
		/*
		double x1 = 123456;
		double y1 = 50;
		double x2 = 123500;
		double y2 = 55;
		double testX = 123477;
		
		//interpolateY(double x1, double y1, double x2, double y2, double testX)
		double testY = Num.interpolateY(x1, y1, x2, y2, testX);
		
		IO.pl(testY);*/
				
	}	
	
}	