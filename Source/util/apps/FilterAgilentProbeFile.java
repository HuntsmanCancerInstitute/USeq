package util.apps;
import util.bio.annotation.*;
import java.io.*;
import java.util.*;

/**Assumes one chromosome.*/
public class FilterAgilentProbeFile {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//min hit
		int minHit = 22;
		
		//load file
		File probeFile = new File (args[0]);
		Coordinate[] probes = Coordinate.parseFile(probeFile, 0, 1);
		
		//lengthen 3' stop to make to 60mer, 3' stop is buried in matrix
		int min = probes[0].getStart();
		int max = probes[0].getStop();
		
		for (int i=0; i< probes.length; i++){
			//reset length, if needed, to 59bp
			int diff = 59 - (probes[i].getStop()-probes[i].getStart());
			if (diff !=0) probes[i].setStop(probes[i].getStop()+diff);

			//subtract 20bp from 3' so only presenting region that actually can hybridize
			probes[i].setStop(probes[i].getStop() - 20);
			
			//set min max
			if (probes[i].getStart() < min) min = probes[i].getStart();
			if (probes[i].getStop() > max) max = probes[i].getStop();
		}
		
		//put into zero coordinates by subtracting min from both start and stop
		for (int i=0; i< probes.length; i++){
			probes[i].setStart(probes[i].getStart() - min);
			probes[i].setStop(probes[i].getStop() - min);
		}
		
		//sort
		Arrays.sort(probes);
		
		//make boolean to hold tiled regions, default is false;
		boolean[] tiled = new boolean[max - min + 1];
		
		//identify uniques
		ArrayList<Coordinate> unique = new ArrayList <Coordinate>();
		ArrayList<Coordinate> intersecting = new ArrayList <Coordinate>();
		
		//assumes sorted
		for (int i=0; i< probes.length; i++){
			boolean uni = true;
			for (int j=i+1; j< probes.length; j++){
				if (probes[i].intersects(probes[j])) {
					intersecting.add(probes[i]);
					uni = false;
					break;
				}
			}
			if (uni) {
				unique.add(probes[i]);
				//load into boolean
				int start = probes[i].getStart();
				int stop = probes[i].getStop() +1;
				for (int x = start; x < stop; x++) tiled[x] = true;
			}
		}
		
		//for each intersecting
		ArrayList<Coordinate> toAdd = new ArrayList <Coordinate>();
		int num = intersecting.size();
		for (int i=0; i< num; i++){
			Coordinate test = intersecting.get(i);
			//count coverage
			int start = test.getStart();
			int stop = test.getStop() +1;
			int numNotHit = 0;
			for (int x = start; x < stop; x++) {
				if (tiled[x]==false) numNotHit++;
			}
			//good coverage?
			if (numNotHit > minHit) {
				//add it
				toAdd.add(test);
				for (int x = start; x < stop; x++) tiled[x]=true;
			}
		}
		
		
		System.out.println("Num uni 1st "+unique.size());
		System.out.println("Num int 1st "+intersecting.size());
		System.out.println("Num 2Add 1st "+toAdd.size());
		
		
		System.out.println("Num total "+(toAdd.size()+unique.size()));
//System.exit(0);
		//combine
		unique.addAll(toAdd);
		//unique.addAll(intersecting);
		probes = new Coordinate[unique.size()];
		unique.toArray(probes);
		
		//add  21 to 3' stop
		for (int i=0; i< probes.length; i++) {
			probes[i].setStart(min + probes[i].getStart());
			probes[i].setStop( min + probes[i].getStop()+ 21);
		}
		Arrays.sort(probes);
		
		//print
		for (int i=0; i< probes.length; i++) System.out.println(probes[i]);
	}

}
