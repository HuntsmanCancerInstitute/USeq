package trans.graphics;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import org.jfree.chart.*;
import org.jfree.data.xy.*;
import org.jfree.chart.plot.*;
import org.jfree.chart.renderer.xy.*;

import trans.misc.*;
import trans.qc.ChartFrame;
import util.bio.seq.*;
import util.gen.*;

public class PlotFloatArrays {
	
	public static void main(String[] args){
		if (args.length == 0){
			Misc.printExit("\nEnter the full path directory text for a set of float[] files to graph. Will invert float[] whose file text begins with 'as....'\n");
		}
		
		//load float[][]s
		
		File directory = new File(args[0]);
		File[] floatFiles = IO.extractFiles(directory);

		float[][] floats = new float[floatFiles.length][];
		for (int i=0; i< floatFiles.length; i++){
			floats[i] = (float[])IO.fetchObject(floatFiles[i]);
			if (floatFiles[i].getName().startsWith("as")) {
				System.out.println("\tInverting "+floatFiles[i].getName());
				Num.invertArray(floats[i]);
			}
		}
		
		//print every 10th base
		int halfLength = (int)Math.round(((double)floats[0].length-1) / 2);
		for (int a=0; a< floats.length; a++){
			//for every 10th score
			System.out.println("\n\n"+floatFiles[a].getName());
			for (int b=0; b<floats[0].length; b+=5){
				//print score
				if (floats[a][b] !=0){
					int index = b- halfLength;
					System.out.println(index+"\t"+floats[a][b]);
				}
			}
		}
		
		//plot
		//make series
		XYSeries[] series = new XYSeries[floatFiles.length];

		for (int i=0; i< floatFiles.length; i++){
			series[i] = new XYSeries( Misc.removeExtension(floatFiles[i].getName()));
		}
		
		//load series
		//for each float
		for (int a=0; a< floats.length; a++){
			//for each score
			for (int b=0; b<floats[0].length; b++){
				//print score
				if (floats[a][b] !=0){
					int index = b- halfLength;
					series[a].add(index, floats[a][b]);
				}
			}
		}
		
		//make collection
		XYSeriesCollection dataset = new XYSeriesCollection(); 
		for (int i=0; i< series.length; i++) dataset.addSeries(series[i]); 
		
		JFreeChart chart = ChartFactory.createXYLineChart( 
				null, // chart title 
				null, // x axis label 
				"Win Score", // y axis label 
				dataset, // data 
				PlotOrientation.VERTICAL, 
				true, // include legend 
				true, // tooltips 
				false // urls 
		); 
		
		XYPlot plot = (XYPlot) chart.getPlot(); 
		XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer(); 
		for (int i=0; i< series.length; i++){
			renderer.setSeriesShapesVisible(i, false);
			renderer.setSeriesLinesVisible(i, true);
			renderer.setSeriesShapesFilled(i, false);
		}
		
		try {
			File png = new File(floatFiles[0].getParentFile(), "mergedPlot.png");
			System.out.println("\tSaving aggregate plot "+png);
			ChartUtilities.saveChartAsPNG(png,chart, 1440,960);
		} catch (Exception e){
			e.printStackTrace();
		}

			JFreeChart[] charts = {chart};
			new ChartFrame(charts, "Merged Aggregate Plot", 1,1);
		
		
	}

}
