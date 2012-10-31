package edu.utah.bass;
import java.awt.Component;
import java.io.*;
import java.util.List;

import org.freehep.util.export.ExportFileType;


public class SequenceGraph {


	//fields
	private char[] sequenceTop;
	private double[] scoresTop;
	private char[] sequenceBottom;
	private double[] scoresBottom;
	private boolean savePNG = false;
	private File pngFile;
	private int width;
	private String matrixName;
	
	/*For testing*/
	public SequenceGraph(){
		sequenceTop =    "GAAUUAACCAAGGAAAAUAACAAGGACAGGGACCAGG".toCharArray();
		sequenceBottom = "CguAAuuGcuuCCuutuAuuGuuCCuGuCCCuGGuCC".toCharArray();
		
		scoresTop = new double[]{0.0,-1.0,-1.0,0.0,0.0,73.50612661354447,23.544396431897276,0.0,0.0,5.358792097779904,25.791687925929406,0.0,0.0,0.878771610792509,21.13183509393328,26.143738808601142,16.244837901374595,0.0,39.31562053645046,36.04382134793665,0.0,8.080348543012285,27.302727520142792,0.0,0.0,1.0411735504800348,0.0,9.572691267025657,0.0,0.0,0.0,0.894568655733588,0.0,0.0,-1.0,0.0,0.0};
		scoresBottom = scoresTop;
		
		savePNG = true;
		pngFile = new File ("/Users/davidnix/Desktop/test.png");

		System.setProperty("java.awt.headless","true");
		Component c = new SequenceGraphPanel (this);
		if (savePNG == false){
			File f = new File ("/Users/davidnix/Desktop/testC.eps");
			List<ExportFileType> efts = ExportFileType.getExportFileTypes("eps");
			ComponentWriter.exportComponent(f, c, efts.get(0));
		}
	}
	
	/**Set bottom fields null if you want to draw just the top*/
	public SequenceGraph(double[] scoresTop, double[] scoresBottom, char[] sequenceTop, char[] sequenceBottom, File pngFile, File epsFile, String matrixName){
		this.sequenceTop = sequenceTop;
		this.sequenceBottom = sequenceBottom;
		this.scoresTop = scoresTop;
		this.scoresBottom = scoresBottom;
		this.pngFile = pngFile;
		this.matrixName = matrixName;
		
		//save png file
		savePNG = true;
		System.setProperty("java.awt.headless","true");
		SequenceGraphPanel sgp = new SequenceGraphPanel (this);
		width = (int)sgp.getPanelWidth();
		
		//save eps
		savePNG = false;
		Component c = new SequenceGraphPanel (this);
		List<ExportFileType> efts = ExportFileType.getExportFileTypes("eps");
		ComponentWriter.exportComponent(epsFile, c, efts.get(0));
	}
	

	
	public static void main(String[] args) {
		new SequenceGraph();
	}

	public char[] getSequenceTop() {
		return sequenceTop;
	}
	public double[] getScoresTop() {
		return scoresTop;
	}
	public char[] getSequenceBottom() {
		return sequenceBottom;
	}
	public double[] getScoresBottom() {
		return scoresBottom;
	}
	public boolean isSavePNG() {
		return savePNG;
	}
	public void setSavePNG(boolean savePNG) {
		this.savePNG = savePNG;
	}
	public File getPngFile() {
		return pngFile;
	}
	public void setPngFile(File pngFile) {
		this.pngFile = pngFile;
	}

	public int getWidth() {
		return width;
	}

	public String getMatrixName() {
		return matrixName;
	}

}
